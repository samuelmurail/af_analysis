"""Flask GUI for af_analysis with interactive Plotly + Mol* panels."""

from __future__ import annotations

import json
import queue
import traceback
from pathlib import Path
from threading import Lock, Thread
from typing import Any

import numpy as np
import scipy.cluster.hierarchy as _sch
import scipy.spatial.distance as _ssd
from flask import Flask, Response, jsonify, render_template, request, send_from_directory

import af_analysis

APP_DIR = Path(__file__).resolve().parent


def _normalize_for_json(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _normalize_for_json(val) for key, val in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_normalize_for_json(item) for item in value]
    if hasattr(value, "tolist"):
        try:
            return _normalize_for_json(value.tolist())
        except Exception:
            return str(value)
    if isinstance(value, float):
        # NaN and Inf are not valid JSON — convert to null so JSON.parse succeeds.
        import math

        if math.isnan(value) or math.isinf(value):
            return None
        return value
    if value is None or isinstance(value, (str, int, bool)):
        return value
    return str(value)


class AppState:
    def __init__(self) -> None:
        self.lock = Lock()
        self.af_data: af_analysis.Data | None = None
        self.last_error = ""
        self.cluster_universes: dict = {}  # {query: (universe, [pdb_files])}


class _TqdmStream:
    """Drop-in tqdm replacement that pushes progress events to a queue.

    It mirrors the subset of the tqdm API used by analysis.py so it can be
    injected as a module-level monkey-patch without changing library code.
    """

    def __init__(self, iterable=None, *, total=None, disable=False, desc=None, **_kw):
        self._iter = iter(iterable) if iterable is not None else None
        self._total = (
            total
            if total is not None
            else (len(iterable) if hasattr(iterable, "__len__") else None)
        )
        self._n = 0
        self._disable = disable
        self._desc = desc or ""
        self._q: queue.Queue | None = _PROGRESS_QUEUE

    # ── context manager ──────────────────────────────────────────────────────
    def __enter__(self):
        return self

    def __exit__(self, *_):
        if not self._disable:
            self._push(done=True)

    # ── iteration ────────────────────────────────────────────────────────────
    def __iter__(self):
        return self

    def __next__(self):
        try:
            value = next(self._iter)
        except StopIteration:
            if not self._disable:
                self._push(done=True)
            raise
        self._n += 1
        if not self._disable:
            self._push()
        return value

    # ── manual update (not used by analysis.py but keeps API complete) ────────
    def update(self, n=1):
        self._n += n
        if not self._disable:
            self._push()

    def close(self):
        pass

    # ── internal ─────────────────────────────────────────────────────────────
    def _push(self, done=False):
        if self._q is None:
            return
        try:
            self._q.put_nowait(
                {
                    "desc": self._desc,
                    "n": self._n,
                    "total": self._total,
                    "done": done,
                }
            )
        except queue.Full:
            pass  # never block the computation thread


# One queue per server process; computations stream into it.
# Bounded so a slow client can't grow memory unboundedly.
_PROGRESS_QUEUE: queue.Queue = queue.Queue(maxsize=512)


STATE = AppState()

app = Flask(
    __name__,
    template_folder=str(APP_DIR / "templates"),
    static_folder=str(APP_DIR / "static"),
)
# Prevent NaN/Infinity literals in JSON responses (not valid JSON, breaks browser JSON.parse).
app.config["JSON_SORT_KEYS"] = False


def _require_data() -> af_analysis.Data:
    if STATE.af_data is None:
        raise RuntimeError("Dataset is not loaded")
    return STATE.af_data


def _preferred_table_columns(df) -> list[str]:
    # Always keep these identifier columns when present.
    id_cols = ["query", "model", "seed"]
    # Add every column whose dtype is numeric (int or float).
    numeric_cols = list(df.select_dtypes(include="number").columns)
    # Preserve id_cols first (in order), then numeric cols not already included.
    seen = set(id_cols)
    columns = [c for c in id_cols if c in df.columns]
    columns += [c for c in numeric_cols if c not in seen]
    return columns


@app.get("/")
def index():
    return render_template("flask_index.html")


@app.get("/favicon.ico")
def favicon():
    return send_from_directory(str(APP_DIR / "static" / "img"), "favicon.ico", mimetype="image/x-icon")


@app.post("/api/load")
def api_load_dataset():
    import sys
    import tqdm as _tqdm_mod

    payload = request.get_json(silent=True) or {}
    directory = str(payload.get("directory", "")).strip()
    format_name = str(payload.get("format", "auto")).strip() or "auto"

    if not directory:
        return jsonify({"error": "Results directory is required"}), 400

    # Drain stale progress events from any previous operation.
    while not _PROGRESS_QUEUE.empty():
        try:
            _PROGRESS_QUEUE.get_nowait()
        except queue.Empty:
            break

    # Patch every af_analysis.* module that has its own `tqdm` binding
    # (each does `from tqdm.auto import tqdm`, creating an independent reference).
    _patched: dict = {}
    for _mod_name, _mod in list(sys.modules.items()):
        if _mod_name.startswith("af_analysis") and hasattr(_mod, "tqdm"):
            _patched[_mod_name] = _mod.tqdm
            _mod.tqdm = _TqdmStream
    _orig_auto = _tqdm_mod.auto.tqdm
    _tqdm_mod.auto.tqdm = _TqdmStream

    try:
        format_arg = None if format_name == "auto" else format_name
        data = af_analysis.Data(directory=directory, format=format_arg, verbose=True)
        with STATE.lock:
            STATE.af_data = data
            STATE.last_error = ""
    except Exception:
        with STATE.lock:
            STATE.last_error = traceback.format_exc()
        _PROGRESS_QUEUE.put_nowait(
            {"desc": "__end__", "n": 0, "total": 0, "done": True}
        )
        return (
            jsonify({"error": "Failed to load dataset", "details": STATE.last_error}),
            500,
        )
    finally:
        for _mod_name, _orig in _patched.items():
            if _mod_name in sys.modules:
                sys.modules[_mod_name].tqdm = _orig
        _tqdm_mod.auto.tqdm = _orig_auto

    _PROGRESS_QUEUE.put_nowait({"desc": "__end__", "n": 0, "total": 0, "done": True})
    return jsonify({"ok": True, "rows": int(len(data.df))})


@app.get("/api/browse")
def api_browse():
    import os

    raw = request.args.get("path", "~")
    try:
        p = Path(raw).expanduser().resolve()
        if not p.is_dir():
            p = p.parent
        entries = sorted(
            [e.name for e in p.iterdir() if e.is_dir() and not e.name.startswith(".")]
        )
        return jsonify(
            {
                "path": str(p),
                "parent": str(p.parent) if p != p.parent else None,
                "dirs": entries,
            }
        )
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.get("/api/table")
def api_table():
    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    df = data.df.reset_index(drop=True)
    columns = _preferred_table_columns(df)
    rows = []
    for idx, row in df.iterrows():
        item = {"row": int(idx)}
        for column in columns:
            item[column] = _normalize_for_json(row[column])
        rows.append(item)

    has_lis = "LIS" in df.columns and bool(
        df["LIS"].apply(lambda v: isinstance(v, (list, np.ndarray))).any()
    )
    has_lia = "LIA" in df.columns and bool(
        df["LIA"].apply(lambda v: isinstance(v, (list, np.ndarray))).any()
    )
    has_iptm_d0_matrix = "ipTM_d0_matrix" in df.columns and bool(
        df["ipTM_d0_matrix"].apply(lambda v: isinstance(v, (list, np.ndarray))).any()
    )
    has_ipsae_matrix = "ipSAE_matrix" in df.columns and bool(
        df["ipSAE_matrix"].apply(lambda v: isinstance(v, (list, np.ndarray))).any()
    )
    return jsonify(
        {
            "columns": ["row", *columns],
            "rows": rows,
            "total": int(len(df)),
            "has_lis": has_lis,
            "has_lia": has_lia,
            "has_iptm_d0_matrix": has_iptm_d0_matrix,
            "has_ipsae_matrix": has_ipsae_matrix,
        }
    )


@app.get("/api/plddt")
def api_plddt():
    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    index = request.args.get("index", default=0, type=int)
    if index < 0 or index >= len(data.df):
        return jsonify({"error": "Model index out of range"}), 400

    plddt_array = data.get_plddt(index)
    if plddt_array is None:
        return jsonify({"error": "No pLDDT data for this model"}), 404

    plddt = np.asarray(plddt_array).reshape(-1).tolist()
    residues = list(range(1, len(plddt) + 1))

    row = data.df.iloc[index]
    query = row.get("query")
    chain_ids = [str(c) for c in data.chains.get(query, [])]
    chain_lengths = [int(x) for x in data.chain_length.get(query, [])]

    boundaries = []
    running = 0
    for length in chain_lengths[:-1]:
        running += int(length)
        boundaries.append(running)

    chain_types = [str(t) for t in data.chain_type.get(query, [])]

    return jsonify(
        {
            "residues": residues,
            "plddt": plddt,
            "chain_ids": chain_ids,
            "chain_lengths": chain_lengths,
            "chain_boundaries": boundaries,
            "chain_types": chain_types,
        }
    )


@app.get("/api/pae")
def api_pae():
    from af_analysis.analysis import get_pae

    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    index = request.args.get("index", default=0, type=int)
    if index < 0 or index >= len(data.df):
        return jsonify({"error": "Model index out of range"}), 400

    row = data.df.iloc[index]
    data_file = row.get("data_file")
    if not data_file or (
        hasattr(data_file, "__class__") and str(data_file) in ("", "nan")
    ):
        return jsonify({"error": "No data file for this model"}), 404

    pae_matrix = get_pae(str(data_file))
    if pae_matrix is None:
        return jsonify({"error": "No PAE data for this model"}), 404

    pae_matrix = np.asarray(pae_matrix)
    n = pae_matrix.shape[0]
    residues = list(range(1, n + 1))

    query = row.get("query")
    chain_ids = [str(c) for c in data.chains.get(query, [])]
    chain_lengths = [int(x) for x in data.chain_length.get(query, [])]

    boundaries = []
    running = 0
    for length in chain_lengths[:-1]:
        running += int(length)
        boundaries.append(running)

    return jsonify(
        {
            "pae": pae_matrix.tolist(),
            "residues": residues,
            "chain_ids": chain_ids,
            "chain_lengths": chain_lengths,
            "chain_boundaries": boundaries,
        }
    )


@app.get("/api/lis")
def api_lis():
    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    index = request.args.get("index", default=0, type=int)
    if index < 0 or index >= len(data.df):
        return jsonify({"error": "Model index out of range"}), 400

    row = data.df.iloc[index]
    lis = row.get("LIS")
    if not isinstance(lis, (list, np.ndarray)):
        return jsonify({"error": "No LIS data for this model"}), 404

    query = row.get("query")
    chain_ids = [str(c) for c in data.chains.get(query, [])]
    return jsonify(
        {"lis": _normalize_for_json(lis), "chain_ids": chain_ids, "label": "LIS"}
    )


@app.get("/api/lia")
def api_lia():
    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    index = request.args.get("index", default=0, type=int)
    if index < 0 or index >= len(data.df):
        return jsonify({"error": "Model index out of range"}), 400

    row = data.df.iloc[index]
    lia = row.get("LIA")
    if not isinstance(lia, (list, np.ndarray)):
        return jsonify({"error": "No LIA data for this model"}), 404

    query = row.get("query")
    chain_ids = [str(c) for c in data.chains.get(query, [])]
    return jsonify(
        {"lis": _normalize_for_json(lia), "chain_ids": chain_ids, "label": "cLIS (LIA)"}
    )


@app.get("/api/iptm_d0")
def api_iptm_d0():
    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    index = request.args.get("index", default=0, type=int)
    if index < 0 or index >= len(data.df):
        return jsonify({"error": "Model index out of range"}), 400

    row = data.df.iloc[index]
    matrix = row.get("ipTM_d0_matrix")
    if not isinstance(matrix, (list, np.ndarray)):
        return jsonify({"error": "No ipTM_d0 matrix data for this model"}), 404

    query = row.get("query")
    chain_ids = [str(c) for c in data.chains.get(query, [])]
    return jsonify(
        {"lis": _normalize_for_json(matrix), "chain_ids": chain_ids, "label": "ipTM d0"}
    )


@app.get("/api/ipsae")
def api_ipsae():
    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    index = request.args.get("index", default=0, type=int)
    if index < 0 or index >= len(data.df):
        return jsonify({"error": "Model index out of range"}), 400

    row = data.df.iloc[index]
    matrix = row.get("ipSAE_matrix")
    if not isinstance(matrix, (list, np.ndarray)):
        return jsonify({"error": "No ipSAE matrix data for this model"}), 404

    query = row.get("query")
    chain_ids = [str(c) for c in data.chains.get(query, [])]
    return jsonify(
        {"lis": _normalize_for_json(matrix), "chain_ids": chain_ids, "label": "ipSAE"}
    )


@app.get("/api/structure")
def api_structure():
    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    index = request.args.get("index", default=0, type=int)
    if index < 0 or index >= len(data.df):
        return jsonify({"error": "Model index out of range"}), 400

    row = data.df.iloc[index]
    pdb_path = row.get("pdb")
    if not pdb_path:
        return jsonify({"error": "No structure file for this model"}), 404

    path_obj = Path(str(pdb_path))
    if not path_obj.is_file():
        return jsonify({"error": f"Structure file not found: {path_obj}"}), 404

    ext = path_obj.suffix.lower()
    if ext in {".cif", ".mmcif"}:
        structure_format = "mmcif"
    elif ext == ".pdb":
        structure_format = "pdb"
    else:
        structure_format = "mmcif"

    try:
        structure_text = path_obj.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return (
            jsonify(
                {
                    "error": "Failed to read structure file",
                    "details": traceback.format_exc(),
                }
            ),
            500,
        )

    return jsonify(
        {"structure_text": structure_text, "structure_format": structure_format}
    )


@app.get("/api/health")
def api_health():
    loaded = STATE.af_data is not None
    rows = int(len(STATE.af_data.df)) if loaded else 0
    directory = (
        str(STATE.af_data.dir) if loaded and getattr(STATE.af_data, "dir", None) else ""
    )
    return jsonify({"loaded": loaded, "rows": rows, "directory": directory})


@app.get("/api/progress/stream")
def api_progress_stream():
    """Server-Sent Events stream of tqdm progress from the current computation."""

    def generate():
        while True:
            try:
                event = _PROGRESS_QUEUE.get(timeout=60)  # wait up to 60 s
            except queue.Empty:
                # Send a keep-alive comment so the connection stays open.
                yield ": keep-alive\n\n"
                continue
            payload = json.dumps(event)
            yield f"data: {payload}\n\n"
            if event.get("done") and event.get("desc") == "__end__":
                return

    return Response(
        generate(),
        mimetype="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


@app.post("/api/compute")
def api_compute():
    """Run one of the supported scoring functions on the loaded dataset.

    Body JSON: { "score": "pdockq2" | "LIS_LIA" | "iptm_d0" }
    Returns:   { "ok": true, "columns_added": [...] }
    """
    import af_analysis.analysis as _ana_mod
    import af_analysis.data as _data_mod
    import tqdm as _tqdm_mod

    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    payload = request.get_json(silent=True) or {}
    score = str(payload.get("score", "")).strip()
    pae_cutoff = float(payload.get("pae_cutoff", 12.0))
    dist_cutoff = float(payload.get("dist_cutoff", 8.0))

    if score not in {"pdockq2", "LIS_LIA", "iptm_d0"}:
        return jsonify({"error": f"Unknown score '{score}'"}), 400

    # Drain any stale events from a previous run.
    while not _PROGRESS_QUEUE.empty():
        try:
            _PROGRESS_QUEUE.get_nowait()
        except queue.Empty:
            break

    # Monkey-patch tqdm in both the analysis module and the tqdm.auto namespace
    # so that all tqdm() calls in analysis.py emit progress events.
    _orig_ana = _ana_mod.tqdm
    _orig_data = _data_mod.tqdm
    _orig_auto = _tqdm_mod.auto.tqdm
    try:
        _ana_mod.tqdm = _TqdmStream
        _data_mod.tqdm = _TqdmStream
        _tqdm_mod.auto.tqdm = _TqdmStream

        data.verbose = True
        cols_before = set(data.df.columns)
        from af_analysis import analysis as ana

        if score == "pdockq2":
            ana.pdockq2(data)
        elif score == "LIS_LIA":
            ana.LIS_matrix(data, pae_cutoff=pae_cutoff)
            ana.cLIS_matrix(data, pae_cutoff=pae_cutoff, dist_cutoff=dist_cutoff)
        elif score == "iptm_d0":
            ana.ipTM_d0(data)
            ana.ipSAE(data, pae_cutoff=pae_cutoff)
        cols_added = [c for c in data.df.columns if c not in cols_before]
    except Exception:
        _PROGRESS_QUEUE.put_nowait(
            {"desc": "__end__", "n": 0, "total": 0, "done": True}
        )
        return (
            jsonify({"error": "Computation failed", "details": traceback.format_exc()}),
            500,
        )
    finally:
        _ana_mod.tqdm = _orig_ana
        _data_mod.tqdm = _orig_data
        _tqdm_mod.auto.tqdm = _orig_auto

    # Signal the SSE stream that computation is finished.
    _PROGRESS_QUEUE.put_nowait({"desc": "__end__", "n": 0, "total": 0, "done": True})
    return jsonify({"ok": True, "columns_added": cols_added})


@app.post("/api/cluster")
def api_cluster():
    """Run hierarchical clustering on all queries.

    Body JSON: { "threshold": 2.0 }
    Returns:   { "queries": [ { "query": str, "n_clusters": int,
                                "points": [{"row", "x", "y", "cluster"}] } ] }
    """
    import pandas as pd
    from af_analysis import clustering as _clust

    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    payload = request.get_json(silent=True) or {}
    threshold = float(payload.get("threshold") or 2.0)
    align_selection    = payload.get("align_selection")    or "backbone"
    distance_selection = payload.get("distance_selection") or "backbone"

    try:
        universes = _clust.hierarchical(
            data.df,
            threshold=threshold,
            align_selection=align_selection,
            distance_selection=distance_selection,
            show_dendrogram=False,
            MDS_coors=True,
            return_universe=True,
        )
        if isinstance(universes, dict):
            STATE.cluster_universes = universes
    except Exception:
        return (
            jsonify({"error": "Clustering failed", "details": traceback.format_exc()}),
            500,
        )

    # Build per-query response using positional (reset) row indices so they
    # match the `row` field used by the table and structure endpoints.
    df = data.df.reset_index(drop=True)
    queries = []
    for query in df["query"].unique():
        qdf = df[df["query"] == query]
        points = []
        for pos_idx, row in qdf.iterrows():
            mds1 = row.get("MDS 1")
            mds2 = row.get("MDS 2")
            clust = row.get("cluster")
            if pd.isnull(mds1) or pd.isnull(mds2):
                continue
            points.append({
                "row": int(pos_idx),
                "x": float(mds1),
                "y": float(mds2),
                "cluster": int(clust) if not pd.isnull(clust) else None,
            })
        n_clusters = int(qdf["cluster"].dropna().nunique()) if "cluster" in qdf else 0
        dendrogram_data = None
        if len(points) >= 2:
            try:
                mds_xy = np.array([[p["x"], p["y"]] for p in points])
                Y = _ssd.pdist(mds_xy)
                Z = _sch.linkage(Y, method="average")
                dend = _sch.dendrogram(Z, no_plot=True)
                row_list = [p["row"] for p in points]
                dendrogram_data = {
                    "icoord": dend["icoord"],
                    "dcoord": dend["dcoord"],
                    "leaves": [int(row_list[i]) for i in dend["leaves"]],
                    "n": len(dend["leaves"]),
                    "threshold": threshold,
                }
            except Exception:
                pass
        queries.append({"query": str(query), "points": points, "n_clusters": n_clusters, "dendrogram": dendrogram_data})

    return jsonify({"queries": queries})


@app.get("/api/superpose")
def api_superpose():
    """Return a multi-model PDB of selected rows, already aligned.

    Query params: query=<str>  rows=<comma-separated int indices>
    """
    import pandas as pd
    import MDAnalysis as mda

    query = request.args.get("query", "")
    rows_str = request.args.get("rows", "")
    try:
        row_indices = [int(r) for r in rows_str.split(",") if r.strip()]
    except ValueError:
        return jsonify({"error": "Invalid rows parameter"}), 400

    if not row_indices:
        return jsonify({"error": "No rows specified"}), 400

    if query not in STATE.cluster_universes:
        return jsonify({"error": "No clustering universe for this query. Run Clusterize first."}), 404

    u, files = STATE.cluster_universes[query]

    try:
        data = _require_data()
    except RuntimeError as exc:
        return jsonify({"error": str(exc)}), 404

    df = data.df.reset_index(drop=True)

    # Map row indices → frame indices in the universe.
    frame_indices = []
    for row_idx in row_indices:
        if row_idx >= len(df):
            continue
        pdb_path = df.iloc[row_idx].get("pdb")
        if pd.isnull(pdb_path):
            continue
        try:
            frame_indices.append(files.index(str(pdb_path)))
        except ValueError:
            pass

    if not frame_indices:
        return jsonify({"error": "No valid frames found for specified rows"}), 404

    # Write aligned frames as a multi-model PDB via a temp file (MDAnalysis PDB
    # writer requires a real filename path — it cannot write to a StringIO).
    import tempfile, os
    tmp_fd, tmp_path = tempfile.mkstemp(suffix='.pdb')
    os.close(tmp_fd)
    try:
        with mda.Writer(tmp_path, n_atoms=u.atoms.n_atoms, multiframe=True) as w:
            for fi in frame_indices:
                u.trajectory[fi]
                w.write(u.atoms)
        with open(tmp_path) as fh:
            pdb_text = fh.read()
    except Exception:
        return jsonify({"error": "Failed to write superposed PDB", "details": traceback.format_exc()}), 500
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass

    return jsonify({"structure_text": pdb_text, "structure_format": "pdb",
                    "n_frames": len(frame_indices)})


def main(argv=None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="af_analysis_gui",
        description="Launch the af_analysis web GUI.",
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=None,
        help="Path to the AlphaFold results directory to load on startup.",
    )
    parser.add_argument(
        "--host", default="127.0.0.1", help="Host to bind (default: 127.0.0.1)."
    )
    parser.add_argument(
        "--port", type=int, default=5000, help="Port to listen on (default: 5000)."
    )
    parser.add_argument(
        "--format",
        default=None,
        dest="fmt",
        help="Force a specific input format (default: auto-detect).",
    )
    args = parser.parse_args(argv)

    if args.directory:
        directory = str(Path(args.directory).expanduser().resolve())
        try:
            data = af_analysis.Data(directory=directory, format=args.fmt, verbose=True)
            with STATE.lock:
                STATE.af_data = data
                STATE.last_error = ""
            print(f"Loaded {len(data.df)} models from {directory}")
        except Exception:
            print(
                f"Warning: could not load directory '{directory}':\n{traceback.format_exc()}"
            )

    app.run(host=args.host, port=args.port, debug=False, threaded=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
