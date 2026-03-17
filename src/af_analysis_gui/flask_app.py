"""Flask GUI for af_analysis with interactive Plotly + Mol* panels."""

from __future__ import annotations

import traceback
from pathlib import Path
from threading import Lock
from typing import Any

import numpy as np
from flask import Flask, jsonify, render_template, request

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
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


class AppState:
    def __init__(self) -> None:
        self.lock = Lock()
        self.af_data: af_analysis.Data | None = None
        self.last_error = ""


STATE = AppState()

app = Flask(
    __name__,
    template_folder=str(APP_DIR / "templates"),
    static_folder=str(APP_DIR / "static"),
)


def _require_data() -> af_analysis.Data:
    if STATE.af_data is None:
        raise RuntimeError("Dataset is not loaded")
    return STATE.af_data


def _preferred_table_columns(df) -> list[str]:
    preferred = [
        "query",
        "model",
        "name",
        "rank",
        "ranking_score",
        "ptm",
        "iptm",
        "mean_plddt",
        "pdockq",
        "pdockq2",
        "LIS_pep",
    ]
    selected = [column for column in preferred if column in df.columns]
    if len(selected) >= 3:
        return selected[:8]
    return list(df.columns[:8])


@app.get("/")
def index():
    return render_template("flask_index.html")


@app.post("/api/load")
def api_load_dataset():
    payload = request.get_json(silent=True) or {}
    directory = str(payload.get("directory", "")).strip()
    format_name = str(payload.get("format", "auto")).strip() or "auto"

    if not directory:
        return jsonify({"error": "Results directory is required"}), 400

    try:
        format_arg = None if format_name == "auto" else format_name
        data = af_analysis.Data(directory=directory, format=format_arg, verbose=False)
        with STATE.lock:
            STATE.af_data = data
            STATE.last_error = ""
        return jsonify({"ok": True, "rows": int(len(data.df))})
    except Exception:
        with STATE.lock:
            STATE.last_error = traceback.format_exc()
        return jsonify({"error": "Failed to load dataset", "details": STATE.last_error}), 500


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

    return jsonify({"columns": ["row", *columns], "rows": rows, "total": int(len(df))})


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

    return jsonify(
        {
            "residues": residues,
            "plddt": plddt,
            "chain_ids": chain_ids,
            "chain_lengths": chain_lengths,
            "chain_boundaries": boundaries,
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
    if not data_file or (hasattr(data_file, "__class__") and str(data_file) in ("", "nan")):
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
        return jsonify({"error": "Failed to read structure file", "details": traceback.format_exc()}), 500

    return jsonify({"structure_text": structure_text, "structure_format": structure_format})


@app.get("/api/health")
def api_health():
    loaded = STATE.af_data is not None
    rows = int(len(STATE.af_data.df)) if loaded else 0
    return jsonify({"loaded": loaded, "rows": rows})


def main() -> int:
    app.run(host="127.0.0.1", port=5000, debug=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
