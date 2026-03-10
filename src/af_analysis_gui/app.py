"""Streamlit GUI for af_analysis."""

from __future__ import annotations

import json
import traceback
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

import af_analysis
from af_analysis import analysis, docking


def _init_state() -> None:
    if "af_data" not in st.session_state:
        st.session_state.af_data = None
    if "last_error" not in st.session_state:
        st.session_state.last_error = ""
    if "results_directory" not in st.session_state:
        st.session_state.results_directory = ""


def _load_data(directory: str, format_name: str) -> None:
    format_arg = None if format_name == "auto" else format_name
    st.session_state.af_data = af_analysis.Data(
        directory=directory,
        format=format_arg,
        verbose=False,
    )
    st.session_state.last_error = ""


def _safe_metric_call(label: str, fn) -> None:
    with st.spinner(f"Running {label}..."):
        try:
            fn()
            st.success(f"{label} completed")
        except Exception:  # pragma: no cover
            st.session_state.last_error = traceback.format_exc()
            st.error(f"{label} failed")


def _browse_directory() -> None:
    """Open a native directory picker when available."""

    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        selected = filedialog.askdirectory(
            title="Select AlphaFold results directory",
            initialdir=st.session_state.results_directory or ".",
        )
        root.destroy()

        if selected:
            st.session_state.results_directory = selected
    except Exception:
        st.session_state.last_error = traceback.format_exc()
        st.warning("Directory picker is unavailable in this environment. Enter the path manually.")


def _render_load_panel() -> None:
    st.subheader("Load AlphaFold results")
    col_dir, col_browse = st.columns([6, 1])
    with col_dir:
        directory = st.text_input(
            "Results directory",
            key="results_directory",
        )
    with col_browse:
        st.write("")
        st.write("")
        st.button("Browse", on_click=_browse_directory)

    format_name = st.selectbox(
        "Format",
        options=[
            "auto",
            "default",
            "colabfold_1.5",
            "AF3_local",
            "AF3_webserver",
            "alphapulldown",
            "alphapulldown_full",
            "boltz1",
            "chai1",
            "massivefold",
            "full_massivefold",
        ],
        index=0,
    )

    if st.button("Load dataset", type="primary"):
        if not directory.strip():
            st.warning("Please provide a directory path")
            return
        try:
            _load_data(directory.strip(), format_name)
            st.success("Dataset loaded")
        except Exception:  # pragma: no cover
            st.session_state.last_error = traceback.format_exc()
            st.error("Failed to load dataset")


def _render_dataset_panel() -> None:
    data = st.session_state.af_data
    if data is None:
        st.info("Load a dataset first")
        return

    st.subheader("Dataset")
    st.write(f"Rows: {len(data.df)}")
    st.dataframe(data.df)

    csv_bytes = data.df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "Download CSV",
        data=csv_bytes,
        file_name="af_analysis_results.csv",
        mime="text/csv",
    )


def _render_metrics_panel() -> None:
    data = st.session_state.af_data
    if data is None:
        st.info("Load a dataset first")
        return

    st.subheader("Metrics")

    c1, c2, c3 = st.columns(3)
    if c1.button("Compute pDockQ"):
        _safe_metric_call("pDockQ", lambda: analysis.pdockq(data, verbose=False))
    if c2.button("Compute pDockQ2"):
        _safe_metric_call("pDockQ2", lambda: analysis.pdockq2(data, verbose=False))
    if c3.button("Compute LIS_pep"):
        _safe_metric_call("LIS_pep", lambda: docking.LIS_pep(data, verbose=False))

    st.write("Updated columns")
    st.dataframe(data.df.tail(min(20, len(data.df))))


def _render_model_panel() -> None:
    data = st.session_state.af_data
    if data is None:
        st.info("Load a dataset first")
        return

    st.subheader("Model plots")
    if len(data.df) == 0:
        st.warning("Dataset is empty")
        return

    index = st.number_input(
        "Model index",
        min_value=0,
        max_value=max(0, len(data.df) - 1),
        value=0,
        step=1,
    )

    try:
        fig_plddt, _ = data.plot_plddt([int(index)])
        st.pyplot(fig_plddt)
    except Exception:
        st.warning("Could not render pLDDT for this model")

    try:
        pae_result = data.plot_pae(int(index))
        if pae_result is None:
            st.info("No PAE data available for this model")
        else:
            fig_pae, _ = pae_result
            st.pyplot(fig_pae)
    except Exception:
        st.warning("Could not render PAE for this model")


def _render_molstar_panel() -> None:
    data = st.session_state.af_data
    if data is None:
        st.info("Load a dataset first")
        return

    st.subheader("3D viewer (Mol*)")
    if len(data.df) == 0:
        st.warning("Dataset is empty")
        return

    index = st.number_input(
        "Model index for 3D",
        min_value=0,
        max_value=max(0, len(data.df) - 1),
        value=0,
        step=1,
        key="molstar_model_index",
    )
    color_mode = st.radio(
        "Color",
        options=["pLDDT", "Chain"],
        horizontal=True,
        key="molstar_color_mode",
    )
    height = st.slider("Viewer height", min_value=450, max_value=1000, value=700)

    row = data.df.iloc[int(index)]
    pdb_path = row.get("pdb")
    if not pdb_path:
        st.info("No PDB file available for this model")
        return

    path_obj = Path(str(pdb_path))
    if not path_obj.is_file():
        st.warning(f"PDB file not found: {path_obj}")
        return

    structure_ext = path_obj.suffix.lower()
    if structure_ext in {".cif", ".mmcif"}:
        structure_format = "mmcif"
    elif structure_ext == ".pdb":
        structure_format = "pdb"
    else:
        # Fallback: let Mol* try mmCIF first for unknown extensions.
        structure_format = "mmcif"

    try:
        structure_text = path_obj.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        st.session_state.last_error = traceback.format_exc()
        st.error("Failed to read structure file")
        return

    # JSON encoding safely escapes file content for JavaScript injection.
    structure_payload = json.dumps(structure_text)
    structure_format_payload = json.dumps(structure_format)
    color_theme = "uncertainty" if color_mode == "pLDDT" else "chain-id"
    html = f"""
<div id=\"molstar-viewer\" style=\"width: 100%; height: {height}px; position: relative;\"></div>
<link rel=\"stylesheet\" href=\"https://unpkg.com/molstar/build/viewer/molstar.css\" />
<script src=\"https://unpkg.com/molstar/build/viewer/molstar.js\"></script>
<script>
(async function() {{
  const target = document.getElementById('molstar-viewer');
  if (!target) return;

  target.innerHTML = '';
    const structureText = {structure_payload};
    const structureFormat = {structure_format_payload};
    const colorTheme = {json.dumps(color_theme)};

  try {{
    const viewer = await molstar.Viewer.create('molstar-viewer', {{
      layoutIsExpanded: false,
            layoutShowControls: false,
            layoutShowSequence: false,
      layoutShowLog: false,
            layoutShowLeftPanel: false,
            viewportShowExpand: false,
            viewportShowSelectionMode: false,
      viewportShowAnimation: false,
      pdbProvider: 'pdbe',
      emdbProvider: 'pdbe'
    }});

        const data = await viewer.plugin.builders.data.rawData({{ data: structureText, label: 'model' }});
        let trajectory = null;
        try {{
            trajectory = await viewer.plugin.builders.structure.parseTrajectory(data, structureFormat);
        }} catch (_err) {{
            // Fallback for mismatched extension/content.
            trajectory = await viewer.plugin.builders.structure.parseTrajectory(data, structureFormat === 'pdb' ? 'mmcif' : 'pdb');
        }}
        const model = await viewer.plugin.builders.structure.createModel(trajectory);
        const structure = await viewer.plugin.builders.structure.createStructure(model);
        const polymer = await viewer.plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');

        if (polymer) {{
            await viewer.plugin.builders.structure.representation.addRepresentation(polymer, {{
                type: 'cartoon',
                color: colorTheme
            }});
            viewer.plugin.managers.camera.reset();
        }}
  }} catch (err) {{
    target.innerText = 'Mol* could not render this structure: ' + err;
  }}
}})();
</script>
"""

    components.html(html, height=height + 20)


def main() -> None:
    st.set_page_config(page_title="af_analysis GUI", layout="wide")
    _init_state()

    st.title("af_analysis GUI")
    st.caption("Interactive analysis for AlphaFold outputs")

    tabs = st.tabs(["Load", "Dataset", "Metrics", "Model", "3D"])
    with tabs[0]:
        _render_load_panel()
    with tabs[1]:
        _render_dataset_panel()
    with tabs[2]:
        _render_metrics_panel()
    with tabs[3]:
        _render_model_panel()
    with tabs[4]:
        _render_molstar_panel()

    if st.session_state.last_error:
        with st.expander("Last error", expanded=False):
            st.code(st.session_state.last_error)


if __name__ == "__main__":
    main()
