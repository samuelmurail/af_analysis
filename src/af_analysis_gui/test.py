import streamlit as st
import numpy as np
import plotly.graph_objects as go
import json
from pathlib import Path
import subprocess
import traceback
from typing import Optional

import streamlit.components.v1 as components

import af_analysis


def _browse_directory() -> None:
    """Open a native directory picker with modern Linux dialog support."""

    initialdir = st.session_state.results_directory or str(Path.cwd())

    # Prefer modern desktop-native dialogs on Linux.
    result = subprocess.run(
        [
            "zenity",
            "--file-selection",
            "--directory",
            "--title=Select AlphaFold results directory",
            f"--filename={Path(initialdir).as_posix().rstrip('/')}/",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode == 0 and result.stdout.strip():
        st.session_state.results_directory = result.stdout.strip()
        return

def _load_data(directory: str, format_name: str) -> None:
    format_arg = None if format_name == "auto" else format_name
    st.session_state.af_data = af_analysis.Data(
        directory=directory,
        format=format_arg,
        verbose=False,
    )
    st.session_state.last_error = ""
    st.session_state.selected_model_index = 0
    st.session_state.selected_residue = None
    st.session_state.selected_residues = []

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

def _render_molstar_panel(
    selected_index: Optional[int] = None,
    selected_residue: Optional[int] = None,
) -> None:
    data = st.session_state.af_data
    chain_id = None
    chain_residue = None

    if selected_residue is not None:
        row = data.df.iloc[selected_index]
        query = row.get("query")
        chain_lengths = data.chain_length.get(query, []) if query is not None else []
        chain_ids = data.chains.get(query, []) if query is not None else []
        if chain_lengths and chain_ids and len(chain_lengths) == len(chain_ids):
            residue_left = int(selected_residue)
            for cid, clen in zip(chain_ids, chain_lengths):
                if residue_left <= int(clen):
                    chain_id = str(cid)
                    chain_residue = int(residue_left)
                    break
                residue_left -= int(clen)

    pdb_path = st.session_state.af_data.df.iloc[selected_index]["pdb"]
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
        structure_format = "mmcif"

    try:
        structure_text = path_obj.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        st.session_state.last_error = traceback.format_exc()
        st.error("Failed to read structure file")
        return

    viewer_height = 700
    viewer_id = f"molstar-viewer-{selected_index}"

    if selected_residue is not None:
        if chain_id is not None and chain_residue is not None:
            st.caption(
                f"Selected residue from graph: {selected_residue} "
                f"(mapped to chain {chain_id}, residue {chain_residue})"
            )
        else:
            st.caption(f"Selected residue from graph: {selected_residue}")

    html = f"""
<div id=\"{viewer_id}-root\" style=\"width: 100%; height: {viewer_height}px; position: relative;\"></div>
<link rel=\"stylesheet\" href=\"https://unpkg.com/molstar/build/viewer/molstar.css\" />
<script crossorigin src=\"https://unpkg.com/react@18/umd/react.development.js\"></script>
<script crossorigin src=\"https://unpkg.com/react-dom@18/umd/react-dom.development.js\"></script>
<script src=\"https://unpkg.com/molstar/build/viewer/molstar.js\"></script>
<script>
(function() {{
    const payload = {{
        viewerId: {json.dumps(viewer_id)},
        structureText: {json.dumps(structure_text)},
        structureFormat: {json.dumps(structure_format)},
        selectedResidue: {json.dumps(selected_residue)},
        chainId: {json.dumps(chain_id)},
        chainResidue: {json.dumps(chain_residue)},
        debugCommand: {{
            enabled: true,
            chainId: 'A',
            residue: 52,
            residueName: 'THR'
        }}
    }};

    const rootEl = document.getElementById(payload.viewerId + '-root');
    if (!rootEl) return;

    const e = React.createElement;

    function getMolstarNS() {{
        return window.molstar || {{}};
    }}

    function buildResidueSelection(plugin, residueNumber, useAuthSeq, chainId, useAuthChain) {{
        const ns = getMolstarNS();
        const MS = ns.molScript?.MolScriptBuilder || ns.MolScriptBuilder;
        const Script = ns.script?.Script || ns.Script;
        const StructureSelection = ns.structure?.StructureSelection || ns.StructureSelection;
        const data = plugin.managers.structure.hierarchy.current.structures[0]?.cell?.obj?.data;
        if (!MS || !Script || !StructureSelection || !data) return null;

        const seqProp = useAuthSeq
            ? MS.struct.atomProperty.macromolecular.auth_seq_id()
            : MS.struct.atomProperty.macromolecular.label_seq_id();

        const atomGroups = {{
            'residue-test': MS.core.set.has([MS.set(residueNumber), seqProp]),
            'group-by': MS.struct.atomProperty.macromolecular.residueKey()
        }};

        if (chainId) {{
            const chainProp = useAuthChain
                ? MS.struct.atomProperty.macromolecular.auth_asym_id()
                : MS.struct.atomProperty.macromolecular.label_asym_id();
            atomGroups['chain-test'] = MS.core.rel.eq([chainProp, chainId]);
        }}

        const query = MS.struct.generator.atomGroups(atomGroups);
        const selection = Script.getStructureSelection(query, data);
        return StructureSelection.toLociWithSourceUnits(selection);
    }}

    function getFirstLocationFromLoci(loci) {{
        const ns = getMolstarNS();
        const StructureElement = ns.structure?.StructureElement || ns.StructureElement;
        if (!StructureElement?.Loci?.getFirstLocation) return null;
        return StructureElement.Loci.getFirstLocation(loci);
    }}

    function getResolvedResidueInfo(loc) {{
        const ns = getMolstarNS();
        const StructureProperties = ns.structure?.StructureProperties || ns.StructureProperties;
        if (!StructureProperties || !loc) return null;
        return {{
            chain_label: StructureProperties.chain.label_asym_id(loc),
            chain_auth: StructureProperties.chain.auth_asym_id?.(loc),
            residue_label_seq_id: StructureProperties.residue.label_seq_id(loc),
            residue_auth_seq_id: StructureProperties.residue.auth_seq_id?.(loc),
            residue_name: String(StructureProperties.atom.label_comp_id(loc) || '').toUpperCase(),
        }};
    }}

    function highlightResidue(plugin, residueNumber, chainId, chainResidue) {{
        if (!Number.isInteger(residueNumber)) return false;

        const attempts = [];
        if (chainId && Number.isInteger(chainResidue)) {{
            attempts.push([chainResidue, true, chainId, false]);
            attempts.push([chainResidue, true, chainId, true]);
            attempts.push([chainResidue, false, chainId, false]);
            attempts.push([chainResidue, false, chainId, true]);
        }}
        attempts.push([residueNumber, true, null, false]);
        attempts.push([residueNumber, false, null, false]);

        for (const [rid, useAuthSeq, cid, useAuthChain] of attempts) {{
            const loci = buildResidueSelection(plugin, rid, useAuthSeq, cid, useAuthChain);
            if (loci && loci.elements && loci.elements.length > 0) {{
                plugin.managers.interactivity.lociHighlights.highlightOnly({{ loci }});
                plugin.managers.camera.focusLoci(loci);
                return true;
            }}
        }}
        return false;
    }}

    function selectResidueWithDebug(plugin, chainId, residueNumber, expectedResidueName) {{
        if (!Number.isInteger(residueNumber)) return false;

        const attempts = [
            [residueNumber, true, chainId, false],
            [residueNumber, true, chainId, true],
            [residueNumber, false, chainId, false],
            [residueNumber, false, chainId, true],
        ];

        for (const [rid, useAuthSeq, cid, useAuthChain] of attempts) {{
            const loci = buildResidueSelection(plugin, rid, useAuthSeq, cid, useAuthChain);
            if (!loci || !loci.elements || loci.elements.length === 0) continue;

            const loc = getFirstLocationFromLoci(loci);
            const info = getResolvedResidueInfo(loc);
            console.log('[Molstar debug command] select chain/residue', {{
                requested: {{ chainId, residueNumber, expectedResidueName }},
                attempt: {{ rid, useAuthSeq, cid, useAuthChain }},
                resolved: info,
            }});
            if (expectedResidueName && info?.residue_name !== String(expectedResidueName).toUpperCase()) {{
                console.warn('[Molstar debug command] Residue name mismatch', {{
                    expected: String(expectedResidueName).toUpperCase(),
                    resolved: info?.residue_name,
                }});
            }}

            plugin.managers.interactivity.lociSelects.select({{ loci }});
            plugin.managers.interactivity.lociHighlights.highlightOnly({{ loci }});
            plugin.managers.camera.focusLoci(loci);
            return true;
        }}
        console.warn('[Molstar debug command] No loci found', {{ chainId, residueNumber }});
        return false;
    }}

    function MolstarReactBridge(props) {{
        const viewerDivRef = React.useRef(null);
        const pluginRef = React.useRef(null);
        const loadedRef = React.useRef(false);
        const [status, setStatus] = React.useState('');

        React.useEffect(() => {{
            let disposed = false;

            async function initAndLoad() {{
                try {{
                    const viewer = await molstar.Viewer.create(props.viewerId, {{
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

                    pluginRef.current = viewer.plugin;

                    const data = await viewer.plugin.builders.data.rawData({{ data: props.structureText, label: 'model' }});
                    let trajectory = null;
                    try {{
                        trajectory = await viewer.plugin.builders.structure.parseTrajectory(data, props.structureFormat);
                    }} catch (_err) {{
                        trajectory = await viewer.plugin.builders.structure.parseTrajectory(
                            data,
                            props.structureFormat === 'pdb' ? 'mmcif' : 'pdb'
                        );
                    }}

                    const model = await viewer.plugin.builders.structure.createModel(trajectory);
                    const structure = await viewer.plugin.builders.structure.createStructure(model);
                    const polymer = await viewer.plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
                    if (polymer) {{
                        await viewer.plugin.builders.structure.representation.addRepresentation(polymer, {{
                            type: 'cartoon',
                            color: 'chain-id'
                        }});
                        viewer.plugin.managers.camera.reset();
                    }}

                    loadedRef.current = true;
                    if (!disposed) setStatus('');
                }} catch (err) {{
                    if (!disposed) setStatus('Mol* could not render this structure: ' + err);
                }}
            }}

            initAndLoad();

            return () => {{
                disposed = true;
                try {{
                    if (pluginRef.current) pluginRef.current.dispose();
                }} catch (_e) {{}}
                pluginRef.current = null;
                loadedRef.current = false;
            }};
        }}, [props.viewerId, props.structureText, props.structureFormat]);

        React.useEffect(() => {{
            if (!loadedRef.current || !pluginRef.current) return;
            if (!Number.isInteger(props.selectedResidue)) return;

            const ok = highlightResidue(
                pluginRef.current,
                props.selectedResidue,
                props.chainId,
                props.chainResidue
            );
            if (!ok) {{
                setStatus('Residue highlight not found for index ' + props.selectedResidue);
            }} else {{
                setStatus('');
            }}
        }}, [props.selectedResidue, props.chainId, props.chainResidue]);

        // First-check command: execute a direct Mol* selection command after load,
        // without reinitializing the viewer in this render.
        const debugCommandRanRef = React.useRef(false);
        React.useEffect(() => {{
            if (!loadedRef.current || !pluginRef.current) return;
            if (debugCommandRanRef.current) return;

            const cmd = props.debugCommand;
            if (!cmd?.enabled) return;

            debugCommandRanRef.current = true;
            const ok = selectResidueWithDebug(
                pluginRef.current,
                String(cmd.chainId || ''),
                Number(cmd.residue),
                cmd.residueName,
            );
            if (!ok) {{
                setStatus('Debug command could not find ' + String(cmd.chainId || '?') + ':' + String(cmd.residue || '?'));
            }}
        }}, [props.debugCommand]);

        return e(
            'div',
            {{ style: {{ width: '100%', height: '100%', position: 'relative' }} }},
            e('div', {{ id: props.viewerId, ref: viewerDivRef, style: {{ width: '100%', height: '100%' }} }}),
            status
                ? e(
                    'div',
                    {{
                        style: {{
                            position: 'absolute',
                            left: '8px',
                            bottom: '8px',
                            background: 'rgba(255,255,255,0.9)',
                            padding: '4px 8px',
                            fontSize: '12px',
                            borderRadius: '6px'
                        }}
                    }},
                    status
                )
                : null
        );
    }}

    const root = ReactDOM.createRoot(rootEl);
    root.render(e(MolstarReactBridge, payload));
}})();
</script>
"""

    components.html(html, height=viewer_height + 20)

def _plot_plddt(selected_index):
    fig_plddt, _ = st.session_state.af_data.plot_plddt([selected_index])
    st.pyplot(fig_plddt)


def _plot_plddt_interactive(selected_index: int) -> Optional[int]:
    plddt_array = st.session_state.af_data.get_plddt(selected_index)
    if plddt_array is None or len(plddt_array) == 0:
        st.warning("No pLDDT values available for this model")
        return None

    residues = np.arange(1, len(plddt_array) + 1)
    fig = go.Figure(
        data=[
            go.Scatter(
                x=residues,
                y=plddt_array,
                mode="lines+markers",
                marker={"size": 4},
                line={"width": 2},
                name="pLDDT",
            )
        ]
    )
    fig.update_layout(
        margin={"l": 20, "r": 20, "t": 10, "b": 30},
        xaxis_title="Residue",
        yaxis_title="predicted LDDT",
        yaxis={"range": [0, 100]},
        dragmode="pan",
        clickmode="event+select",
    )

    st.caption("Click one or more points to select residues")

    event = st.plotly_chart(
        fig,
        use_container_width=True,
        key=f"plddt_interactive_{selected_index}",
        on_select="rerun",
        selection_mode=("points",),
    )

    selected_residues = st.session_state.get("selected_residues", [])

    if event is None or event.selection is None:
        st.write("Selected residues:", selected_residues)
        return st.session_state.get("selected_residue")

    points = event.selection.get("points", [])
    if not points:
        st.write("Selected residues:", selected_residues)
        return st.session_state.get("selected_residue")

    try:
        selected_residues = sorted(
            {
                int(point["x"])
                for point in points
                if point.get("x") is not None
            }
        )
        st.session_state.selected_residues = selected_residues
        st.write("Selected residues:", selected_residues)
        if selected_residues:
            return selected_residues[0]
        return st.session_state.get("selected_residue")
    except (KeyError, TypeError, ValueError):
        return st.session_state.get("selected_residue")

def _plot_pae(selected_index):
    fig_pae, _ = st.session_state.af_data.plot_pae(selected_index)
    st.pyplot(fig_pae)

def _plot_graph_model(score, selected_index) -> Optional[int]:
    if score == "plddt":
        return _plot_plddt_interactive(selected_index)
    if score == "pae":
        _plot_pae(selected_index)
        return None
    return None

def _render_graph_panel(selected_index: Optional[int] = None) -> Optional[int]:
    graph_name = st.selectbox(
        "Score",
        options=[
            "plddt",
            "pae",
        ],
        index=0,
        label_visibility = "collapsed",#("visible", "hidden", or "collapsed")
    )

    return _plot_graph_model(graph_name, selected_index)

def _render_explorer_panel() -> None:
    data = st.session_state.af_data
    if data is None:
        st.info("Load a dataset first")
        return
    if len(data.df) == 0:
        st.warning("Dataset is empty")
        return

    col_left, col_right = st.columns([4, 2])
    with col_right:
        config = {
            "Preview": st.column_config.ImageColumn(),
            "Progress": st.column_config.ProgressColumn(),
        }
        col_right.dataframe(st.session_state.af_data.df, column_config=config, use_container_width=True)
        # selected_index = _render_table_selector(data)
    with col_left:
        graph_col, molstar_col = st.columns([2, 2])

        selected_index = 0

        with graph_col:
            selected_residue = _render_graph_panel(selected_index=selected_index)
            st.session_state.selected_residue = selected_residue
        with molstar_col:
            _render_molstar_panel(
                selected_index=selected_index,
                selected_residue=st.session_state.get("selected_residue"),
            )


def main():
    # Header:
    st.set_page_config(page_title="af_analysis GUI", layout="wide")


    with st.container(key="app_title"):
        st.title("Af Analysis GUI")

    logo_path = "https://raw.githubusercontent.com/samuelmurail/af_analysis/master/docs/source/logo.jpeg"

    st.markdown(
        f"""
        <style>
        .stDeployButton,  .stAppDeployButton{{
            visibility: hidden;
        }}

        header.stAppHeader {{
            background: rgba(0, 0, 0, 0);
        }}

        div.stMainBlockContainer {{
            padding-top: 0.2rem;
        }}

        h1#af-analysis-gui {{
            padding-left: 100px;
            text-align: left;
            font-size: 1.8em;
            background-image: url("{logo_path}");
            background-size: contain;
            background-repeat: no-repeat;
            background-position: left 10px center;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

    data_path = "/home/murail/Documents/Code/af_analysis/tests/inputs/fold_2024_07_01_12_14_prot_dna_zn/"
    if "af_data" not in st.session_state:
        st.session_state.af_data = af_analysis.Data(data_path)
    if "selected_residue" not in st.session_state:
        st.session_state.selected_residue = None
    if "selected_residues" not in st.session_state:
        st.session_state.selected_residues = []


    tab1, tab2 = st.tabs(["Load Data", "Explore"])


    with tab1:
        _render_load_panel()
    with tab2:
        _render_explorer_panel()


if __name__ == "__main__":
    main()