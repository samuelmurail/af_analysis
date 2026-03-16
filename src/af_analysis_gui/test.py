import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import base64
import subprocess
import traceback

from streamlit.components.v1 import html

import af_analysis
from af_analysis import analysis, docking


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
    st.session_state.af_data = af_analysis.Data(data_path)

    config = {
        "Preview": st.column_config.ImageColumn(),
        "Progress": st.column_config.ProgressColumn(),
    }


    tab1, tab2 = st.tabs(["Load Data", "Dataframe"])

    tab2.dataframe(st.session_state.af_data.df, column_config=config, use_container_width=True)

    with tab1:
        _render_load_panel()


if __name__ == "__main__":
    main()