"""CLI launcher for af_analysis Streamlit GUI."""

from __future__ import annotations

import sys
from pathlib import Path

from streamlit.web.cli import main as st_main


def main() -> int:
    app_path = str(Path(__file__).with_name("app.py"))
    sys.argv = ["streamlit", "run", app_path]
    return st_main()


if __name__ == "__main__":
    raise SystemExit(main())
