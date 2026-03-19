"""CLI launcher for af_analysis Flask GUI."""

from __future__ import annotations
import sys

from af_analysis_gui.flask_app import main


def main_cli():
    raise SystemExit(main(sys.argv[1:]))


if __name__ == "__main__":
    main_cli()
