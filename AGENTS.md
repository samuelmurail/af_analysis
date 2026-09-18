# AGENTS.md

Practical instructions for AI coding agents working in this repository.

## Purpose

`af-analysis` is a Python package for analyzing AlphaFold/ColabFold/AF3 outputs, with an optional Flask GUI.

## Start Here

1. Install in editable mode:
   - `pip install -e .`
2. Run tests before and after changes:
   - `pytest`
   - or `make test`
3. Format touched Python files:
   - `black src/af_analysis src/af_analysis_gui tests`

## High-Value Commands

- `pytest` - run all tests.
- `pytest tests/test_data.py -k <pattern>` - run focused tests.
- `make test` - test shortcut.
- `make lint` - runs Black (currently on `src/af_analysis`).
- `make simple-docs` - rebuild Sphinx docs.
- `pip install "af-analysis[gui]"` then `af_analysis_gui` - run GUI locally.

## Codebase Map

- `src/af_analysis/` - core library.
  - `data.py` defines the `Data` class and format auto-detection.
  - `analysis.py`, `docking.py`, `clustering.py`, `plot.py`, `sequence.py` contain metrics and analysis helpers.
  - `format/` contains format-specific parsers (AF2/AF3/ColabFold/Boltz1/Chai1/MassiveFold/AlphaPulldown).
- `src/af_analysis_gui/` - optional Flask web app (`flask_app.py`, templates, static assets).
- `tests/` - pytest suite and real input fixtures in `tests/inputs/`.
- `docs/source/` - Sphinx source docs.

## Working Rules

1. Edit source under `src/`, not `build/lib/`.
2. Keep parser behavior backwards compatible: many tests depend on exact dataframe columns and expected metrics.
3. When adding support for a new result format, mirror patterns in `src/af_analysis/format/default.py` and add tests in `tests/test_format.py` with representative fixtures.
4. Prefer small, targeted changes; avoid broad refactors unless requested.
5. If API behavior changes, update docs and tests in the same change.

## Testing Expectations

- Always run `pytest` for full validation.
- For parser or IO work, run at least:
  - `pytest tests/test_data.py tests/test_format.py`
- For metric changes, run at least:
  - `pytest tests/test_analysis.py tests/test_docking.py tests/test_clustering.py`

## Common Pitfalls

- `Data.read_directory()` uses file/folder heuristics for format detection. Preserve existing detection order unless explicitly changing behavior.
- Test fixtures are large and format-specific; avoid renaming fixture files/directories unless tests are updated accordingly.
- `make lint` currently targets `src/af_analysis`; if you modify GUI/tests, run Black directly on those paths too.
- `pyproject.toml` metadata has some historical inconsistencies (for example classifier/license mismatch). Do not "fix" packaging metadata unless the task asks for packaging updates.

## Existing Docs (Link, Do Not Duplicate)

- Project overview and usage: [README.md](README.md)
- Installation details: [docs/source/installation.md](docs/source/installation.md)
- Contribution process and style: [docs/source/contributing.md](docs/source/contributing.md)
- Code of conduct: [docs/source/code_of_conduct.md](docs/source/code_of_conduct.md)

## PR Readiness Checklist

1. Tests pass: `pytest`
2. Formatting is clean: `black <touched paths>`
3. Docs updated when user-facing behavior changes
4. No edits in generated/build artifacts unless explicitly requested