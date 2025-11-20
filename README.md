# GEOS LDAS Analysis

This repository now contains the curated code extracted from the historical `GEOSldas_diagnostics` repo. The migration (Nov 2024) intentionally preserved only source artefacts (notebooks, Python scripts, Matlab/FORTRAN utilities, configs) that were touched within the last two years, and split them into a reusable layout:

| Location | Purpose |
| --- | --- |
| `common/python/{io,plotting,stats}` | Shared Python helpers (readers, plotting utilities, statistical functions). |
| `common/matlab`, `common/fortran` | Core Matlab/FORTRAN processing utilities reused by multiple projects. |
| `projects/<name>` | Domain-specific workspaces (e.g., `cygnss_da`, `ascat_da`, `M21C_ls`, `snow_da`, `era5_land`, `matlab2python`, `utils`). Each contains `notebooks`, `scripts`, and language-specific subfolders as needed. |
| `data/` | Documentation on how to access the archived data. No large datasets live in this repo. |

## Data handling

The heavyweight `test_data` tree from the archived repo stays outside version control. To keep notebooks runnable, create a local symlink that points back to the archived directory (update the path if the archive moves):

```bash
ln -s /Users/amfox/Desktop/GEOSldas_diagnostics/test_data data/external/test_data
```

All notebooks/scripts should reference inputs through `data/external/test_data/...`. The `data/.gitignore` prevents accidental commits of raw data; keep only lightweight docs under `data/`.

## Notebook hygiene

All `.ipynb` files have been scrubbed of execution outputs and execution counts so future diffs stay readable. Before committing new work, re-run the clearing step or use `jupyter nbconvert --ClearOutputPreprocessor.enabled=True --inplace <notebook>` to maintain the same standard.

## Future contributions

1. Place shared logic in `common/` so it can be imported across projects.
2. Drop notebooks/scripts into the relevant `projects/<name>/notebooks` or `scripts` directory and keep data-dependent paths relative to `data/`.
3. Document new workflows via per-project `README.md` files so future Codex/ChatGPT runs understand the context without reopening the archived repository.
