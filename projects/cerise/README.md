# CERISE

Workspace for CERISE-related analysis, notebooks, scripts, and reporting.

## Public Deliverables Archive

The CERISE public deliverables archive is stored in:

`data/public_deliverables/`

It was built from the public CERISE deliverables page:

`https://www.cerise-project.eu/deliverables`

Run or refresh the archive with:

```bash
/Users/amfox/mamba/envs/regrid/bin/python projects/cerise/scripts/cerise_public_deliverables.py \
  --out projects/cerise/data/public_deliverables
```

Use `--dry-run` first to update the manifest without downloading files:

```bash
/Users/amfox/mamba/envs/regrid/bin/python projects/cerise/scripts/cerise_public_deliverables.py \
  --out projects/cerise/data/public_deliverables \
  --dry-run
```

The archive includes:

- `00_manifest/cerise_deliverables_manifest.csv` - one row per deliverable, including confidential/skipped items.
- `00_manifest/cerise_deliverables.bib` - starter BibTeX entries using `@techreport`.
- `00_manifest/README.md` - archive summary for the current run.
- `00_manifest/download_log.txt` - link discovery and run summary.
- `WP*/D*/` - downloaded public deliverable PDFs organized by work package and deliverable.

Manifest status values:

- `downloaded` - public deliverable file downloaded during the run.
- `exists` - local file already existed and was reused.
- `planned` - dry-run path for a file that would be downloaded.
- `no_file_link_found` - public deliverable listed on the page, but no downloadable file was found.
- `skipped_confidential` - confidential deliverable retained in the manifest but not downloaded.
- `error: ...` - download failed or returned unexpected content.

Unlike most project data directories in this repository, the CERISE public
deliverable PDFs and manifest files are intentionally committed. These are
public grey-literature technical reports without obvious DOI records, and the
local archive preserves a citable project-document snapshot.

## Layout

- `notebooks/` - exploratory and analysis notebooks.
- `scripts/` - reusable processing or plotting scripts.
- `report/` - figures, notes, and report-ready outputs.
- `data/` - local inputs or derived data products. The public CERISE
  deliverables archive is an intentional exception to the usual data-ignore
  pattern.
