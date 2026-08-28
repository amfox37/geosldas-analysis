# H121 / H139 data availability on the H SAF FTP mirror

Investigated 2026-08-28 while extending `scripts/download_hsaf_ascat.sh` to
support H139. Server: `ftphsaf.meteoam.it` (credentials in `~/.netrc`,
`machine ftphsaf.meteoam.it`).

**Correction (same day):** the first pass at this only checked the 2021
year directory for each satellite and wrongly concluded Metop-B/C stopped
in Dec 2021. A full listing of every year directory shows they continue
through Dec 2024. Table below and the gap section have been corrected.

## H121 (CDR)

Full reprocessed archive, organized `/h121/netcdf/<sat>/<year>/<month>/`.
Confirmed start and end of record per satellite by listing every year
directory and the first/last file in the boundary months (directory
listing, not documentation):

| Satellite | First observation | Last observation |
|---|---|---|
| Metop-A | 2007-01-01 | 2021-11-15 (matches its real decommissioning) |
| Metop-B | 2013-06-01 | 2024-12-31 |
| Metop-C | 2019-04-01 | 2024-12-31 |

No H121 data exists outside these ranges on this server.

## H139 (near-real-time counterpart)

**Not a continuous archive.** Server layout is completely different from
H121: a single flat directory `/h139/h139_cur_mon_data/`, no
year/month/satellite subfolders, all satellites mixed together by filename.
It holds only a rolling window of the current + prior month. Checked
2026-08-28: window covered **2026-06-22 → 2026-08-23**, Metop-B and Metop-C
only (no Metop-A, consistent with H121's Metop-A cutoff).

## The gap

**Metop-A:** nothing reachable past 2021-11-15 on either product - H121
stops there (decommissioned) and H139 has never carried Metop-A.

**Metop-B/C:** H121 runs through 2024-12-31; H139's rolling window only
reaches back to 2026-06-22. So **Jan 2025 → Jun 2026 (~18 months) is not
reachable through this FTP path on either product** for these two
satellites. Checked and ruled out before concluding this:

- Server's own `README_H-SAF_Product_List.txt` and `Product User Manuals/` -
  stale, doesn't list H121 or H139 at all
- Local Discover storage - no H139 data has ever been downloaded here
- `iv_tc` manifests referencing "H121_H139" - just per-day test fixtures from
  Oct 2018 (within the H121 era), not evidence of a broader H139 archive

If a fuller H139 ICDR archive exists (e.g. via a separate H SAF distribution
channel or the EUMETSAT Data Store), it is not exposed through this FTP
mirror's `h139_cur_mon_data` folder. Andrew is pursuing this directly with
H SAF staff as of 2026-08-28.

## Script support

`scripts/download_hsaf_ascat.sh` now handles both products (`-p h121` /
`-p h139`) but they behave differently - see the script's `--help` for
the exact layout, filtering, and known-limitation details. H139 downloads
route through `curl` per file rather than `lftp`, because `lftp`'s own
`get`/`mirror --include` hang or fail with a `550` on this product's
filenames (they contain embedded commas, e.g.
`W_IT-HSAF-ROME,SAT,SSM-ASCAT-METOPB-...`); `curl` handles them fine.
