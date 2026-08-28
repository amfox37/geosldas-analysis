# ASCAT H121 / H139 data inventory (Metop-A, B, C)

Compiled 2026-08-28 from a live investigation of the H SAF FTP mirror
(`ftphsaf.meteoam.it`) — see `h121_h139_data_availability.md` for the
narrative writeup this draws on. Fields marked **TBD** were not
independently verifiable today (web documentation for these products was
either unreachable or gave numbers that conflicted with each other) and
need an OPS follow-up rather than a guess.

**Correction (same day):** an earlier version of this table had H121
Metop-B/C ending Dec 2021 — that was wrong, from only checking the 2021
year directory. Start/end dates below are from a full listing of every
year directory plus the first/last file in the boundary months. A
**"Downloaded to Discover"** column has been added, checked against the
existing local archive at `/gpfsm/dnb06/projects/p284/ASCAT_SSM_CDR/H121/`
(see its own `README.md` for the download/tarring history).

| Observation Type | Product ID | DOI | Start Date (H SAF) | End Date (H SAF) | Known Gaps (H SAF itself) | Downloaded to Discover | Latency (availability) | Data Volume | Format | OPS action needed | Current Path | Source URL | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ASCAT Surface Soil Moisture, Metop-A (CDR) | H121 | **TBD** — conflicting values found online, not confirmed against an authoritative H SAF doc | **2007-01-01** | **2021-11-15** | Nothing reachable past 2021-11-15 on either product — H139 has never carried Metop-A | **Complete** — 2007-01-01 → 2021-11-15, matches the server exactly (md5-verified) | N/A — static reprocessed archive, not updated | Not measured | NetCDF, `/h121/netcdf/metop_a/<year>/<month>/` | Confirm DOI from an official PUM/ATBD | `ASCAT_SSM_CDR/H121/metop_a_archive/metop_a_Y<2007-2014>.tar` (tarred) + `ASCAT_SSM_CDR/H121/metop_a/Y<2015-2021>/` (live) | `ftp://ftphsaf.meteoam.it/h121/netcdf/metop_a/` | Metop-A ASCAT is decommissioned — this end date is a hard archive limit, not a mirror gap |
| ASCAT Surface Soil Moisture, Metop-B (CDR) | H121 | **TBD** | **2013-06-01** | **2024-12-31** | Jan 2025 → Jun 2026 (~18 months) not reachable on this server, either product | **Partial** — 2015-01-01 → 2022-12-31 + Jan 2023 only. Missing **2013-06 → 2014-12** (18 mo, never attempted — original backfill array only covered 2015-2021) and **2023-02 → 2024-12** (23 mo, exists on server, just never downloaded) | N/A — static archive | Not measured | NetCDF, `/h121/netcdf/metop_b/<year>/<month>/` | Confirm DOI; **backfill 2013-06→2014-12 and 2023-02→2024-12 using `download_h121_metopbc_year.sh`/`_month.sh`** (proven 8-parallel-connection limit, see that project's README); decide whether/how to bridge the true Jan2025-Jun2026 HSAF-side gap | `ASCAT_SSM_CDR/H121/metop_b/Y<2015-2022>/`, `Y2023/M01/` | `ftp://ftphsaf.meteoam.it/h121/netcdf/metop_b/` | |
| ASCAT Surface Soil Moisture, Metop-C (CDR) | H121 | **TBD** | **2019-04-01** | **2024-12-31** | Same Jan 2025 → Jun 2026 gap | **Partial** — 2019-04-01 → 2022-12-31 + Jan 2023 only (2015-2018 dirs exist locally but are empty, correctly — no product exists before Apr 2019). Missing **2023-02 → 2024-12** (23 mo) | N/A — static archive | Not measured | NetCDF, `/h121/netcdf/metop_c/<year>/<month>/` | **Backfill 2023-02→2024-12**, same script/method as Metop-B; confirm DOI | `ASCAT_SSM_CDR/H121/metop_c/Y<2019-2022>/`, `Y2023/M01/` | `ftp://ftphsaf.meteoam.it/h121/netcdf/metop_c/` | |
| ASCAT Surface Soil Moisture, Metop-A (NRT/ICDR) | H139 | N/A | N/A | N/A | Not produced — no Metop-A files ever appear in the H139 window (consistent with the 2021-11-15 CDR cutoff) | N/A — product doesn't exist for this satellite | N/A | N/A | N/A | None — confirm this stays true if H SAF later reprocesses Metop-A | N/A | `ftp://ftphsaf.meteoam.it/h139/h139_cur_mon_data/` | |
| ASCAT Surface Soil Moisture, Metop-B (NRT/ICDR) | H139 | **TBD** — one search hit suggested `10.15770/EUM_SAF_H_0012` but a second source attributed that same DOI to a different resolution/product, so treat as unconfirmed | Rolling — window observed 2026-06-22 on 2026-08-28 (≈2 months back from "today") | Rolling — window observed 2026-08-23 on 2026-08-28 | **Jan 2025 → Jun 2026 (~18 months) not reachable on this server, either product**, for Metop-B/C | **Nothing retained** — only test-downloaded to scratch and deleted today, to verify the curl path works | ~1–2 days (est. from one sample: files timestamped 2026-08-24 for obs from 2026-08-22/23; not a measured statistic) | ~2.1–2.3 MB/file; roughly hourly granules per satellite in the one day sampled (24 files for 2026-08-22→23) → order ~50 MB/day/satellite, not rigorously measured | NetCDF, flat (no year/month subfolders); filenames contain embedded commas | Confirm DOI; decide whether to poll this window on a schedule to build a forward archive, since H SAF itself retains nothing older than ~2 months | None yet | `ftp://ftphsaf.meteoam.it/h139/h139_cur_mon_data/` | Must use `curl` per file, not `lftp mirror` — `lftp` hangs/`550`s on the embedded-comma filenames (see `download_hsaf_ascat.sh`) |
| ASCAT Surface Soil Moisture, Metop-C (NRT/ICDR) | H139 | **TBD** — same caveat as Metop-B | Rolling — window observed 2026-06-22 on 2026-08-28 | Rolling — window observed 2026-08-23 on 2026-08-28 | Same Jan 2025 → Jun 2026 gap | **Nothing retained** — same as Metop-B | Same as Metop-B (est., not measured) | Same order of magnitude as Metop-B (not separately sampled) | Same as Metop-B | Same as Metop-B | None yet | `ftp://ftphsaf.meteoam.it/h139/h139_cur_mon_data/` | |

## Open OPS items (rollup)

1. **Backfill H121 Metop-B/C locally** — real, available data sitting
   undownloaded: **2013-06→2014-12** (Metop-B only, 18 mo) and
   **2023-02→2024-12** (both satellites, 23 mo each). Existing scripts
   (`download_h121_metopbc_year.sh`/`_month.sh`) already handle this
   pattern; just needs new array jobs. Metop-A is already complete.
2. **DOI confirmation** — no authoritative DOI was confirmed for either
   product today; the H SAF product-detail pages and PDF docs either 404'd
   or were unreachable from here, and web search results conflicted with
   each other. Needs a direct check against the H SAF PUM/ATBD PDFs or a
   request to H SAF support.
3. **The Jan 2025 → Jun 2026 gap (Metop-B/C only)** — this is the one gap
   that's real on H SAF's own server, not just locally. H121 ends
   2024-12-31, H139's rolling window only reaches back to 2026-06-22.
   Metop-A has a separate, permanent gap from 2021-11-15 onward
   (decommissioned, never in H139). Being pursued with H SAF staff
   directly (per prior conversation).
4. **H139 has no standing local archive at all** — if a forward-looking
   record matters, this needs a scheduled/periodic pull (H SAF itself
   retains nothing older than ~2 months).
5. **Data volume** — only sampled from 2 downloaded H139 files; nothing
   measured for H121. If a real archive/staging decision depends on
   volume, worth a proper `du` after a full-month pull.
