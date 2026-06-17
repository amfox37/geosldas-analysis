#!/usr/bin/env python3
"""
Download and organize public CERISE deliverables.

What it does
------------
1. Reads the CERISE deliverables page.
2. Uses a curated deliverable metadata table matching the public CERISE page.
3. Downloads only public deliverables with discoverable document links.
4. Organizes files by Work Package and Deliverable.
5. Writes:
   - 00_manifest/cerise_deliverables_manifest.csv
   - 00_manifest/cerise_deliverables.bib
   - 00_manifest/README.md
   - 00_manifest/download_log.txt

Notes
-----
CERISE deliverables appear to be project reports / grey literature, not DOI-based papers.
Some public deliverables have future due dates or may not yet have a file posted.
Those are retained in the manifest with status "no_file_link_found".

Usage
-----
    python cerise_public_deliverables.py --dry-run
    python cerise_public_deliverables.py
    python cerise_public_deliverables.py --out CERISE

Dependencies
------------
    python -m pip install requests beautifulsoup4
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import re
import sys
from dataclasses import dataclass, asdict
from datetime import date
from pathlib import Path
from urllib.parse import urljoin, urlparse

import requests
from bs4 import BeautifulSoup


BASE_URL = "https://www.cerise-project.eu"
DELIVERABLES_URL = "https://www.cerise-project.eu/deliverables"

WP_FOLDER_NAMES = {
    "WP1": "WP1_Land_DA_Methodology",
    "WP2": "WP2_Coupled_Surface_Atmosphere_DA",
    "WP3": "WP3_Seasonal_Initialization",
    "WP4": "WP4_Reanalysis_Prototypes",
    "WP5": "WP5_Seasonal_Forecast_Demonstrators",
    "WP6": "WP6_Demonstrator_Evaluation",
    "WP7": "WP7_Land_Atmosphere_Interface_Observations",
    "WP8": "WP8_Project_Management_Dissemination",
}


@dataclass(frozen=True)
class Deliverable:
    did: str
    wp: str
    wp_title: str
    title: str
    dtype: str
    dissemination: str
    est_date: str


# Metadata from the CERISE deliverables page as checked 2026-06-17.
# Keep the confidential ones here too so the manifest can explicitly say they were skipped.
DELIVERABLES: list[Deliverable] = [
    Deliverable("D1.1", "WP1", "Land Data Assimilation Methodology for reanalysis",
                "Preliminary assessment of ensemble perturbation methods for the land-surface assimilation systems",
                "REPORT", "Public", "Dec-2023"),
    Deliverable("D1.2", "WP1", "Land Data Assimilation Methodology for reanalysis",
                "Unified, ensemble-based global land data assimilation system and documentation",
                "REPORT", "Public", "Dec-2024"),
    Deliverable("D1.3", "WP1", "Land Data Assimilation Methodology for reanalysis",
                "Unified, ensemble-based regional land data assimilation system and documentation",
                "REPORT", "Public", "Dec-2024"),
    Deliverable("D1.4", "WP1", "Land Data Assimilation Methodology for reanalysis",
                "Report on observation operator methodology ready for implementation in coupled global and regional systems",
                "REPORT", "Public", "Dec-2025"),

    Deliverable("D2.1", "WP2", "Coupled surface-atmosphere assimilation for global and reanalysis systems",
                "Documentation of coupled assimilation infrastructure and methodology and preliminary assessment towards optimal degrees of coupling for coupled global reanalysis",
                "REPORT", "Public", "Dec-2024"),
    Deliverable("D2.2", "WP2", "Coupled surface-atmosphere assimilation for global and reanalysis systems",
                "Documentation of coupled assimilation methodology and preliminary assessment towards optimal degrees of coupling for regional reanalysis",
                "REPORT", "Public", "Dec-2024"),
    Deliverable("D2.3", "WP2", "Coupled surface-atmosphere assimilation for global and reanalysis systems",
                "Documentation on coupled skin temperature assimilation for coupled reanalysis",
                "REPORT", "Public", "Dec-2025"),
    Deliverable("D2.4", "WP2", "Coupled surface-atmosphere assimilation for global and reanalysis systems",
                "Documentation on coupled skin temperature assimilation for regional reanalyses",
                "REPORT", "Public", "Dec-2025"),
    Deliverable("D2.5", "WP2", "Coupled surface-atmosphere assimilation for global and reanalysis systems",
                "Documentation on next reanalysis generation coupled assimilation systems",
                "REPORT", "Public", "Dec-2026"),

    Deliverable("D3.1", "WP3", "Balanced initialization of land surface for seasonal forecasts",
                "Documentation of the intermediate set of land surface initialisation systems",
                "REPORT", "Public", "Sep-2025"),
    Deliverable("D3.2", "WP3", "Balanced initialization of land surface for seasonal forecasts",
                "One or more sets of land surface initial conditions for 1993-2022 for use in seasonal forecast demonstrators",
                "DATA", "Confidential", "Sep-2025"),
    Deliverable("D3.3", "WP3", "Balanced initialization of land surface for seasonal forecasts",
                "Monthly mean land surface data from initialisation systems used for D3.2, to allow intercomparison and further assessment",
                "DATA", "Confidential", "Sep-2025"),
    Deliverable("D3.4", "WP3", "Balanced initialization of land surface for seasonal forecasts",
                "Monthly mean land surface data from the real-time initialisation systems - for further assessment",
                "DATA", "Confidential", "Mar-2026"),
    Deliverable("D3.5", "WP3", "Balanced initialization of land surface for seasonal forecasts",
                "Documentation of the final version of land surface initialisation systems, including real-time prototypes",
                "REPORT", "Public", "Jun-2026"),

    Deliverable("D4.1", "WP4", "Reanalysis prototypes",
                "Deliver ERA6-Land-Pv2 Land prototype to provide a basis for ERA6-Land",
                "DATA", "Confidential", "Sep-2025"),
    Deliverable("D4.2", "WP4", "Reanalysis prototypes",
                "Deliver ERA7-Pv2 second prototype as preparation for ERA7",
                "DATA", "Confidential", "Dec-2026"),
    Deliverable("D4.3", "WP4", "Reanalysis prototypes",
                "Deliver CARRA-land-Pv2, second prototype of an offline Arctic land system",
                "DATA", "Confidential", "Aug-2025"),
    Deliverable("D4.4", "WP4", "Reanalysis prototypes",
                "Deliver CERRA-land-Pv1 reanalysis over Europe",
                "DATA", "Confidential", "Dec-2025"),
    Deliverable("D4.5", "WP4", "Reanalysis prototypes",
                "Deliver CARRA3-Pv1, prototype for CARRA3",
                "DATA", "Confidential", "Dec-2025"),
    Deliverable("D4.6", "WP4", "Reanalysis prototypes",
                "Deliver CERRA2-Pv1, prototype for CERRA2",
                "DATA", "Confidential", "Sep-2026"),

    Deliverable("D5.1", "WP5", "Seasonal forecast demonstrators",
                "Output of CMCC demonstrators in the data archive",
                "DATA", "Confidential", "Dec-2025"),
    Deliverable("D5.2", "WP5", "Seasonal forecast demonstrators",
                "Output of DWD demonstrators in the data archive",
                "DATA", "Confidential", "Dec-2025"),
    Deliverable("D5.3", "WP5", "Seasonal forecast demonstrators",
                "Output of ECMWF demonstrators in the data archive",
                "DATA", "Confidential", "Dec-2025"),
    Deliverable("D5.4", "WP5", "Seasonal forecast demonstrators",
                "Output of MF demonstrators in the data archive",
                "DATA", "Confidential", "Dec-2025"),
    Deliverable("D5.5", "WP5", "Seasonal forecast demonstrators",
                "Output of MetO demonstrators in the data archive",
                "DATA", "Confidential", "Dec-2025"),
    Deliverable("D5.6", "WP5", "Seasonal forecast demonstrators",
                "Conclusions of the common assessment of the impact of improvements in land initialization on forecast performance",
                "REPORT", "Public", "Oct-2026"),

    Deliverable("D6.1", "WP6", "Evaluation and exploitation of demonstrator results for future C3S implementations",
                "Report providing a protocol for assessing the improvement in the quality in the demonstrators",
                "REPORT", "Public", "Jun-2025"),
    Deliverable("D6.2", "WP6", "Evaluation and exploitation of demonstrator results for future C3S implementations",
                "Report providing feedback on the impact of time varying vegetation in reanalysis on seasonal forecasts",
                "REPORT", "Public", "Dec-2025"),
    Deliverable("D6.3", "WP6", "Evaluation and exploitation of demonstrator results for future C3S implementations",
                "Report providing feedback on the impact of land assimilation in reanalysis on seasonal forecasts",
                "REPORT", "Public", "Dec-2025"),
    Deliverable("D6.4", "WP6", "Evaluation and exploitation of demonstrator results for future C3S implementations",
                "Recommendations report on next C3S system developments and future reanalysis and seasonal forecast assessment",
                "REPORT", "Public", "Dec-2026"),
    Deliverable("D6.5", "WP6", "Evaluation and exploitation of demonstrator results for future C3S implementations",
                "Recommendations report of new land products for C3S",
                "REPORT", "Public", "Dec-2026"),

    Deliverable("D7.1", "WP7", "Land-atmosphere interface observations",
                "Satellite snow and soil moisture datasets in the CERISE verification database",
                "DATA", "Confidential", "Dec-2024"),
    Deliverable("D7.2", "WP7", "Land-atmosphere interface observations",
                "Albedo, vegetation and LST satellite datasets in the CERISE verification database",
                "DATA", "Confidential", "Dec-2024"),
    Deliverable("D7.3", "WP7", "Land-atmosphere interface observations",
                "Preliminary in situ observations dataset in the CERISE verification database",
                "DATA", "Confidential", "Dec-2023"),
    Deliverable("D7.4", "WP7", "Land-atmosphere interface observations",
                "Time-varying lake cover, Land Cover and LAI, and extension of CONFESS vegetation data back to 1925 datasets",
                "DATA", "Public", "Apr-2025"),

    Deliverable("D8.1", "WP8", "Project management, coordination, dissemination and oversight",
                "Risk and Quality Management Plan",
                "REPORT", "Confidential", "Mar-2023"),
    Deliverable("D8.2", "WP8", "Project management, coordination, dissemination and oversight",
                "Dissemination and Exploitation Plan",
                "REPORT", "Public", "May-2023"),
    Deliverable("D8.3", "WP8", "Project management, coordination, dissemination and oversight",
                "Media and Communication Plan",
                "REPORT", "Public", "Jun-2023"),
    Deliverable("D8.4", "WP8", "Project management, coordination, dissemination and oversight",
                "Data Management Plan",
                "REPORT", "Public", "Jun-2023"),
    Deliverable("D8.5", "WP8", "Project management, coordination, dissemination and oversight",
                "Interim technical review report",
                "REPORT", "Confidential", "Dec-2024"),
    Deliverable("D8.6", "WP8", "Project management, coordination, dissemination and oversight",
                "Mid-Term Dissemination and Exploitation Report",
                "REPORT", "Public", "Dec-2024"),
    Deliverable("D8.7", "WP8", "Project management, coordination, dissemination and oversight",
                "Final Dissemination and Exploitation Report",
                "REPORT", "Public", "Dec-2026"),
]


def slugify(text: str, max_len: int = 95) -> str:
    text = html.unescape(text)
    text = text.replace("&", "and")
    text = re.sub(r"[^\w\s.-]+", "", text, flags=re.UNICODE)
    text = re.sub(r"\s+", "_", text.strip())
    text = re.sub(r"_+", "_", text)
    return text[:max_len].strip("_")


def wp_folder(deliv: Deliverable) -> str:
    return WP_FOLDER_NAMES.get(deliv.wp, f"{deliv.wp}_{slugify(deliv.wp_title, 70)}")


def deliverable_folder(deliv: Deliverable) -> str:
    return f"{deliv.did}_{slugify(deliv.title, 75)}"


def is_doc_url(url: str) -> bool:
    u = url.lower()
    return (
        "/sites/default/files/" in u
        or u.endswith((".pdf", ".doc", ".docx", ".xls", ".xlsx", ".zip"))
    )


def guess_ext(url: str, content_type: str | None = None) -> str:
    suffix = Path(urlparse(url).path).suffix
    if suffix:
        return suffix

    ct = (content_type or "").lower()
    if "pdf" in ct:
        return ".pdf"
    if "word" in ct or "officedocument.wordprocessingml" in ct:
        return ".docx"
    if "spreadsheet" in ct or "excel" in ct:
        return ".xlsx"
    if "zip" in ct:
        return ".zip"
    return ".bin"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def expected_year(deliv: Deliverable) -> str:
    match = re.search(r"(\d{4})", deliv.est_date)
    return match.group(1) if match else "unknown_year"


def discover_version(url: str) -> str:
    """Return a compact version string from a URL/filename, e.g. v1.1."""
    path = Path(urlparse(url).path).name
    match = re.search(r"(?i)(?:^|[-_.])v(?:ersion)?[-_.]?(\d+(?:[._]\d+)*)", path)
    if not match:
        return ""
    return "v" + match.group(1).replace("_", ".")


def looks_like_html(data: bytes) -> bool:
    head = data[:512].lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html")


def file_metadata(path: Path) -> tuple[str, str]:
    data = path.read_bytes()
    return str(len(data)), sha256_bytes(data)


def discover_links_from_page(html_text: str) -> dict[str, list[str]]:
    """
    Return {deliverable_id: [url, ...]} discovered from the CERISE deliverables page.

    This is deliberately defensive because Drupal pages can change markup.
    Strategy:
      1. Look at all anchors whose href points to likely files.
      2. Try to associate each file with a D-number from nearby row/div text.
      3. Also infer D-number from filename if it contains patterns like D1-2 or D1.2.
    """
    soup = BeautifulSoup(html_text, "html.parser")
    out: dict[str, list[str]] = {}

    for a in soup.find_all("a", href=True):
        href = a["href"]
        url = urljoin(BASE_URL, href)
        if not is_doc_url(url):
            continue

        did: str | None = None

        # Prefer nearby visible container text.
        for parent_name in ["tr", "li", "div", "section", "article"]:
            parent = a.find_parent(parent_name)
            if parent:
                text = " ".join(parent.get_text(" ", strip=True).split())
                m = re.search(r"\bD(\d+)[.-](\d+)\b", text, flags=re.IGNORECASE)
                if m:
                    did = f"D{m.group(1)}.{m.group(2)}"
                    break

        # Fallback to URL/file name: CERISE-D1-2-V1.1.pdf => D1.2
        if did is None:
            m = re.search(r"\bD(\d+)[._-](\d+)\b", url, flags=re.IGNORECASE)
            if m:
                did = f"D{m.group(1)}.{m.group(2)}"

        if did:
            out.setdefault(did, [])
            if url not in out[did]:
                out[did].append(url)

    return out


def bib_key(deliv: Deliverable) -> str:
    # These are project reports, so use CERISE + deliverable ID + expected year.
    year = re.search(r"(\d{4})", deliv.est_date)
    year_str = year.group(1) if year else "nd"
    return f"CERISE_{deliv.did.replace('.', '')}_{year_str}"


def bibtex_entry(deliv: Deliverable, url: str | None) -> str:
    year_match = re.search(r"(\d{4})", deliv.est_date)
    year = year_match.group(1) if year_match else "n.d."
    month = deliv.est_date.split("-")[0] if "-" in deliv.est_date else ""

    fields = [
        f"  title       = {{{deliv.title}}}",
        "  institution = {{CERISE Project}}",
        f"  number      = {{{deliv.did}}}",
        f"  year        = {{{year}}}",
    ]
    if month:
        fields.append(f"  month       = {{{month}}}")
    if url:
        fields.append(f"  url         = {{{url}}}")
    fields.append("  note        = {{Copernicus Climate Change Service Evolution (CERISE) deliverable}}")

    return "@techreport{" + bib_key(deliv) + ",\n" + ",\n".join(fields) + "\n}\n"


def write_readme(outdir: Path, downloaded: int, no_link: int, skipped: int) -> None:
    readme = f"""# CERISE public deliverables archive

Source page: {DELIVERABLES_URL}

Downloaded on: {date.today().isoformat()}

## Contents

- `cerise_deliverables_manifest.csv` — one row per CERISE deliverable, including confidential/skipped items.
- `cerise_deliverables.bib` — starter BibTeX entries for public deliverables.
- `download_log.txt` — run log.
- `../WP*/D*/` — downloaded files organized by work package and deliverable.

## Notes

CERISE deliverables are treated here as project reports / grey literature, not DOI-bearing journal papers.

Some public deliverables have future expected delivery dates or may not have a posted file yet. These remain in the manifest with status `no_file_link_found`.

## This run

- Downloaded public files: {downloaded}
- Public deliverables with no file link found: {no_link}
- Confidential/skipped deliverables: {skipped}
"""
    (outdir / "00_manifest" / "README.md").write_text(readme, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="CERISE", help="Output directory")
    parser.add_argument("--dry-run", action="store_true", help="Do not download files; just write manifests and planned paths")
    parser.add_argument("--overwrite", action="store_true", help="Redownload files even if the local file already exists")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds")
    args = parser.parse_args()

    outdir = Path(args.out)
    manifest_dir = outdir / "00_manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 CERISE public deliverables archiver"
    })

    log_lines: list[str] = []
    log_lines.append(f"Source: {DELIVERABLES_URL}")
    log_lines.append(f"Date: {date.today().isoformat()}")
    log_lines.append(f"Dry run: {args.dry_run}")
    log_lines.append(f"Overwrite: {args.overwrite}")
    log_lines.append("")

    print(f"Reading {DELIVERABLES_URL}")
    page = session.get(DELIVERABLES_URL, timeout=args.timeout)
    page.raise_for_status()

    discovered = discover_links_from_page(page.text)
    log_lines.append(f"Discovered document links for {len(discovered)} deliverable IDs from page markup.")
    for did, urls in sorted(discovered.items()):
        log_lines.append(f"  {did}: {len(urls)} link(s)")
    log_lines.append("")

    manifest_rows = []
    bib_entries = []

    downloaded = 0
    no_link = 0
    skipped = 0

    for deliv in DELIVERABLES:
        base = {
            "deliverable_id": deliv.did,
            "work_package": deliv.wp,
            "wp_title": deliv.wp_title,
            "title": deliv.title,
            "type": deliv.dtype,
            "dissemination": deliv.dissemination,
            "expected_date": deliv.est_date,
            "source_page": DELIVERABLES_URL,
            "file_url": "",
            "local_file": "",
            "download_status": "",
            "download_date": date.today().isoformat(),
            "byte_size": "",
            "sha256": "",
        }

        if deliv.dissemination.lower() != "public":
            base["download_status"] = "skipped_confidential"
            manifest_rows.append(base)
            skipped += 1
            continue

        urls = discovered.get(deliv.did, [])
        if not urls:
            base["download_status"] = "no_file_link_found"
            manifest_rows.append(base)
            bib_entries.append(bibtex_entry(deliv, None))
            no_link += 1
            continue

        # Usually one PDF per deliverable; handle multiple attached docs if present.
        for idx, url in enumerate(urls, start=1):
            row = dict(base)
            row["file_url"] = url

            ext = guess_ext(url)
            version = discover_version(url)
            filename_parts = [
                "CERISE",
                deliv.did,
                expected_year(deliv),
                deliv.est_date,
                slugify(deliv.title),
            ]
            if version:
                filename_parts.append(version)
            if len(urls) > 1:
                filename_parts.append(f"part{idx}")
            filename = "_".join(filename_parts) + ext

            folder = outdir / wp_folder(deliv) / deliverable_folder(deliv)
            local_path = folder / filename
            row["local_file"] = str(local_path)

            if args.dry_run:
                row["download_status"] = "planned"
                manifest_rows.append(row)
                continue

            folder.mkdir(parents=True, exist_ok=True)

            try:
                if local_path.exists() and not args.overwrite:
                    row["download_status"] = "exists"
                    row["byte_size"], row["sha256"] = file_metadata(local_path)
                    downloaded += 1
                    manifest_rows.append(row)
                    continue

                print(f"Downloading {deliv.did}: {url}")
                resp = session.get(url, timeout=args.timeout)
                resp.raise_for_status()
                data = resp.content
                if looks_like_html(data):
                    raise ValueError("response looks like HTML, not a deliverable file")

                # If extension wasn't obvious, update after response.
                ext2 = guess_ext(url, resp.headers.get("Content-Type"))
                if local_path.suffix == ".bin" and ext2 != ".bin":
                    local_path = local_path.with_suffix(ext2)
                    row["local_file"] = str(local_path)
                    if local_path.exists() and not args.overwrite:
                        row["download_status"] = "exists"
                        row["byte_size"], row["sha256"] = file_metadata(local_path)
                        downloaded += 1
                        manifest_rows.append(row)
                        continue

                local_path.write_bytes(data)
                row["download_status"] = "downloaded"
                row["byte_size"] = str(len(data))
                row["sha256"] = sha256_bytes(data)
                downloaded += 1

            except Exception as exc:
                row["download_status"] = f"error: {exc}"

            manifest_rows.append(row)

        # Use first URL in starter BibTeX.
        bib_entries.append(bibtex_entry(deliv, urls[0]))

    manifest_path = manifest_dir / "cerise_deliverables_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "deliverable_id", "work_package", "wp_title", "title", "type",
            "dissemination", "expected_date", "source_page", "file_url",
            "local_file", "download_status", "download_date", "byte_size",
            "sha256",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    bib_path = manifest_dir / "cerise_deliverables.bib"
    bib_path.write_text("\n".join(bib_entries), encoding="utf-8")

    write_readme(outdir, downloaded=downloaded, no_link=no_link, skipped=skipped)

    log_lines.append(f"Manifest: {manifest_path}")
    log_lines.append(f"BibTeX: {bib_path}")
    log_lines.append(f"Downloaded: {downloaded}")
    log_lines.append(f"No file link found: {no_link}")
    log_lines.append(f"Skipped confidential: {skipped}")

    (manifest_dir / "download_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    print()
    print(f"Manifest written: {manifest_path}")
    print(f"BibTeX written:   {bib_path}")
    print(f"README written:   {manifest_dir / 'README.md'}")
    print(f"Downloaded:       {downloaded}")
    print(f"No file links:    {no_link}")
    print(f"Skipped conf.:    {skipped}")

    if args.dry_run:
        print("\nDry run only. Re-run without --dry-run to download files.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
