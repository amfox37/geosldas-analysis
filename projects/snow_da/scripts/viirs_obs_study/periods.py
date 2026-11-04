"""
Study periods and products (see docs/viirs_obs_study_plan.md; confirmed 2026-11-04).

Shared by all viirs_obs_study scripts so that every step uses the same days.
"""

from datetime import date, timedelta
from pathlib import Path

# 7 consecutive days per period
PERIODS = {
    'P1': ('2020-01-13', '2020-01-19', 'mid-winter; overlaps the CF0360 reader-test runs'),
    'P2': ('2020-04-11', '2020-04-17', 'spring melt'),
    'P3': ('2024-01-13', '2024-01-19', 'mid-winter, J1 + J2 + S-NPP'),
    'P4': ('2024-03-23', '2024-03-29', 'early melt'),
    'P5': ('2024-05-18', '2024-05-24', 'late melt'),
    'P6': ('2024-07-13', '2024-07-19', 'Southern Hemisphere winter'),
    'P7': ('2024-11-09', '2024-11-15', 'early accumulation'),
}

VIIRS = ('VJ110C1', 'VJ210C1', 'VNP10C1')      # JPSS-1, JPSS-2 (from 2023-02-10), S-NPP
MODIS = ('MYD10C1', 'MOD10C1')                 # Aqua (13:30, main reference), Terra (10:30)

JPSS2_FIRST_DAY = date(2023, 2, 10)

# downloaded files: <DATA_DIR>/<product>/<YYYY>/<product>.AYYYYDDD.<version>.<ext> (renamed, as for GEOSldas)
DATA_DIR = Path(__file__).resolve().parents[4] / 'data' / 'viirs_obs_study' / 'cmg'


def days(period):
    """List of datetime.date for a period key."""
    first, last, _ = PERIODS[period]
    d0, d1 = date.fromisoformat(first), date.fromisoformat(last)
    return [d0 + timedelta(n) for n in range((d1 - d0).days + 1)]


def products(period):
    """Products available for a period (VJ210C1 only after JPSS-2 data start)."""
    first = days(period)[0]
    viirs = tuple(p for p in VIIRS if p != 'VJ210C1' or first >= JPSS2_FIRST_DAY)
    return viirs + MODIS


def cmg_file(product, day):
    """Path of the renamed daily file of a product (as written by get_snow_cmg_obs.py)."""
    version, ext = ('002', 'h5') if product.startswith('V') else ('061', 'hdf')
    return DATA_DIR / product / f'{day.year}' / f'{product}.A{day.year}{day.timetuple().tm_yday:03d}.{version}.{ext}'
