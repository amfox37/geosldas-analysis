# IV / TC Python Tools

Early Python port of the independent validation (IV), triple collocation (TC),
and R-diff workflow used for GEOSldas soil-moisture evaluation.

This project is intentionally small and plain at the start. The first useful
piece is a Discover fixture collector that can dry-run the known data locations
and, once checked, copy a few representative files into `test_data/inputs/`.

## Layout

- `iv_tc/`: importable Python package.
- `scripts/collect_discover_fixtures.py`: dry-run/copy helper for small Discover
  fixtures.
- `tests/`: local tests using fake files, so path logic can be checked without
  Discover access.

## Fixture Dry Run

From the repository root:

```bash
python projects/iv_tc/scripts/collect_discover_fixtures.py \
  --date 2020-05-15 \
  --run OLv7_M36_MULTI_type_13_H121=/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121
```

The default mode is dry-run: it prints found/missing source files and the
portable destination path it would use. Add `--copy` only after reviewing the
dry-run output.

