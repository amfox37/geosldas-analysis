from datetime import date
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from iv_tc.config import ProductRoots, RunConfig, parse_run_spec  # noqa: E402
from iv_tc.fixtures import collect_fixture_files, parse_date  # noqa: E402


def touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x")


def test_parse_date_accepts_two_common_formats():
    assert parse_date("2020-05-15") == date(2020, 5, 15)
    assert parse_date("20200515") == date(2020, 5, 15)


def test_parse_run_spec_requires_name_and_path():
    run = parse_run_spec("OL=/tmp/ol")
    assert run.name == "OL"
    assert run.root == Path("/tmp/ol")


def test_collect_fixture_files_uses_known_paths(tmp_path):
    source = tmp_path / "source"
    out = tmp_path / "out"
    day = date(2020, 5, 15)

    roots = ProductRoots(
        ascat_cdr_root=source / "ASCAT_SSM_CDR",
        legacy_ascat_root=source / "ASCAT_HSAF",
        smosic_root=source / "SMOS_IC" / "preprocessed_m36_daily",
        smap_l3_root=source / "SPL3SMP_v009",
    )
    run = RunConfig("RUN_A", source / "runs" / "RUN_A")

    manifest = (
        roots.ascat_cdr_root
        / "flists"
        / "Y2020"
        / "M05"
        / "D15"
        / "H121_METOPA.txt"
    )
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text("h121_a_1.nc\nh121_a_2.nc\n")
    touch(roots.ascat_cdr_root / "H121" / "metop_a" / "Y2020" / "M05" / "h121_a_1.nc")
    touch(roots.ascat_cdr_root / "H121" / "metop_a" / "Y2020" / "M05" / "h121_a_2.nc")

    touch(
        roots.legacy_ascat_root
        / "H119_H120_processed"
        / "Y2020"
        / "M05"
        / "ASCAT_HSAF_H119_SM_20200515_AD.mat"
    )
    touch(roots.smosic_root / "smos_ic_sm_m36_20200515.nc")
    touch(roots.smap_l3_root / "Y2020" / "SMAP_L3_SM_P_20200515_R19240_001.h5")
    touch(
        run.root
        / "output"
        / roots.domain
        / "cat"
        / "ens_avg"
        / "Y2020"
        / "M05"
        / "RUN_A.tavg24_1d_lnd_Nt.20200515_1200z.nc4"
    )
    touch(run.root / "output" / roots.domain / "rc_out" / "RUN_A.ldas_tilecoord.bin")

    files = collect_fixture_files([day], out, roots, [run])
    labels = [item.label for item in files]
    sources = {item.source.name for item in files}
    dests = {str(item.dest.relative_to(out)) for item in files}

    assert "H121 manifest metop_a 2020-05-15" in labels
    assert "h121_a_1.nc" in sources
    assert "ASCAT_HSAF_H119_SM_20200515_AD.mat" in sources
    assert "smos_ic_sm_m36_20200515.nc" in sources
    assert "SMAP_L3_SM_P_20200515_R19240_001.h5" in sources
    assert "RUN_A.tavg24_1d_lnd_Nt.20200515_1200z.nc4" in sources
    assert "RUN_A.ldas_tilecoord.bin" in sources

    assert "ASCAT_SSM_CDR/flists/Y2020/M05/D15/H121_METOPA.txt" in dests
    assert "ASCAT_SSM_CDR/H121/metop_a/Y2020/M05/h121_a_1.nc" in dests
    assert "SPL3SMP_v009/Y2020/SMAP_L3_SM_P_20200515_R19240_001.h5" in dests
    assert (
        "runs/RUN_A/output/SMAP_EASEv2_M36_GLOBAL/rc_out/RUN_A.ldas_tilecoord.bin"
        in dests
    )


def test_missing_h121_manifest_does_not_invent_raw_files(tmp_path):
    roots = ProductRoots(
        ascat_cdr_root=tmp_path / "ASCAT_SSM_CDR",
        legacy_ascat_root=tmp_path / "ASCAT_HSAF",
        smosic_root=tmp_path / "SMOS_IC",
        smap_l3_root=tmp_path / "SPL3SMP_v009",
    )
    files = collect_fixture_files([date(2020, 5, 15)], tmp_path / "out", roots, [])

    h121_files = [item for item in files if item.label.startswith("H121")]
    assert len(h121_files) == 3
    assert all("manifest" in item.label for item in h121_files)
