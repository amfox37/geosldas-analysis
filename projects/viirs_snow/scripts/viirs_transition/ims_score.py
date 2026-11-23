"""Re-run the analysis-C IMS scoring (skill_score.run_ims, unchanged logic) on the flip-corrected per-cell counts from
ims_extract.py. Outputs IMS_skill_by_arm_fixed.csv and IMS_skill_contrasts_fixed.csv (the original C_* files are untouched)."""
import os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import skill_common as sc
src = sc.OUT / "imsfix"; links = sc.OUT / "imsfix_links"; links.mkdir(exist_ok=True)
for a in sc.load_arms()["arms"]:
    l = links / f"skill_{a}.npz"
    if l.is_symlink() or l.exists(): l.unlink()
    l.symlink_to(src / f"ims_{a}.npz")
sc.OUT = links                      # skill_score reads its inputs from sc.OUT at import time
import skill_score as ss
if __name__ == "__main__":
    a, c = ss.run_ims()
    a.to_csv(sc.FIG / "IMS_skill_by_arm_fixed.csv", index=False); c.to_csv(sc.FIG / "IMS_skill_contrasts_fixed.csv", index=False)
    print(a.shape, c.shape)
