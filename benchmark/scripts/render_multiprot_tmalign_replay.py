# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///

import os
from pathlib import Path

os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol  # pytype: disable=import-error

pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd  # pytype: disable=import-error


ROOT = Path("/scratch/rshadi25/GitHub/PRISM-prescript")
OUTPUT = ROOT / "tmp/agent/20260728-multiprot-native-gate-replay/results"
FILES = {
    "native_template": ROOT / "tmp/agent/20260722-pipeline-validation/tm_external/processed/pdbs/1ahw.pdb",
    "tm_transform_left": ROOT / "tmp/agent/20260722-pipeline-validation/tm_external/processed/transformation/1ahwBC_1ahwB_1ahwC_o1_L.pdb",
    "tm_transform_right": ROOT / "tmp/agent/20260722-pipeline-validation/tm_external/processed/transformation/1ahwBC_1ahwB_1ahwC_o1_R.pdb",
    "mp_transform_left": OUTPUT / "processed/transformation/1ahwBC_1fgnHL_1tfhA_o1_L.pdb",
    "mp_transform_right": OUTPUT / "processed/transformation/1ahwBC_1fgnHL_1tfhA_o1_R.pdb",
    "tm_rosetta": ROOT / "tmp/agent/20260722-pipeline-validation/tm_external/processed/rosetta_refinement/structures/1ahwBC_1ahwB_1ahwC_o1_L_1ahwBC_1ahwB_1ahwC_o1_R_rosetta_0001_0001.pdb",
    "tm_fiberdock": ROOT / "tmp/agent/20260722-pipeline-validation/tm_fiberdock/processed/fiberdock_refinement/1ahwBC_1ahwB_1ahwC_o1/fiberdock_energies_1.ref.pdb",
    "mp_fiberdock": OUTPUT / "processed/fiberdock_refinement/1ahwBC_1fgnHL_1tfhA_o1/fiberdock_energies_1.ref.pdb",
    "mp_rosetta_raw": OUTPUT / "processed/rosetta_refinement/structures/1ahwBC_1fgnHL_1tfhA_o1_L_1ahwBC_1fgnHL_1tfhA_o1_R_rosetta_0001_0001.pdb",
}

for name, path in FILES.items():
    if not path.is_file():
        print(f"ERROR missing {name}: {path}")
        cmd.quit()
        raise SystemExit(2)
    cmd.load(str(path), name)
    atoms = cmd.count_atoms(name)
    print(f"loaded {name} atoms={atoms} path={path}")
    if atoms == 0:
        print(f"ERROR zero atoms for {name}")
        cmd.quit()
        raise SystemExit(2)

cmd.hide("everything", "all")
cmd.show("cartoon", "native_template")
cmd.color("gray70", "native_template")
cmd.show("cartoon", "tm_transform_left or tm_transform_right")
cmd.color("yellow", "tm_transform_left")
cmd.color("orange", "tm_transform_right")
cmd.show("cartoon", "mp_transform_left or mp_transform_right")
cmd.color("green", "mp_transform_left")
cmd.color("lime", "mp_transform_right")
cmd.show("cartoon", "tm_rosetta or tm_fiberdock")
cmd.color("cyan", "tm_rosetta")
cmd.color("magenta", "tm_fiberdock")
cmd.show("cartoon", "mp_fiberdock or mp_rosetta_raw")
cmd.color("salmon", "mp_fiberdock")
cmd.color("purple", "mp_rosetta_raw")
cmd.set("cartoon_transparency", 0.65, "native_template")
cmd.set("cartoon_transparency", 0.45, "tm_transform_left or tm_transform_right or mp_transform_left or mp_transform_right")
cmd.set("cartoon_transparency", 0.1, "tm_rosetta or tm_fiberdock or mp_fiberdock or mp_rosetta_raw")
cmd.orient("native_template or tm_rosetta or tm_fiberdock or mp_fiberdock or mp_rosetta_raw")
cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")
cmd.png(str(OUTPUT / "1ahwBC_multiprot_tmalign_replay.png"), width=1800, height=1300, dpi=150)
cmd.save(str(OUTPUT / "1ahwBC_multiprot_tmalign_replay.pse"))
cmd.quit()
