
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
OLD = ROOT / "tmp/agent/20260728-matched-align-refiner-replay/results-v2"
OUT = ROOT / "tmp/agent/20260728-fiberdock-corrected-replay/results"
FILES = {
    "native_template": ROOT / "tmp/agent/20260722-pipeline-validation/tm_external/processed/pdbs/1ahw.pdb",
    "tmalign_transform_left": OLD / "tmalign/processed/transformation/1ahwAF_1fgnHL_1tfhA_o1_L.pdb",
    "tmalign_transform_right": OLD / "tmalign/processed/transformation/1ahwAF_1fgnHL_1tfhA_o1_R.pdb",
    "multiprot_transform_left": OLD / "multiprot/processed/transformation/1ahwBC_1fgnHL_1tfhA_o1_L.pdb",
    "multiprot_transform_right": OLD / "multiprot/processed/transformation/1ahwBC_1fgnHL_1tfhA_o1_R.pdb",
    "tmalign_fiberdock": OUT / "tmalign_1ahwAF_o1/fiberdock_energies_1.ref.pdb",
    "multiprot_fiberdock": OUT / "multiprot_1ahwBC_o1/fiberdock_energies_1.ref.pdb",
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
cmd.show("cartoon", "tmalign_transform_left or tmalign_transform_right")
cmd.color("yellow", "tmalign_transform_left")
cmd.color("orange", "tmalign_transform_right")
cmd.show("cartoon", "multiprot_transform_left or multiprot_transform_right")
cmd.color("green", "multiprot_transform_left")
cmd.color("lime", "multiprot_transform_right")
cmd.show("cartoon", "tmalign_fiberdock")
cmd.color("cyan", "tmalign_fiberdock")
cmd.show("cartoon", "multiprot_fiberdock")
cmd.color("salmon", "multiprot_fiberdock")
cmd.set("cartoon_transparency", 0.70, "native_template")
cmd.set("cartoon_transparency", 0.45, "tmalign_transform_left or tmalign_transform_right or multiprot_transform_left or multiprot_transform_right")
cmd.set("cartoon_transparency", 0.10, "tmalign_fiberdock or multiprot_fiberdock")
cmd.orient("native_template or tmalign_fiberdock or multiprot_fiberdock")
cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")
cmd.png(str(OUT / "corrected_fiberdock_replay.png"), width=1800, height=1300, dpi=150)
cmd.save(str(OUT / "corrected_fiberdock_replay.pse"))
cmd.quit()
