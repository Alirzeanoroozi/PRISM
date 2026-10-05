# Hotspot
In this protocol, we use the HotPoint web server48 to predict hot spots in template protein interfaces, which is freely available to all users. HotPoint labels a residue as a hot spot if it is buried (its relative accessibility is less than 20%, as determined by Naccess) and its total contact potential (with respect to its neighbors within a radius of 7.0 Å) is large (more than 18.0).

# NAccess Download
https://www.bioinf.manchester.ac.uk/naccess/nac_intro.html

Surface extraction supports both NACCESS and FreeSASA. NACCESS is the default
for legacy-compatible runs; select FreeSASA explicitly with:

```bash
python prism.py --surface_backend freesasa \
  --freesasa_python /path/to/python-with-freesasa
```

FreeSASA uses the Lee-Richards algorithm, a 1.4 Angstrom probe, and the same
PRISM residue reference table for the surface threshold. Use NACCESS for the
legacy baseline and compare backend outputs before interpreting prediction
count differences.

# Steps to run
1. Get the `templates.zip` from github and unzip that

# Running commads
``` bash
cd external_tools
g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp
chmod 755 TMalign
cd naccess
csh install.scr 
cd ../..
module load rosetta/2022.42
python prism.py --generate_templates True
```

# GTalign (Optional)
PRISM can also run the alignment stage with GTalign (instead of TMalign) using the `--aligner gtalign` option.

## Install GTalign (recommended: Conda, from official GTalign repo instructions)
CPU (multiprocessing) version:
```bash
conda install minmarg::gtalign_mp
```

GPU version:
```bash
conda install minmarg::gtalign_gpu
```

## Example PRISM run with GTalign
```bash
python prism.py \
  --aligner gtalign \
  --gtalign_path /path/to/gtalign \
  --generate_templates True
```

Notes:
- The default PRISM behavior remains `TMalign` (`--aligner tmalign`).
- When `--aligner gtalign` is used, alignment JSON files are written to `processed/alignment_gtalign/` so they are distinguishable from TMalign outputs.
