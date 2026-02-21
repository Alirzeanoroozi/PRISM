# Hotspot
In this protocol, we use the HotPoint web server48 to predict hot spots in template protein interfaces, which is freely available to all users. HotPoint labels a residue as a hot spot if it is buried (its relative accessibility is less than 20%, as determined by Naccess) and its total contact potential (with respect to its neighbors within a radius of 7.0 Å) is large (more than 18.0).

# NAccess Download
https://www.bioinf.manchester.ac.uk/naccess/nac_intro.html

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
