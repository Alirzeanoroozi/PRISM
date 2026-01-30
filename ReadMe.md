# Hotspot
In this protocol, we use the HotPoint web server48 to predict hot spots in template protein interfaces, which is freely available to all users. HotPoint labels a residue as a hot spot if it is buried (its relative accessibility is less than 20%, as determined by Naccess) and its total contact potential (with respect to its neighbors within a radius of 7.0 Å) is large (more than 18.0).

# STEPS to develpe
1. Read paper to see how template generator works?
2. What is difference between contact and hotspot and interfaces?
3. update templates and pdb

# NAccess Download
https://www.bioinf.manchester.ac.uk/naccess/nac_intro.html

# paper link
https://www.nature.com/articles/nprot.2011.367
https://pmc.ncbi.nlm.nih.gov/articles/PMC7384353/pdf/nihms-1607029.pdf

# Steps to run
1. Get the `template.zip` from github and unzip that
2. Get the `pdb.zip` from email
https://drive.google.com/file/d/1kZfTXfaip0B-J43E6HX8o7DsO_Ckg_7v/view

# DB	
$my_host = '172.20.31.92';
$my_user = 'prismUser';
$my_pass = 'prismPass';
$my_db = 'prismDatabaseTMalignRosettaRemote';

# Running commads
``` bash
cd external_tools
# for mac
g++ -O3 -ffast-math -lm -o TMalign TMalign.cpp
# for Linux
g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp
chmod 755 TMalign
cd naccess
csh install.scr 
cd ../..
module load rosetta/2022.42
python prism.py
```

