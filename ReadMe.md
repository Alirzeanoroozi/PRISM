# Steps to run
1. Get the `template.zip` from github and unzip that
2. Get the `pdb.zip` from email


# DB	
$my_host = '172.20.31.92';
$my_user = 'prismUser';
$my_pass = 'prismPass';
$my_db = 'prismDatabaseTMalignRosettaRemote';

# Running commads
``` bash
mkdir jobs
mkdir pdb
cd external_tools
chmod 755 TMalign
cd naccess
csh install.scr 
cd ../BeEM-master
g++ -O3 BeEM.cpp -o BeEM
cd ../..
```

