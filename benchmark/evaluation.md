# Benchmark
downloaded from [zlab.wenglab.org/benchmark](https://zlab.wenglab.org/benchmark/) and  
![Benchmark image](image.png)


python irmsd.py first_structure.pdb chainA chainB second_structure.pdb
chainA chainB

ChainA's refer to the receptor part of the structures.
ChainB's refer to the ligand part of the structures.

make sure you also run the reciprocal cases:

python irmsd.py first_structure.pdb chainA chainB second_structure.pdb
chainB chainA

There are two scripts in the link, one of them uses the backbone atoms
and the other uses all.

If the resulting iRMSD is greater than 4 we consider two structures to
be different from each other. Reference for this threshold can be found
here or Table III of the attached paper.