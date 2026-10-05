#! /bin/sh

../src/pops --pdb 1aki.pdb --traj 1aki.sdtraj.gro.2 --compositionOut --typeOut --topologyOut --atomOut --residueOut --chainOut || exit 1

