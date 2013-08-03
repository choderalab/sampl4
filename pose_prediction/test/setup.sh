#!/bin/tcsh

source ~/code/amber/amber.csh

antechamber -i AVX-15988_0.mol2 -fi mol2 -o AVX-15988_0.gaff.mol2 -fo mol2 -c bcc
parmchk -i AVX-15988_0.gaff.mol2 -o AVX-15988_0.frcmod -f mol2
tleap -f setup.leap.in

cat complex.crd | ambpdb -p complex.prmtop > complex.pdb
cat ligand.crd | ambpdb -p ligand.prmtop > ligand.pdb
cat receptor.crd | ambpdb -p receptor.prmtop > receptor.pdb
