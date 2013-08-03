#!/bin/env python

import glob
import os, os.path
import shutil
import commands

# Setup systems for YANK using Antechamber

# Get ligand names.
ligand_directory = '../ligands'
ligand_mol2_filenames = glob.glob(ligand_directory + '/*.mol2')
ligand_names = list()
for filename in ligand_mol2_filenames:
    [head, tail] = os.path.split(filename)
    [root, ext] = os.path.splitext(tail)
    ligand_name = root
    ligand_names.append(ligand_name)
print ligand_names

# Parameterize ligands.
systems_directory = 'systems'
for ligand_name in ligand_names:
    print ligand_name

    ligand_mol2_filename = os.path.join(ligand_directory, ligand_name + '.mol2')
    system_directory = os.path.join(systems_directory, ligand_name)
    if not os.path.exists(system_directory): os.makedirs(system_directory)
    shutil.copyfile(ligand_mol2_filename, os.path.join(system_directory, 'ligand.mol2'))
    
    # Change directory.
    cwd = os.getcwd()
    os.chdir(system_directory)

    # Run antechamber.
    #command = 'antechamber -i ligand.mol2 -fi mol2 -o ligand.gaff.mol2 -fo mol2 -c bcc -nc 0'
    #output = commands.getoutput(command)
    #print output

    # Run parmchk.
    #command = 'parmchk -i ligand.gaff.mol2 -o ligand.frcmod -f mol2'
    #output = commands.getoutput(command)
    #print output

    import ligandtools
    molecule = ligandtools.readMolecule('ligand.mol2')
    gaff_mol2_filename = 'ligand.gaff.mol2'
    frcmod_filename = 'ligand.frcmod'
    ligandtools.parameterizeForAmber(molecule, charge_model='bcc', verbose=True, resname='MOL', ligand_obj_name='molecule', frcmod_filename=frcmod_filename, gaff_mol2_filename=gaff_mol2_filename)
    # Restore original directory.
    os.chdir(cwd)
    
    
