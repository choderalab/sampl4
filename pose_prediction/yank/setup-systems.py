#!/bin/env python

import glob
import os, os.path
import shutil
import commands

#===================================================================================================================================
# MAIN
#===================================================================================================================================

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

    # Charge and parameterize ligand.
    import ligandtools
    molecule = ligandtools.readMolecule('ligand.mol2')
    gaff_mol2_filename = 'ligand.gaff.mol2'
    frcmod_filename = 'ligand.frcmod'
    nconfs = 10
    charged_molecule = ligandtools.assignPartialCharges(molecule, charge_model = 'am1bcc', multiconformer = nconfs, minimize_contacts = False, verbose = False)
    ligandtools.parameterizeForAmber(charged_molecule, charge_model=None, verbose=True, resname='MOL', ligand_obj_name='molecule', frcmod_filename=frcmod_filename, gaff_mol2_filename=gaff_mol2_filename)

    # Set up complex.
    leap_template = """\
# Set up complex for GBSA simulation with OBC model.

# Load AMBER 99SB-ILDN forcefield.
source leaprc.ff99SBildn

# Load GAFF parameters.
source leaprc.gaff

# Set GB radii to recommended values for OBC.
set default PBRadii mbondi2 

# Load in protein.
receptor_A = loadPdb ../../../test/3NF8_A.modeller.pdb
receptor_B = loadPdb ../../../test/3NF8_B.modeller.pdb
receptor = combine { receptor_A receptor_B }

# Load parameters for ligand.
loadAmberParams ligand.frcmod

# Load ligand.
ligand = loadMol2 ligand.gaff.mol2

# Create complex.
complex = combine { receptor ligand }

# Check complex.
check complex

# Report on net charge.
charge complex

# Write parameters.
saveAmberParm ligand ligand.prmtop ligand.crd
saveAmberParm receptor receptor.prmtop receptor.crd
saveAmberParm complex complex.prmtop complex.crd

# Exit
quit
    """
    outfile = open('setup.leap.in', 'w')
    outfile.write(leap_template)
    outfile.close()
    command = 'tleap -f setup.leap.in >& setup.leap.out'
    output = commands.getoutput(command)

    # Generate PDB files
    for prefix in ['complex', 'ligand', 'receptor']:
        command = 'cat %s.crd | ambpdb -p %s.prmtop > %s.pdb' % (prefix, prefix, prefix)
        output = commands.getoutput(command)
        print output
    
    # Restore original directory.
    os.chdir(cwd)

