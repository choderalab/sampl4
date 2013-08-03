#!/usr/local/bin/env python

# Submit jobs to cluster using slurm.

import os
import os.path
import commands

pbs_template = """\
#!/bin/tcsh
#  Batch script for mpirun job on cbio cluster.
#
#SBATCH --job-name=%(jobname)s
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --output=%(jobname)s-%%j.stdout
#SBATCH --error=%(jobname)s-%%j.stderr
#SBATCH --gres=gpu:4

setenv | grep SLURM

canopy
openmm

date

setenv YANKHOME /cbio/jclab/home/chodera/yank

# Select job directory.
setenv JOBDIR %(jobdir)s

# Set YANK directory.
setenv YANKDIR ${YANKHOME}/src/

date

# Change to job directory
cd $JOBDIR

# Clean up old working files
#rm -f *.nc

# Set PYTHONPATH
setenv PYTHONPATH ${YANKDIR}:${PYTHONPATH}

# Run YANK
mpirun -bootstrap slurm python $YANKDIR/yank.py --mpi --receptor_prmtop receptor.prmtop --ligand_prmtop ligand.prmtop --complex_prmtop complex.prmtop --complex_crd complex.crd --output . --restraints flat-bottom --randomize_ligand --iterations 10000 --verbose

date
"""

def isletter(c):
   if ((c >= 'a') and (c <= 'z')) or ((c >= 'A') and (c <= 'Z')):
      return True
   return False   

import os, os.path
molecules = os.listdir('systems')

for molecule in molecules:
   print molecule

   jobname = molecule
   
   jobdir = os.path.join('systems', molecule)
   jobdir = os.path.abspath(jobdir)

   # Form PBS script
   pbs = pbs_template % vars()
   #print pbs

   # Construct directory.
   filename = os.path.join(jobdir, 'run.slurm')
   outfile = open(filename, 'w')
   outfile.write(pbs)
   outfile.close()

   # Submit to PBS
   output = commands.getoutput('sbatch %(filename)s' % vars());
   print output

