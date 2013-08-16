#!/bin/tcsh

canopy
openmm-git
source /cbio/jclab/share/pymbar/pymbar.csh

setenv YANKDIR ${HOME}/yank/src/
setenv JOBDIR "."

# Simple benzene-toluene test.
#setenv JOBDIR "examples/benzene-toluene"



python analyze.py


