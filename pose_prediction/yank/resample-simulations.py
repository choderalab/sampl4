#!/opt/local/bin/python2.5

#=============================================================================================
# Resample configurations from fully-interacting state from equilibrated phase of simulations
# and compile a PDB to do binding mode clustering.
#=============================================================================================

#=============================================================================================
# REQUIREMENTS
#
# The netcdf4-python module is now used to provide netCDF v4 support:
# http://code.google.com/p/netcdf4-python/
#
# This requires NetCDF with version 4 and multithreading support, as well as HDF5.
#=============================================================================================

#=============================================================================================
# TODO
#=============================================================================================

#=============================================================================================
# CHAGELOG
#=============================================================================================

#=============================================================================================
# VERSION CONTROL INFORMATION
#=============================================================================================

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy
from numpy import *
import netCDF4 as netcdf # netcdf4-python
import os
import sys
import os.path
import math
import gzip
from pymbar import MBAR # multistate Bennett acceptance ratio
from pymbar import timeseries # for statistical inefficiency analysis

import simtk.openmm as openmm
import simtk.unit as units

#=============================================================================================
# PARAMETERS
#=============================================================================================

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

#=============================================================================================
# SUBROUTINES
#=============================================================================================
def write_file(filename, contents):
    """Write the specified contents to a file.

    ARGUMENTS
      filename (string) - the file to be written
      contents (string) - the contents of the file to be written

    """

    outfile = open(filename, 'w')

    if type(contents) == list:
        for line in contents:
            outfile.write(line)
    elif type(contents) == str:
        outfile.write(contents)
    else:
        raise "Type for 'contents' not supported: " + repr(type(contents))

    outfile.close()

    return

def read_file(filename):
    """Read contents of the specified file.

    ARGUMENTS
      filename (string) - the name of the file to be read

    RETURNS
      lines (list of strings) - the contents of the file, split by line
    """

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    return lines


def read_pdb(filename):
    """
    Read the contents of a PDB file.

    ARGUMENTS

    filename (string) - name of the file to be read

    RETURNS

    atoms (list of dict) - atoms[index] is a dict of fields for the ATOM residue

    """
    
    # Read the PDB file into memory.
    pdbfile = open(filename, 'r')

    # Extract the ATOM entries.
    # Format described here: http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html
    atoms = list()
    for line in pdbfile:
        if line[0:6] == "ATOM  ":
            # Parse line into fields.
            atom = dict()
            atom["serial"] = line[6:11]
            atom["atom"] = line[12:16]
            atom["altLoc"] = line[16:17]
            atom["resName"] = line[17:20]
            atom["chainID"] = line[21:22]
            atom["Seqno"] = line[22:26]
            atom["iCode"] = line[26:27]
            atom["x"] = line[30:38]
            atom["y"] = line[38:46]
            atom["z"] = line[46:54]
            atom["occupancy"] = line[54:60]
            atom["tempFactor"] = line[60:66]
            atoms.append(atom)
            
    # Close PDB file.
    pdbfile.close()

    # Return dictionary of present residues.
    return atoms

def write_pdb_replica_trajectories(basepdb, directory, prefix, title, ncfile, trajectory_by_state=True):
    """Write out replica trajectories as multi-model PDB files.

    ARGUMENTS
       basepdb (string) - name of PDB file to read atom names and residue information from
       directory (string) - the directory to write files to
       prefix (string) - prefix for replica trajectory files
       title (string) - the title to give each PDB file
       ncfile (NetCDF) - NetCDF file object for input file       
    """
    niterations = ncfile.variables['positions'].shape[0]
    nstates = ncfile.variables['positions'].shape[1]
    natoms = ncfile.variables['positions'].shape[2]
    
    atom_list=read_pdb(basepdb)
    if (len(atom_list) != natoms):
        print ("Number of atoms in trajectory (%d) differs from number of atoms in reference PDB (%d)." % (natoms, len(atom_list)))
        raise Exception
    
    if trajectory_by_state:
    	for state_index in range(0,nstates):
    		print "Working on state %d / %d" % (state_index,nstates)  
        	file_name= "%s-%03d.pdb" % (prefix,state_index)
		full_filename=directory+'/'+file_name
		outfile = open(full_filename, 'w')
		for iteration in range(niterations):
        		state_indices = ncfile.variables['states'][iteration,:]
        		replica_index = list(state_indices).index(state_index)
    			outfile.write('MODEL     %4d\n' % (iteration+1))                        
			write_pdb(atom_list,outfile,iteration,replica_index,title,ncfile,trajectory_by_state=True)
    			outfile.write('ENDMDL\n')
		
		outfile.close()	

    else:
	for replica_index in range(nstates):
		print "Working on replica %d / %d" % (replica_index,nstates)
		file_name="R-%s-%03d.pdb" % (prefix,replica_index)
		full_filename=directory+'/'+file_name
		outfile = open(full_filename, 'w')
		for iteration in range(niterations):
                    outfile.write('MODEL     %4d\n' % (iteration+1))                                            
                    write_pdb(atom_list,outfile,iteration,replica_index,title,ncfile,trajectory_by_state=False)
                    outfile.write('ENDMDL\n')                    

		outfile.close()
		
    return

def write_pdb(atoms, filename, iteration, replica, title, ncfile):
    """Write out replica trajectories as multi-model PDB files.

    ARGUMENTS
       atoms (list of dict) - parsed PDB file ATOM entries from read_pdb() - WILL BE CHANGED
       filename (string) - name of PDB file to be written
       title (string) - the title to give each PDB file
       ncfile (NetCDF) - NetCDF file object for input file       
    """

    # Extract coordinates to be written.
    coordinates = numpy.array(ncfile.variables['positions'][iteration,replica,:,:])
    coordinates *= 10.0 # convert nm to angstroms

    # Create file.
    #outfile = open(filename, 'w')

    # Write ATOM records.
    for (index, atom) in enumerate(atoms):
        atom["x"] = "%8.3f" % coordinates[index,0]
        atom["y"] = "%8.3f" % coordinates[index,1]
        atom["z"] = "%8.3f" % coordinates[index,2]
        filename.write('ATOM  %(serial)5s %(atom)4s%(altLoc)c%(resName)3s %(chainID)c%(Seqno)5s   %(x)8s%(y)8s%(z)8s\n' % atom)
        
    # Close file.
    #outfile.close()

    return

def write_crd(filename, iteration, replica, title, ncfile):
    """
    Write out AMBER format CRD file.

    """
    # Extract coordinates to be written.
    coordinates = array(ncfile.variables['positions'][iteration,replica,:,:])
    coordinates *= 10.0 # convert nm to angstroms

    # Create file.
    outfile = open(filename, 'w')

    # Write title.
    outfile.write(title + '\n')

    # Write number of atoms.
    natoms = ncfile.variables['positions'].shape[2]
    outfile.write('%6d\n' % natoms)

    # Write coordinates.
    for index in range(natoms):
        outfile.write('%12.7f%12.7f%12.7f' % (coordinates[index,0], coordinates[index,1], coordinates[index,2]))
        if ((index+1) % 2 == 0): outfile.write('\n')
        
    # Close file.
    outfile.close()

def _logsum(a_n):
  """
  Compute the log of a sum of exponentiated terms exp(a_n) in a numerically-stable manner:

    _logsum a_n = max_arg + \log \sum_{n=1}^N \exp[a_n - max_arg]

  where max_arg = max_n a_n.  This is mathematically (but not numerically) equivalent to

    _logsum a_n = \log \sum_{n=1}^N \exp[a_n]

  ARGUMENTS
    a_n (numpy array) - a_n[n] is the nth exponential argument
  
  RETURNS
    log_sum (float) - the log of the sum of exponentiated a_n, log (\sum_n exp(a_n))

  EXAMPLE  

  >>> a_n = numpy.array([0.0, 1.0, 1.2], numpy.float64)
  >>> print '%.3e' % _logsum(a_n)
  1.951e+00
    
  """

  # Compute the maximum argument.
  max_log_term = numpy.max(a_n)

  # Compute the reduced terms.
  terms = numpy.exp(a_n - max_log_term)

  # Compute the log sum.
  log_sum = numpy.log(numpy.sum(terms)) + max_log_term
        
  return log_sum

def write_pdb_resampled(ncfile, atoms, output_pdb_filename, reference_state=0, ndiscard=0, nsamples=1000):
    """
    Resample configurations with probability in specified state, concatenating into PDB file.

    ARGUMENTS
       ncfile (NetCDF) - input YANK netcdf file
       atoms (list) - atom records from reference PDB file
       output_pdb_filename (string) - PDB file of resampled configurations to write

    OPTIONAL ARGUMENTS
       reference_state (int) - state to reweight to
       ndiscard (int) - number of iterations to discard to equilibration
       nsamples (int) - number of sampls to generate

    """

    import numpy

    # Get current dimensions.
    niterations = ncfile.variables['energies'].shape[0]
    nstates = ncfile.variables['energies'].shape[1]
    natoms = ncfile.variables['energies'].shape[2]

    # Extract energies.
    print "Reading energies..."
    energies = ncfile.variables['energies']
    u_kln_replica = zeros([nstates, nstates, niterations], float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    print "Done."

    # Deconvolute replicas to collect by thermodynamic state.
    print "Deconvoluting replicas..."
    u_kln = numpy.zeros([nstates, nstates, niterations], float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    print "Done."

    # Discard initial data to equilibration.
    u_kln = u_kln[:,:,ndiscard:]
    [K,L,N] = u_kln.shape
    N_k = N * numpy.ones([K], numpy.int32)

    # Compute snapshot energies.
    print "Computing snapshot energies..."
    u_kn = numpy.zeros([K,N], numpy.float32)
    # Get temperature.
    temperature = ncfile.groups['thermodynamic_states'].variables['temperatures'][reference_state] * units.kelvin
    kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # Boltzmann constant
    kT = kB * temperature # thermal energy
    # Deserialize system.
    system = openmm.XmlSerializer.deserialize(str(ncfile.groups['thermodynamic_states'].variables['systems'][reference_state]))
    # Create Context.
    integrator = openmm.VerletIntegrator(1.0 * units.femtoseconds)
    context = openmm.Context(system, integrator)
    # Turn off restraints.
    context.setParameter('restraint_lambda', 0.0)
    for n in range(N):
        iteration = ndiscard + n
        for replica_index in range(K):        
            # Recompute energies.
            state_index = ncfile.variables['states'][iteration,replica_index]
            positions = ncfile.variables['positions'][iteration,replica_index]
            context.setPositions(positions)
            openmm_state = context.getState(getEnergy=True)
            u_kn[state_index,n] = openmm_state.getPotentialEnergy() / kT

    #===================================================================================================
    # Initialize MBAR.
    #===================================================================================================   
   
    # Initialize MBAR (computing free energy estimates, which may take a while)
    print "Initializing MBAR (warning: using all correlated data)..."
    mbar = MBAR(u_kln, N_k, verbose=False) 

    # Get snapshot weights.
    #u_kn = numpy.squeeze(u_kln[:,state,:])
    log_w_kn = mbar._computeUnnormalizedLogWeights(u_kn)
    f = _logsum(log_w_kn[mbar.indices])
    w_kn = numpy.exp(log_w_kn - f)
    p_kn = w_kn / w_kn[mbar.indices].sum()

    # Form linear list of potential snapshots to choose.
    [K, N] = p_kn.shape
    snapshots = list()
    probabilities = list()
    for replica in range(K):
        for sample in range(N):
            state = ncfile.variables['states'][ndiscard+sample, replica]
            snapshots.append( (replica, ndiscard+sample) )
            probabilities.append( p_kn[state,sample] )
    probabilities = numpy.array(probabilities)

    # Draw samples.
    state_hist = numpy.zeros([K], numpy.int32)
    import numpy.random
    sample_indices = numpy.random.choice(range(len(snapshots)), size=[nsamples], p=probabilities)
    title = 'generated by resampling'
    outfile = open(output_pdb_filename, 'w')
    for (sample_number, sample_index) in enumerate(sample_indices):
        [replica, iteration] = snapshots[sample_index]      
        state = ncfile.variables['states'][iteration, replica]
        state_hist[state] += 1
        print "sample %8d : replica %3d iteration %6d state %3d" % (sample_number, replica, iteration, state)
        outfile.write('MODEL     %4d\n' % (sample_number+1))
        write_pdb(atoms, outfile, iteration, replica, title, ncfile)
        outfile.write('ENDMDL\n')                        
    outfile.close()

    print "Counts from each state:"
    print "%5s %6s" % ('state', 'counts')
    for k in range(K):
        print "%5d %6d" % (k, state_hist[k])
    print ""

    return

def extract_u_n(ncfile):
    """
    Extract timeseries of u_n = - log q(x_n)

    """

    # Get current dimensions.
    niterations = ncfile.variables['energies'].shape[0]
    nstates = ncfile.variables['energies'].shape[1]
    natoms = ncfile.variables['energies'].shape[2]

    # Extract energies.
    print "Reading energies..."
    energies = ncfile.variables['energies']
    u_kln_replica = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    for n in range(niterations):
        u_kln_replica[:,:,n] = energies[n,:,:]
    print "Done."

    # Deconvolute replicas
    print "Deconvoluting replicas..."
    u_kln = numpy.zeros([nstates, nstates, niterations], numpy.float64)
    for iteration in range(niterations):
        state_indices = ncfile.variables['states'][iteration,:]
        u_kln[state_indices,:,iteration] = energies[iteration,:,:]
    print "Done."

    # Compute total negative log probability over all iterations.
    u_n = numpy.zeros([niterations], numpy.float64)
    for iteration in range(niterations):
        u_n[iteration] = numpy.sum(numpy.diagonal(u_kln[:,:,iteration]))

    return u_n

#=============================================================================================
# MAIN
#=============================================================================================

data_directory = 'systems/'

# Store molecule data.
molecule_data = dict()

# Generate list of files in this directory.
import commands
molecules = commands.getoutput('ls -1 %s' % data_directory).split()
for molecule in molecules:
    #try:
        source_directory = os.path.join(data_directory, molecule)
        print source_directory

        # Process complex NetCDF file.
        fullpath = os.path.join(source_directory, 'complex.nc')

        # Skip if the file doesn't exist.
        if (not os.path.exists(fullpath)): continue
    
        # Open NetCDF file for reading.
        print "Opening NetCDF trajectory file '%(fullpath)s' for reading..." % vars()
        ncfile = netcdf.Dataset(fullpath, 'r')

        # DEBUG
        print "dimensions:"
        for dimension_name in ncfile.dimensions.keys():
            print "%16s %8d" % (dimension_name, len(ncfile.dimensions[dimension_name]))
    
        # Read dimensions.
        niterations = ncfile.variables['positions'].shape[0]
        nstates = ncfile.variables['positions'].shape[1]
        natoms = ncfile.variables['positions'].shape[2]
        print "Read %(niterations)d iterations, %(nstates)d states" % vars()

        # Read reference PDB file.
        reference_pdb_filename = os.path.join(source_directory, "complex.pdb")
        atoms = read_pdb(reference_pdb_filename)

        # Choose number of samples to discard to equilibration
        u_n = extract_u_n(ncfile)
        if numpy.any(numpy.isnan(u_n)): continue
        nskip = int(len(u_n) / 100.0)
        [nequil, g_t, Neff_max] = timeseries.detectEquilibration(u_n, nskip)
        print [nequil, Neff_max]

        # Resample configurations for state 0.
        state = 0
        nsamples = 5000
        output_pdb_filename = os.path.join(source_directory, 'resampled.pdb')
        write_pdb_resampled(ncfile, atoms, output_pdb_filename, state, nequil, nsamples)

        # Close input NetCDF file.
        ncfile.close()
    #except Exception as e:
    #    print str(e)
    #    pass
