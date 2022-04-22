# Macaque auditory thalamocortical circuits model
## Description
Multiscale model of the macaque auditory thalamocortical circuits developed using NetPyNE (www.netpyne.org).

The model is described in the following preprint paper (currently under review in Neuron):

Dura-Bernal S, Griffith EY, Barczak A, O’Connell MN, McGinnis T, Schroeder C, Lytton WW, Lakatos P, Neymotin SA. (2022) **Data-driven multiscale model of macaque auditory thalamocortical circuits reproduces in vivo dynamics.** BioRxiv 2022.02.03.479036


## Setup and execution

Requires NEURON with Python and MPI support. 

1. From /sim run `nrnivmodl ../mod`. This should create a directory called x86_64. 
2. To run type: `./runsim [num_proc]' or the equivalent `mpiexec -np [num_proc] nrniv -python -mpi init.py`

## Overview of file structure:

* /init.py: Main executable; calls functions from other modules. Sets what parameter file to use.

* /netParams.py: Network parameters

* /cfg.py: Simulation configuration

* /cells: source .py, .hoc, .json or .pkl files for the different cell types used in the model; these will be imported into netpyne

* /conns: source .py, .json or .pkl files for the network connectivity; these will be imported into netpyne

* /mod: NMODL files containing the ionic channel and synaptic mechanisms used in the model 

* /data: where the model and simulation data is stored 


For further information please contact: salvador.dura-bernal@downstate.edu 

