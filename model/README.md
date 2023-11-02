# Macaque auditory thalamocortical circuits model (/model folder)
## Description
Multiscale model of the macaque auditory thalamocortical circuits developed using NetPyNE (www.netpyne.org).

The model is described in the following publication:

Dura-Bernal S, Griffith EY, Barczak A, O’Connell MN, McGinnis T, Moreira JV, Schroeder C, Lytton WW, Lakatos P, Neymotin SA. (2023) **Data-driven multiscale model of macaque auditory thalamocortical circuits reproduces in vivo dynamics.** Cell Reports (In Press)



## Setup and execution
This folder includes all the source code of the NetPyNE model of macaque thalamocortical circuits. 
It requires NEURON with Python and MPI support, and NetPyNE. See here installation instructions: http://netpyne.org/install.html#  

The first required step is compiling the NMODL files. This requires running the following command:
From /sim run `nrnivmodl ../mod`. This should create a directory called x86_64. 

The next step requires downloading the model data, see instructions in data/README.md.

To run the full model version, run the following command: `./runsim [num_proc]` or the equivalent `mpiexec -np [num_proc] nrniv -python -mpi init.py`

To run a a batch of simulations via SLURM on high performance computing (HPC) systems use: `python batch.py`

To run a demo version use `python init_demo.py`. This will run a reduced toy auditory thalamocortical model with 5% cell density and simulation duration of 500ms. It will generate an output spiking raster plot and LFP plot figures. Expected run time for demo on a "normal" desktop computer is 150 seconds.

## Overview of file structure:

* /init.py: Main executable to create, simulate and analyze the full model.

* /netParams.py: Network model parameters following the NetPyNE specification. For details see: http://netpyne.org

* /cfg.py: Simulation configuration following the NetPyNE specification. For details see: http://netpyne.org

* /init_demo.py: Demo executable to create, simulate and analyze a toy model with 5% cell density.

* /cfg_demo.py: Simulation configuration for the toy demo model.

* /cells: source .py, .hoc, .json or .pkl files for the different cell types used in the model; these will be imported into NetPyNE.

* /conns: source .py, .json or .pkl files for the network connectivity; these will be imported into NetPyNE.

* /mod: NMODL files containing the ionic channel and synaptic mechanisms used in the model.


For further information please contact: salvador.dura-bernal@downstate.edu, erica.griffith@downstate.edu and/or samuel.neymotin@nki.rfmh.org

