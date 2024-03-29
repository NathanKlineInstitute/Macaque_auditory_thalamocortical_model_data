# Macaque auditory thalamocortical circuits model
## Description
Multiscale model of the macaque auditory thalamocortical circuits developed using NetPyNE (www.netpyne.org).

The model is described in the following publication:

Dura-Bernal S, Griffith EY, Barczak A, O’Connell MN, McGinnis T, Moreira JV, Schroeder C, Lytton WW, Lakatos P, Neymotin SA. (2023) **Data-driven multiscale model of macaque auditory thalamocortical circuits reproduces in vivo dynamics.** Cell Reports (In Press)


## Overview of file structure:

This repository contains all the source code and data required to generate the paper results. Each subfolder has a separate README.md with instructions:

* /model: Python source code of the auditory thalamocortical circuit model. Implemented in NetPyNE/NEURON.

* /analysis: Python source code to perform the analyses used to generate the paper results. 

* /data: Experimental data and data generated by model used to generate the paper results. 

## System requirements
The computational model has the following software dependencies:
- NetPyNE v1.0.2 
- NEURON v8.0
- UR_EAR 2.0 

The experimental and model data analysis has the following software dependencies:
- MATLAB v9.3 (64bit) 2017b
- OEvent v51cc9b

The computational model and analysis code has been tested on Mac OS 12.1 and Linux Ubuntu 16.04

## Installation guide

Detailed installation instructions for the software dependencies are available in the links below:
- NetPyNE: http://netpyne.org/install.html
- NEURON v8.0: https://neuron.yale.edu/neuron/ 
- UR_EAR 2.0: https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/carney-lab/codes/UR_EAR_2020b.zip 

Typical install time on a "normal" desktop computer is 1 hour.

## Demo
Instructions, expected output and expected run time to simulate a demo version of the thalamocortical circuit model can be found in model/README.md.

## Instructions for use
Instructions on how to run the full thalamocortical model and analyses can be found in the model/README.md, analysis/README.md and data/README.md files.


For further information please contact: salvador.dura-bernal@downstate.edu, erica.griffith@downstate.edu and/or samuel.neymotin@nki.rfmh.org

