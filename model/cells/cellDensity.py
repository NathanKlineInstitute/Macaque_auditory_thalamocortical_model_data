#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 13:22:51 2019

@author: ricky
"""

'''


Script to obtain cell densities for different layers and cell types 
BASED ON M1 model (mouse) [cells/cellDensity.py]

'''

import numpy as np
from scipy.io import loadmat, savemat
from pprint import pprint
from scipy import interpolate
from pylab import *
from pprint import pprint
from netpyne import specs
import pickle

# --------------------------------------------------------------------------------------------- #
# MAIN SCRIPT
# --------------------------------------------------------------------------------------------- #

density = {}

## cell types
cellTypes = ['IT', 'PT', 'CT', 'PV', 'SOM', 'VIP', 'nonVIP'] # added VIP and non-VIP interneuron classes # based on Tremblay et al., 2016 ('GABAergic Interneurons in the Neocortex: From Cellular Properties to Circuits')
## MAY NEED TO CHANGE nonVIP --> NGF (at least for L1?) 
  
# ------------------------------------------------------------------------------------------------------------------
# 1) Use neuron density profile from 3D Quantification paper (Kelly & Hawken 2017) --> neurons/mm3
# Avg for L1, L2/3, L4, L5A, L5B, L6 from fig 6d
# ------------------------------------------------------------------------------------------------------------------
density['nrn_density'] = [0.536*(10**5),2.043*(10**5),1.932*(10**5),3.24*(10**5),2.511*(10**5),1.57*(10**5)] # see: https://docs.google.com/spreadsheets/d/1rXU6ujzg6TBG59XEFuyE1HTJ6VWM-Jr9XIjMMOxvU1g/edit#gid=460197992



# ------------------------------------------------------------------------------------------------------------------
# 2) Percentage of Excitatory Cells [from Lefort09 (mouse S1)]
# Avg for L2/3, L5A, L5B, L6 from fig 2D 
# overall 85:15 ratio consistent with Markram 2015 (87% +- 1% and 13% +- 1%) -->> Again is this percentage or ratio? 
# This is taken from M1 model cellDensity.py & modified 
# -------------------------------------------------------------------------------------------------------------------
#ratioEI = {}
#ratioEI['Lefort09'] = [(0.193+0.11)/2, 0.09, 0.21, 0.21, 0.10]
percentE = {}
percentE['Lefort09'] = [0, 0.88, 0.92, 0.84, 0.83, 0.91] # See Schuman et al., 2019 ("Four Unique Interneuron Populations...") : "Although L1 is devoid of excitatory cells..." 
density[('A1','E')] = [round(density['nrn_density'][i]) * (percentE['Lefort09'][i]) for i in range(len(density['nrn_density']))] # keep in mind this is number of excitatory cells / (mm^3)
density[('A1','I')] = [round(density['nrn_density'][i]) * (1-percentE['Lefort09'][i]) for i in range(len(density['nrn_density']))] 

# ------------------------------------------------------------------------------------------------------------------
# 3) Use interneuron proportions from 'GABAergic interneurons in neocortex' (Tremblay et al., 2016)
# Avg for PV, SOM, VIP, non-VIP in each layer of mouse somatosensory cortex (fig 2)
# ------------------------------------------------------------------------------------------------------------------
PV = 	[0.007, 0.29, 0.641, 0.54, 0.465, 0.424]       # L1: 0.7% (0.007)
SOM = 	[0.04, 0.116, 0.169, 0.319, 0.389, 0.318]      # L1: 4% (0.04)
VIP = 	[0.052, 0.347, 0.092, 0.078, 0.06, 0.064]      # L1: 5.2% (0.052)
nonVIP = [0.9, 0.247, 0.099, 0.064, 0.085, 0.193]      # L1: 90% (0.9)

# Keep in mind again that all of these numbers are /mm3
density[('A1','PV')] =     [(density[('A1','I')][i])*(PV[i]) for i in range(len(PV))] 
density[('A1','SOM')] =    [(density[('A1','I')][i])*(SOM[i]) for i in range(len(SOM))]
density[('A1','VIP')] =    [(density[('A1','I')][i])*(VIP[i]) for i in range(len(VIP))]
density[('A1','nonVIP')] = [(density[('A1','I')][i])*(nonVIP[i]) for i in range(len(nonVIP))]



# save density data in pickle object 
savePickle = 1

data = {'density': density}

if savePickle:
    with open('cellDensity.pkl', 'wb') as fileObj:        
        pickle.dump(data, fileObj)

