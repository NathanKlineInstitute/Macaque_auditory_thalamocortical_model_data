"""
cfg.py 

Simulation configuration for A1 model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: ericaygriffith@gmail.com, salvadordura@gmail.com
"""


from netpyne import specs
import pickle

cfg = specs.SimConfig()

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 0.5*1e3			## Duration of the sim, in ms
cfg.dt = 0.05                   ## Internal Integration Time Step  
cfg.verbose = 0         	## Show detailed messages
cfg.hParams['celsius'] = 37
cfg.createNEURONObj = 1
cfg.createPyStruct = 1
cfg.printRunTime = 0.1

cfg.connRandomSecFromList = False  # set to false for reproducibility 
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1
cfg.oneSynPerNetcon = False
cfg.includeParamsLabel = False
cfg.printPopAvgRates = [0, cfg.duration]

cfg.validateNetParams = True

#------------------------------------------------------------------------------
# Recording 
#------------------------------------------------------------------------------
cfg.allpops = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM', 'IC']
cfg.allCorticalPops = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6']
cfg.allThalPops = ['TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM', 'IC']

alltypes = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'ITS4', 'PT5B', 'TC', 'HTC', 'IRE', 'TI']

cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc': 0.5, 'var':'v'}}  ## Dict with traces to record 
cfg.recordStim = False			
cfg.recordTime = False  		 
cfg.recordStep = 0.1            ## Step size (in ms) to save data  

cfg.recordLFP = [[100, y, 100] for y in range(0, 2000, 100)]
cfg.recordLFP = [[x, 1000, 100] for x in range(100, 2200, 200)]
#cfg.saveLFPPops =  cfg.allCorticalPops #, "IT3", "SOM3", "PV3", "VIP3", "NGF3", "ITP4", "ITS4", "IT5A", "CT5A", "IT5B", "PT5B", "CT5B", "IT6", "CT6"]

# cfg.recordDipole = True
# cfg.saveDipoleCells = ['all']
# cfg.saveDipolePops = cfg.allpops


#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------

cfg.simLabel = 'demo' 
cfg.saveFolder = '.'                	## Set file output name
cfg.savePickle = True         	## Save pkl file
cfg.saveJson = False           	## Save json file
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net'] 
cfg.backupCfgFile = None 		
cfg.gatherOnlySimData = False	 
cfg.saveCellSecs = True		 
cfg.saveCellConns = False		 

#------------------------------------------------------------------------------
# Analysis and plotting 
#----------------------------------------------------------------------------- 
#

# cfg.analysis['plotTraces'] = {'include': [(pop, 0) for pop in cfg.allpops], 'oneFigPer': 'trace', 'overlay': True, 'saveFig': True, 'showFig': False, 'figSize':(12,8)} #[(pop,0) for pop in alltypes]		## Seen in M1 cfg.py (line 68) 
cfg.analysis['plotRaster'] = {'include': cfg.allpops, 'saveFig': True, 'showFig': False, 'popRates': True, 'orderInverse': True, 'timeRange': [0,cfg.duration], 'figSize': (14,12), 'lw': 0.3, 'markerSize': 3, 'marker': '.', 'dpi': 300}      	## Plot a raster
# cfg.analysis['plotSpikeStats'] = {'stats': ['rate'], 'figSize': (6,12), 'timeRange': [0, 2500], 'dpi': 300, 'showFig': 0, 'saveFig': 1}

cfg.analysis['plotLFP'] = {'plots': ['timeSeries'], 'electrodes': [10], 'maxFreq': 80, 'figSize': (8,4), 'saveData': False, 'saveFig': True, 'showFig': False} # 'PSD', 'spectrogram'
#cfg.analysis['plotDipole'] = {'saveFig': True}
#cfg.analysis['plotEEG'] = {'saveFig': True}


#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
cfg.weightNormThreshold = 5.0  # maximum weight normalization factor with respect to the soma
cfg.weightNormScaling = {'NGF_reduced': 1.0, 'ITS4_reduced': 1.0}


#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------
cfg.AMPATau2Factor = 1.0
cfg.synWeightFractionEE = [0.5, 0.5] # E->E AMPA to NMDA ratio
cfg.synWeightFractionEI = [0.5, 0.5] # E->I AMPA to NMDA ratio
cfg.synWeightFractionSOME = [0.9, 0.1] # SOM -> E GABAASlow to GABAB ratio
cfg.synWeightFractionNGF = [0.5, 0.5] # NGF GABAA to GABAB ratio
cfg.synWeightFractionENGF = [0.834, 0.166] # NGF AMPA to NMDA ratio


#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
cfg.singleCellPops = False
cfg.singlePop = ''
cfg.removeWeightNorm = False
cfg.scale = 1.0      
cfg.sizeY = 2000.0 
cfg.sizeX = 200.0   
cfg.sizeZ = 200.0
cfg.scaleDensity = 0.05 # Should be 1.0 unless need lower cell density for test simulation or visualization


#------------------------------------------------------------------------------
# Connectivity
#------------------------------------------------------------------------------
cfg.synWeightFractionEE = [0.5, 0.5] # E->E AMPA to NMDA ratio
cfg.synWeightFractionEI = [0.5, 0.5] # E->I AMPA to NMDA ratio
cfg.synWeightFractionIE = [0.9, 0.1]  # SOM -> E GABAASlow to GABAB ratio (update this)
cfg.synWeightFractionII = [0.9, 0.1]  # SOM -> E GABAASlow to GABAB ratio (update this)

# Cortical
cfg.addConn = 1

cfg.EEGain = 1.0
cfg.EIGain = 1.0 	
cfg.IEGain = 1.0 
cfg.IIGain = 1.0 

## E/I->E/I layer weights (L1-3, L4, L5, L6)
cfg.EELayerGain = {'1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0 , '5A': 1.0, '5B': 1.0, '6': 1.0}
cfg.EILayerGain = {'1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0 , '5A': 1.0, '5B': 1.0, '6': 1.0}
cfg.IELayerGain = {'1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0 , '5A': 1.0, '5B': 1.0, '6': 1.0}
cfg.IILayerGain = {'1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0 , '5A': 1.0, '5B': 1.0, '6': 1.0}

## E->I by target cell type
cfg.EICellTypeGain= {'PV': 1.0, 'SOM': 1.0, 'VIP': 1.0, 'NGF': 1.0}

## I->E by target cell type
cfg.IECellTypeGain= {'PV': 1.0, 'SOM': 1.0, 'VIP': 1.0, 'NGF': 1.0}

# Thalamic
cfg.addIntraThalamicConn = 1.0
cfg.addIntraThalamicConn = 1.0
cfg.addCorticoThalamicConn = 1.0
cfg.addThalamoCorticalConn = 1.0

cfg.thalamoCorticalGain = 1.0
cfg.intraThalamicGain = 1.0
cfg.corticoThalamicGain = 1.0

cfg.addSubConn = 1

## full weight conn matrix
with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
cfg.wmat = connData['wmat']


#------------------------------------------------------------------------------
# Background inputs
#------------------------------------------------------------------------------
cfg.addBkgConn = 1
cfg.noiseBkg = 1.0  # firing rate random noise
cfg.delayBkg = 5.0  # (ms)
cfg.startBkg = 0  # start at 0 ms

# cfg.weightBkg = {'IT': 12.0, 'ITS4': 0.7, 'PT': 14.0, 'CT': 14.0,
#                 'PV': 28.0, 'SOM': 5.0, 'NGF': 80.0, 'VIP': 9.0,
#                 'TC': 1.8, 'HTC': 1.55, 'RE': 9.0, 'TI': 3.6}
cfg.rateBkg = {'exc': 40, 'inh': 40}

## options to provide external sensory input
#cfg.randomThalInput = True  # provide random bkg inputs spikes (NetStim) to thalamic populations 

cfg.EbkgThalamicGain = 4.0
cfg.IbkgThalamicGain = 4.0

cfg.cochlearThalInput = False #{'numCells': 200, 'freqRange': [9*1e3, 11*1e3], 'toneFreq': 10*1e3, 'loudnessDBs': 50}  # parameters to generate realistic  auditory thalamic inputs using Brian Hears 

# parameters to generate realistic cochlear + IC input ; weight =unitary connection somatic EPSP (mV)
cfg.ICThalInput = {} #'file': 'data/ICoutput/ICoutput_CF_9600_10400_wav_01_ba_peter.mat', 
                    #'startTime': 500, 'weightE': 0.5, 'weightI': 0.5, 'probE': 0.12, 'probI': 0.26, 'seed': 1}  

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = 0

#------------------------------------------------------------------------------
# NetStim inputs 
#------------------------------------------------------------------------------

cfg.addNetStim = 0

## LAYER 1
cfg.NetStim1 = {'pop': 'NGF1', 'ynorm': [0,2.0], 'sec': 'soma', 'loc': 0.5, 'synMech': ['AMPA'], 'synMechWeightFactor': [1.0], 'start': 0, 'interval': 1000.0/60.0, 'noise': 0.0, 'number': 0.0, 'weight': 10.0, 'delay': 0}


cfg.tune = {}





