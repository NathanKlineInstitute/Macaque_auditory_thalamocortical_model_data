"""
batch.py 

Batch simulation for A1 model using NetPyNE

Contributors: salvadordura@gmail.com
"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np


# ----------------------------------------------------------------------------------------------
# Weight Normalization 
# ----------------------------------------------------------------------------------------------
def bkgWeights(pops=[], weights=list(range(50))):

    params = specs.ODict()
    params['singlePop'] = pops
    params['weightBkg'] = weights

    # set initial config
    initCfg = {}
    # sim and recoring params
    initCfg['duration'] = 10.0 * 1e3
    initCfg['singleCellPops'] = True
    initCfg['singlePopForNetstim'] = True
    initCfg['removeWeightNorm'] = False
    initCfg[('analysis','plotTraces','include')] = [0]
    initCfg[('analysis','plotTraces','timeRange')] = [0, 3000]
    initCfg[('analysis', 'plotRaster')] = False

    initCfg[('rateBkg', 'exc')] = 40
    initCfg[('rateBkg', 'inh')] = 40

    ## turn off components not required
    initCfg['addBkgConn'] = True
    initCfg['addConn'] = False
    initCfg['addIntraThalamicConn'] = False
    initCfg['addIntraThalamicConn'] = False
    initCfg['addCorticoThalamicConn'] = False
    initCfg['addCoreThalamoCorticalConn'] = False
    initCfg['addMatrixThalamoCorticalConn'] = False
    initCfg['stimSubConn'] = False
    initCfg['addIClamp'] = False
    initCfg['addNetStim'] = False
 
    b = Batch(params=params, netParamsFile='netParams_bkg.py', cfgFile='cfg_cell.py', initCfg=initCfg)
    b.method = 'grid'

    return b


# ----------------------------------------------------------------------------------------------
# Weight Normalization 
# ----------------------------------------------------------------------------------------------
def weightNorm(pops=[], rule = None, segs = None, allSegs = True, weights=list(np.arange(0.01, 0.2, 0.01)/100.0)):

    # Add params
    from cfg_cell import cfg
    from netParams_cell import netParams

    excludeSegs = ['axon']
    if not segs:
        secs = []
        locs = []
        for secName,sec in netParams.cellParams[rule]['secs'].items():
            if secName not in excludeSegs:
                if allSegs:
                    nseg = sec['geom']['nseg']
                    for iseg in range(nseg):
                        secs.append(secName) 
                        locs.append((iseg+1)*(1.0/(nseg+1)))
                else:
                    secs.append(secName) 
                    locs.append(0.5)

    params = specs.ODict()
    params[('NetStim1', 'pop')] = pops
    params[('NetStim1', 'sec')] = secs
    params[('NetStim1', 'loc')] = locs
    params[('NetStim1', 'weight')] = weights

    groupedParams = [('NetStim1', 'sec'), ('NetStim1', 'loc')] 

    # set initial config
    initCfg = {}
    # sim and recoring params
    initCfg['duration'] = 1.0 * 1e3
    initCfg['singleCellPops'] = True
    initCfg['removeWeightNorm'] = True
    initCfg[('analysis','plotTraces','include')] = []
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]
    
    ## turn off components not required
    #initCfg[('analysis', 'plotRaster')] = False
    initCfg['addConn'] = False
    initCfg['addIntraThalamicConn'] = False
    initCfg['addIntraThalamicConn'] = False
    initCfg['addCorticoThalamicConn'] = False
    initCfg['addCoreThalamoCorticalConn'] = False
    initCfg['addMatrixThalamoCorticalConn'] = False
    initCfg['addBkgConn'] = False
    initCfg['stimSubConn'] = False
    initCfg['addIClamp'] = 0
 
    ## set netstim params
    initCfg['addNetStim'] = True
    initCfg[('NetStim1', 'synMech')] = ['AMPA','NMDA']
    initCfg[('NetStim1','synMechWeightFactor')] = [0.5,0.5]
    initCfg[('NetStim1', 'start')] = 700
    initCfg[('NetStim1', 'interval')] = 1000
    initCfg[('NetStim1','ynorm')] = [0.0, 2.0]
    initCfg[('NetStim1', 'noise')] = 0
    initCfg[('NetStim1', 'number')] = 1
    initCfg[('NetStim1', 'delay')] = 1
    
    
    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)
    b.method = 'grid'

    return b

# ----------------------------------------------------------------------------------------------
# f-I curve
# ----------------------------------------------------------------------------------------------
def fIcurve(pops = [], amps = list(np.arange(0.0, 6.5, 0.5)/10.0) ):
    params = specs.ODict()

    params['singlePop'] = pops
    params[('IClamp1', 'amp')] = amps
    #params['ihGbar'] = [0.0, 1.0, 2.0]
    # params['axonNa'] = [5, 6, 7, 8] 
    # params['gpas'] = [0.6, 0.65, 0.70, 0.75] 
    # params['epas'] = [1.0, 1.05] 
    # params['ihLkcBasal'] = [0.0, 0.01, 0.1, 0.5, 1.0] 

    # initial config
    initCfg = {}
    initCfg['duration'] = 2.0*1e3
    initCfg['addIClamp'] = True
    initCfg['addNetStim'] = False
    initCfg['weightNorm'] = True
    initCfg[('IClamp1','sec')] = 'soma'
    initCfg[('IClamp1','loc')] = 0.5
    initCfg[('IClamp1','start')] = 750
    initCfg[('IClamp1','dur')] = 1000
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [0, 2000]
    initCfg['printPopAvgRates'] = [750,1750]

    initCfg[('hParams', 'celsius')] = 37

    ## turn off components not required
    initCfg['addBkgConn'] = False
    initCfg['addConn'] = False
    initCfg['addIntraThalamicConn'] = False
    initCfg['addIntraThalamicConn'] = False
    initCfg['addCorticoThalamicConn'] = False
    initCfg['addCoreThalamoCorticalConn'] = False
    initCfg['addMatrixThalamoCorticalConn'] = False
    initCfg['stimSubConn'] = False
    initCfg['addNetStim'] = False

    groupedParams = [] 

    b = Batch(params=params, netParamsFile='netParams_cell.py', cfgFile='cfg_cell.py', initCfg=initCfg, groupedParams=groupedParams)
    b.method = 'grid'

    return b



# ----------------------------------------------------------------------------------------------
# Custom
# ----------------------------------------------------------------------------------------------
def custom_spont(filename):
    params = specs.ODict()

    if not filename:
        filename = 'data/v34_batch25/trial_2142/trial_2142_cfg.json'

    # from prev 
    import json
    with open(filename, 'rb') as f:
        cfgLoad = json.load(f)['simConfig']
    cfgLoad2 = cfgLoad

    
    params[('seeds', 'conn')] = list(range(1)) #[4321+(17*i) for i in range(5)]
    params[('seeds', 'stim')] = list(range(1)) #[1234+(17*i) for i in range(5)]

    groupedParams = [] 

    # --------------------------------------------------------
    # initial config
    initCfg = {} # set default options from prev sim
    
    initCfg['duration'] = 11500
    initCfg['printPopAvgRates'] = [1500, initCfg['duration']] 
    initCfg['scaleDensity'] = 1.0
    initCfg['recordStep'] = 0.05

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotTraces', 'timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotSpikeStats', 'timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotLFP', 'timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotCSD', 'timeRange')] = [1500, 1700]

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False
    
    # from prev 
    updateParams = ['EEGain', 'EIGain', 'IEGain', 'IIGain',
                    ('EICellTypeGain', 'PV'), ('EICellTypeGain', 'SOM'), ('EICellTypeGain', 'VIP'), ('EICellTypeGain', 'NGF'),
                    ('IECellTypeGain', 'PV'), ('IECellTypeGain', 'SOM'), ('IECellTypeGain', 'VIP'), ('IECellTypeGain', 'NGF'),
                    ('EILayerGain', '1'), ('IILayerGain', '1'),
                    ('EELayerGain', '2'), ('EILayerGain', '2'),  ('IELayerGain', '2'), ('IILayerGain', '2'), 
                    ('EELayerGain', '3'), ('EILayerGain', '3'), ('IELayerGain', '3'), ('IILayerGain', '3'), 
                    ('EELayerGain', '4'), ('EILayerGain', '4'), ('IELayerGain', '4'), ('IILayerGain', '4'), 
                    ('EELayerGain', '5A'), ('EILayerGain', '5A'), ('IELayerGain', '5A'), ('IILayerGain', '5A'), 
                    ('EELayerGain', '5B'), ('EILayerGain', '5B'), ('IELayerGain', '5B'), ('IILayerGain', '5B'), 
                    ('EELayerGain', '6'), ('EILayerGain', '6'), ('IELayerGain', '6'), ('IILayerGain', '6')] 

    for p in updateParams:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad[p]})

    # good thal params for 100% cell density 
    updateParams2 = ['thalamoCorticalGain', 'intraThalamicGain', 'EbkgThalamicGain', 'IbkgThalamicGain', 'wmat']

    for p in updateParams2:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad2[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad2[p]})


    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)
    b.method = 'grid'

    return b


# ----------------------------------------------------------------------------------------------
# Custom
# ----------------------------------------------------------------------------------------------
def custom_speech(filename):
    params = specs.ODict()

    if not filename:
        filename = 'data/v34_batch25/trial_2142/trial_2142_cfg.json'

    # from prev 
    import json
    with open(filename, 'rb') as f:
        cfgLoad = json.load(f)['simConfig']
    cfgLoad2 = cfgLoad

    params[('seeds', 'conn')] = [4321+(17*i) for i in range(5)]
    params[('seeds', 'stim')] = [1234+(17*i) for i in range(5)]
    
    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {} # set default options from prev sim
    
    initCfg['duration'] = 4500
    initCfg['printPopAvgRates'] = [1500, 4500] 
    initCfg['scaleDensity'] = 1.0
    initCfg['recordStep'] = 0.05

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotTraces', 'timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotSpikeStats', 'timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotLFP', 'timeRange')] = initCfg['printPopAvgRates']
    #initCfg[('analysis', 'plotCSD', 'timeRange')] = [1500, 1700]

    initCfg['ICThalInput'] = {'file': 'data/ICoutput/ICoutput_CF_9600_10400_wav_01_ba_peter.mat', 
                            'startTime': 2500, 
                            'weightE': 0.25,#1.0, 
                            'weightI': 0.25,#1.0, 
                            'probE': 0.12, 
                            'probI': 0.12, #0.25 
                            'seed': 1}  



    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False
    
    # from prev - best of 50% cell density
    updateParams = ['EEGain', 'EIGain', 'IEGain', 'IIGain',
                    ('EICellTypeGain', 'PV'), ('EICellTypeGain', 'SOM'), ('EICellTypeGain', 'VIP'), ('EICellTypeGain', 'NGF'),
                    ('IECellTypeGain', 'PV'), ('IECellTypeGain', 'SOM'), ('IECellTypeGain', 'VIP'), ('IECellTypeGain', 'NGF'),
                    ('EILayerGain', '1'), ('IILayerGain', '1'),
                    ('EELayerGain', '2'), ('EILayerGain', '2'),  ('IELayerGain', '2'), ('IILayerGain', '2'), 
                    ('EELayerGain', '3'), ('EILayerGain', '3'), ('IELayerGain', '3'), ('IILayerGain', '3'), 
                    ('EELayerGain', '4'), ('EILayerGain', '4'), ('IELayerGain', '4'), ('IILayerGain', '4'), 
                    ('EELayerGain', '5A'), ('EILayerGain', '5A'), ('IELayerGain', '5A'), ('IILayerGain', '5A'), 
                    ('EELayerGain', '5B'), ('EILayerGain', '5B'), ('IELayerGain', '5B'), ('IILayerGain', '5B'), 
                    ('EELayerGain', '6'), ('EILayerGain', '6'), ('IELayerGain', '6'), ('IILayerGain', '6')] 

    for p in updateParams:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad[p]})

    # good thal params for 100% cell density 
    updateParams2 = ['thalamoCorticalGain', 'intraThalamicGain', 'EbkgThalamicGain', 'IbkgThalamicGain', 'wmat']

    for p in updateParams2:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad2[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad2[p]})


    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)
    b.method = 'grid'

    return b



# ----------------------------------------------------------------------------------------------
# Evol
# ----------------------------------------------------------------------------------------------
def evolRates():
    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    # bkg inputs
    params['EEGain'] = [0.5, 2.0]
    params['EIGain'] = [0.5, 2.0]

    params[('IELayerGain', '1-3')] = [0.5, 2.0]
    params[('IELayerGain', '4')] = [0.5, 2.0]
    params[('IELayerGain', '5')] = [0.5, 2.0]
    params[('IELayerGain', '6')] = [0.5, 2.0]

    params[('IILayerGain', '1-3')] = [0.5, 2.0]
    params[('IILayerGain', '4')] = [0.5, 2.0]
    params[('IILayerGain', '5')] = [0.5, 2.0]
    params[('IILayerGain', '6')] = [0.5, 2.0]
    
    params['thalamoCorticalGain'] = [0.5, 2.0]  
    params['intraThalamicGain'] = [0.5, 2.0] 
    params['corticoThalamicGain'] = [0.5, 2.0]

    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg = {}
    initCfg['duration'] = 1500
    initCfg['printPopAvgRates'] = [[500, 750], [750, 1000], [1000, 1250], [1250, 1500]]
    initCfg['dt'] = 0.05

    initCfg['scaleDensity'] = 0.5

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = [500,1500]
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [500,1500]

    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}
    
    ## Exc pops
    Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'TCM', 'HTC']  # all layers + thal + IC

    Etune = {'target': 5, 'width': 20, 'min': 0.05}
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    Ipops = ['NGF1',                            # L1
            'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
            'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
            'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
            'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
            'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
            'PV6', 'SOM6', 'VIP6', 'NGF6',       # L6
            'IRE', 'IREM', 'TI']  # Thal 

    Itune = {'target': 10, 'width': 30, 'min': 0.05}
    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 1000
    fitnessFuncArgs['tranges'] = initCfg['printPopAvgRates']


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        tranges = kwargs['tranges']

        popFitnessAll = []

        for trange in tranges:
            popFitnessAll.append([min(np.exp(abs(v['target'] - simData['popRates'][k]['%d_%d'%(trange[0], trange[1])])/v['width']), maxFitness) 
                if simData['popRates'][k]['%d_%d'%(trange[0], trange[1])] > v['min'] else maxFitness for k, v in pops.items()])
        
        popFitness = np.mean(np.array(popFitnessAll), axis=0)
        
        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f' % (p, np.mean(list(simData['popRates'][p].values())), popFitness[i]) for i,p in enumerate(pops)])
        print('  ' + popInfo)

        return fitness

    
    #from IPython import embed; embed()

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.method = 'evol'

    # Set evol alg configuration
    b.evolCfg = {
        'evolAlgorithm': 'custom',
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'pop_size': 100,
        'num_elites': 2,
        'mutation_rate': 0.5,
        'crossover': 0.5,
        'maximize': False, # maximize fitness function?
        'max_generations': 200,
        'time_sleep': 150, # 2.5min wait this time before checking again if sim is completed (for each generation)
        'maxiter_wait': 5, # max number of times to check if sim is completed (for each generation)
        'defaultFitness': 1000, # set fitness value in case simulation time is over
        'scancelUser': 'ext_salvadordura_gmail_com'
    }

    return b


# ----------------------------------------------------------------------------------------------
# Adaptive Stochastic Descent (ASD)
# ----------------------------------------------------------------------------------------------
def optunaRates():

    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    # bkg inputs
    params['EEGain'] = [0.1, 1.0]

    params[('EILayerGain', '1-3')] = [0.1, 3.0]
    params[('EILayerGain', '4')] = [0.1, 3.0]
    params[('EILayerGain', '5')] = [0.1, 3.0]
    params[('EILayerGain', '6')] = [0.1, 3.0]

    params[('IELayerGain', '1-3')] = [0.1, 3.0]
    params[('IELayerGain', '4')] = [0.1, 3.0]
    params[('IELayerGain', '5')] = [0.1, 3.0]
    params[('IELayerGain', '6')] = [0.1, 3.0]

    params[('IILayerGain', '1-3')] = [0.1, 3.0]
    params[('IILayerGain', '4')] = [0.1, 3.0]
    params[('IILayerGain', '5')] = [0.1, 3.0]
    params[('IILayerGain', '6')] = [0.1, 3.0]

    params[('EICellTypeGain', 'PV')] = [max(cfgLoad['EICellTypeGain']['PV']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['PV']+rangeV2, maxV)]
    params[('EICellTypeGain', 'SOM')] = [max(cfgLoad['EICellTypeGain']['SOM']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['SOM']+rangeV2, maxV)]
    params[('EICellTypeGain', 'VIP')] = [max(cfgLoad['EICellTypeGain']['VIP']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['VIP']+rangeV2, maxV)]
    params[('EICellTypeGain', 'NGF')] = [max(cfgLoad['EICellTypeGain']['NGF']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['NGF']+rangeV2, maxV)]

    # params['thalamoCorticalGain'] = [0.25, 2.0]
    # params['intraThalamicGain'] = [0.25, 2.0]
    # params['corticoThalamicGain'] = [0.25, 2.0]


    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg = {}
    initCfg['duration'] = 2000
    initCfg['printPopAvgRates'] = [[1000, 1250], [1250, 1500], [1500, 1750], [1750, 2000]]
    initCfg['dt'] = 0.05

    initCfg['scaleDensity'] = 0.5

    # plotting and saving params
    initCfg[('analysis','plotRaster','timeRange')] = [1000,2000]
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [1000,2000]
    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'
    initCfg['recordLFP'] = None
    initCfg[('analysis', 'plotLFP')] = False

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}
    
    ## Exc pops
    Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'TCM', 'HTC']  # all layers + thal + IC

    Etune = {'target': 5, 'width': 20, 'min': 0.05}
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    Ipops = ['NGF1',                            # L1
            'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
            'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
            'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
            'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
            'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
            'PV6', 'SOM6', 'VIP6', 'NGF6',       # L6
            'IRE', 'IREM', 'TI']  # Thal 

    Itune = {'target': 10, 'width': 30, 'min': 0.05}
    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 1000
    fitnessFuncArgs['tranges'] = initCfg['printPopAvgRates']


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        tranges = kwargs['tranges']

        popFitnessAll = []

        for trange in tranges:
            popFitnessAll.append([min(np.exp(abs(v['target'] - simData['popRates'][k]['%d_%d'%(trange[0], trange[1])])/v['width']), maxFitness) 
                if simData['popRates'][k]['%d_%d'%(trange[0], trange[1])] > v['min'] else maxFitness for k, v in pops.items()])
        
        popFitness = np.mean(np.array(popFitnessAll), axis=0)
        
        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f' % (p, np.mean(list(simData['popRates'][p].values())), popFitness[i]) for i,p in enumerate(pops)])
        print('  ' + popInfo)

        return fitness
        
    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.method = 'optuna'

    b.optimCfg = {
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'maxFitness': fitnessFuncArgs['maxFitness'],
        'maxiters':     1e6,    #    Maximum number of iterations (1 iteration = 1 function evaluation)
        'maxtime':      None,    #    Maximum time allowed, in seconds
        'maxiter_wait': 16,
        'time_sleep': 60,
        'popsize': 1  # unused - run with mpi 
    }

    return b



# ----------------------------------------------------------------------------------------------
# Optuna optimization
# ----------------------------------------------------------------------------------------------
def optunaRatesLayers():

    '''
    # from prev
    import json
    with open('data/v32_batch14/trial_981/trial_981_cfg.json', 'rb') as f:
        cfgLoad = json.load(f)['simConfig']
    '''

    # best from grid search
    import json
    # with open('data/v32_batch4/trial_15057/trial_15057_cfg.json', 'rb') as f:
    #     cfgLoad = json.load(f)['simConfig']

    # good thal params for 100% cell density 
    with open('data/v34_batch5/v34_batch5_0_2_0_0_2_0_2_2_cfg.json', 'rb') as f:
        cfgLoad2 = json.load(f)['simConfig']

    cfgLoad = cfgLoad2

    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    rangeV = 0.25
    rangeV2 = 1.0
    minV = 0.1
    maxV = 4.0


    # params based on v32_batch8 grid search best solutions; plus added L2 and L3 IE specific gains since those were problematic layers
    params['EEGain'] = [0.1, 0.5] 
    params['EIGain'] = [1.3, 1.7] 
    params['IEGain'] = [0.8, 1.7] 
    params['IIGain'] = [0.3, 0.7] 
    params[('EICellTypeGain', 'PV')] = [0.8, 1.7] 
    params[('EICellTypeGain', 'SOM')] = [0.3, 0.7] 
    params[('EICellTypeGain', 'VIP')] = [1.3, 1.7] 
    params[('EICellTypeGain', 'NGF')] = [1.3, 1.7] 

    params[('IELayerGain', '2')] = [max(cfgLoad2['IELayerGain']['2']-rangeV, minV), min(cfgLoad2['IELayerGain']['2']+rangeV, maxV)]
    params[('IELayerGain', '3')] = [max(cfgLoad2['IELayerGain']['3']-rangeV, minV), min(cfgLoad2['IELayerGain']['3']+rangeV, maxV)]
    params[('IELayerGain', '4')] = [max(cfgLoad2['IELayerGain']['4']-rangeV, minV), min(cfgLoad2['IELayerGain']['4']+rangeV, maxV)]
    params[('IELayerGain', '5A')] = [max(cfgLoad2['IELayerGain']['5A']-rangeV, minV), min(cfgLoad2['IELayerGain']['5A']+rangeV, maxV)]
    params[('IELayerGain', '5B')] = [max(cfgLoad2['IELayerGain']['5B']-rangeV, minV), min(cfgLoad2['IELayerGain']['5B']+rangeV, maxV)]

    
    '''

    params[('EICellTypeGain', 'PV')] = [max(cfgLoad['EICellTypeGain']['PV']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['PV']+rangeV2, maxV)]
    params[('EICellTypeGain', 'SOM')] = [max(cfgLoad['EICellTypeGain']['SOM']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['SOM']+rangeV2, maxV)]
    params[('EICellTypeGain', 'VIP')] = [max(cfgLoad['EICellTypeGain']['VIP']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['VIP']+rangeV2, maxV)]
    params[('EICellTypeGain', 'NGF')] = [max(cfgLoad['EICellTypeGain']['NGF']-rangeV2, minV), min(cfgLoad['EICellTypeGain']['NGF']+rangeV2, maxV)]

    params[('IECellTypeGain', 'PV')] = [max(cfgLoad['IECellTypeGain']['PV']-rangeV, minV), min(cfgLoad['IECellTypeGain']['PV']+rangeV, maxV)]
    params[('IECellTypeGain', 'SOM')] = [max(cfgLoad['IECellTypeGain']['SOM']-rangeV, minV), min(cfgLoad['IECellTypeGain']['SOM']+rangeV, maxV)]
    params[('IECellTypeGain', 'VIP')] = [max(cfgLoad['IECellTypeGain']['VIP']-rangeV, minV), min(cfgLoad['IECellTypeGain']['VIP']+rangeV, maxV)]
    params[('IECellTypeGain', 'NGF')] = [max(cfgLoad['IECellTypeGain']['NGF']-rangeV, minV), min(cfgLoad['IECellTypeGain']['NGF']+rangeV, maxV)]

    params[('EILayerGain', '1')] = [minV, maxV] #(cfgLoad['EILayerGain']['1']-rangeV2, minV), min(cfgLoad['EILayerGain']['1']+rangeV, maxV)]
    params[('IILayerGain', '1')] = [minV, maxV] #(cfgLoad['IILayerGain']['1']-rangeV2, minV), min(cfgLoad['IILayerGain']['1']+rangeV, maxV)]

    params[('EELayerGain', '2')] = [max(cfgLoad['EELayerGain']['2']-rangeV, minV), min(cfgLoad['EELayerGain']['2']+rangeV, maxV)]
    params[('EILayerGain', '2')] = [max(cfgLoad['EILayerGain']['2']-rangeV, minV), min(cfgLoad['EILayerGain']['2']+rangeV, maxV)]
    params[('IELayerGain', '2')] = [max(cfgLoad['IELayerGain']['2']-rangeV, minV), min(cfgLoad['IELayerGain']['2']+rangeV, maxV)]
    params[('IILayerGain', '2')] = [max(cfgLoad['IILayerGain']['2']-rangeV, minV), min(cfgLoad['IILayerGain']['2']+rangeV, maxV)]

    params[('EELayerGain', '3')] = [max(cfgLoad['EELayerGain']['3']-rangeV, minV), min(cfgLoad['EELayerGain']['3']+rangeV, maxV)]
    params[('EILayerGain', '3')] = [max(cfgLoad['EILayerGain']['3']-rangeV, minV), min(cfgLoad['EILayerGain']['3']+rangeV, maxV)]
    params[('IELayerGain', '3')] = [max(cfgLoad['IELayerGain']['3']-rangeV, minV), min(cfgLoad['IELayerGain']['3']+rangeV, maxV)]
    params[('IILayerGain', '3')] = [max(cfgLoad['IILayerGain']['3']-rangeV, minV), min(cfgLoad['IILayerGain']['3']+rangeV, maxV)]

    params[('EELayerGain', '4')] = [max(cfgLoad['EELayerGain']['4']-rangeV, minV), min(cfgLoad['EELayerGain']['4']+rangeV, maxV)]
    params[('EILayerGain', '4')] = [max(cfgLoad['EILayerGain']['4']-rangeV, minV), min(cfgLoad['EILayerGain']['4']+rangeV, maxV)]
    params[('IELayerGain', '4')] = [max(cfgLoad['IELayerGain']['4']-rangeV, minV), min(cfgLoad['IELayerGain']['4']+rangeV, maxV)]
    params[('IILayerGain', '4')] = [max(cfgLoad['IILayerGain']['4']-rangeV, minV), min(cfgLoad['IILayerGain']['4']+rangeV, maxV)]

    params[('EELayerGain', '5A')] = [minV, maxV] # [max(cfgLoad['EELayerGain']['5A']-rangeV2, minV), min(cfgLoad['EELayerGain']['5A']+rangeV, maxV)]
    params[('EILayerGain', '5A')] = [minV, maxV] # [max(cfgLoad['EILayerGain']['5A']-rangeV2, minV), min(cfgLoad['EILayerGain']['5A']+rangeV, maxV)]
    params[('IELayerGain', '5A')] = [minV, maxV] # [max(cfgLoad['IELayerGain']['5A']-rangeV2, minV), min(cfgLoad['IELayerGain']['5A']+rangeV, maxV)]
    params[('IILayerGain', '5A')] = [minV, maxV] # [max(cfgLoad['IILayerGain']['5A']-rangeV2, minV), min(cfgLoad['IILayerGain']['5A']+rangeV, maxV)]

    params[('EELayerGain', '5B')] = [minV, maxV] # [max(cfgLoad['EELayerGain']['5B']-rangeV2, minV), min(cfgLoad['EELayerGain']['5B']+rangeV, maxV)]
    params[('EILayerGain', '5B')] = [minV, maxV] # [max(cfgLoad['EILayerGain']['5B']-rangeV2, minV), min(cfgLoad['EILayerGain']['5B']+rangeV, maxV)]
    params[('IELayerGain', '5B')] = [minV, maxV] # [max(cfgLoad['IELayerGain']['5B']-rangeV2, minV), min(cfgLoad['IELayerGain']['5B']+rangeV, maxV)]
    params[('IILayerGain', '5B')] = [minV, maxV] # [max(cfgLoad['IILayerGain']['5B']-rangeV2, minV), min(cfgLoad['IILayerGain']['5B']+rangeV, maxV)]

    params[('EELayerGain', '6')] = [minV, maxV] # [max(cfgLoad['EELayerGain']['6']-rangeV2, minV), min(cfgLoad['EELayerGain']['6']+rangeV, maxV)]
    params[('EILayerGain', '6')] = [minV, maxV] # [max(cfgLoad['EILayerGain']['6']-rangeV2, minV), min(cfgLoad['EILayerGain']['6']+rangeV, maxV)]
    params[('IELayerGain', '6')] = [minV, maxV] # [max(cfgLoad['IELayerGain']['6']-rangeV2, minV), min(cfgLoad['IELayerGain']['6']+rangeV, maxV)]
    params[('IILayerGain', '6')] = [minV, maxV] # [max(cfgLoad['IILayerGain']['6']-rangeV2, minV), min(cfgLoad['IILayerGain']['6']+rangeV, maxV)]
    '''

    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg = {}
    initCfg['duration'] = 2500
    initCfg['printPopAvgRates'] = [[1500, 1750], [1750, 2000], [2000, 2250], [2250, 2500]]
    initCfg['dt'] = 0.05

    initCfg['scaleDensity'] = 1.0

    # plotting and saving params
    initCfg[('analysis','plotRaster','markerSize')] = 10

    initCfg[('analysis','plotRaster','timeRange')] = [1500, 2500]
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [1500, 2500]
    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'
    initCfg['recordLFP'] = None
    initCfg[('analysis', 'plotLFP')] = False

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False
    
    '''
    initCfg['EEGain'] = 1.0 	
    initCfg['EIGain'] = 1.0 	
    initCfg['IEGain'] = 1.0 	
    initCfg['IIGain'] = 1.0 	

    initCfg.update({'thalamoCorticalGain': cfgLoad['thalamoCorticalGain'],
                    'intraThalamicGain': cfgLoad['intraThalamicGain'],
                    'EbkgThalamicGain': cfgLoad['EbkgThalamicGain'],
                    'IbkgThalamicGain': cfgLoad['IbkgThalamicGain']})

    print(initCfg)
    '''

    # from prev - best of 50% cell density
    updateParams = ['EEGain', 'EIGain', 'IEGain', 'IIGain',
                    ('EICellTypeGain', 'PV'), ('EICellTypeGain', 'SOM'), ('EICellTypeGain', 'VIP'), ('EICellTypeGain', 'NGF'),
                    ('IECellTypeGain', 'PV'), ('IECellTypeGain', 'SOM'), ('IECellTypeGain', 'VIP'), ('IECellTypeGain', 'NGF'),
                    ('EILayerGain', '1'), ('IILayerGain', '1'),
                    ('EELayerGain', '2'), ('EILayerGain', '2'),  ('IELayerGain', '2'), ('IILayerGain', '2'), 
                    ('EELayerGain', '3'), ('EILayerGain', '3'), ('IELayerGain', '3'), ('IILayerGain', '3'), 
                    ('EELayerGain', '4'), ('EILayerGain', '4'), ('IELayerGain', '4'), ('IILayerGain', '4'), 
                    ('EELayerGain', '5A'), ('EILayerGain', '5A'), ('IELayerGain', '5A'), ('IILayerGain', '5A'), 
                    ('EELayerGain', '5B'), ('EILayerGain', '5B'), ('IELayerGain', '5B'), ('IILayerGain', '5B'), 
                    ('EELayerGain', '6'), ('EILayerGain', '6'), ('IELayerGain', '6'), ('IILayerGain', '6')] 

    for p in updateParams:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad[p]})

    # good thal params for 100% cell density 
    updateParams2 = ['thalamoCorticalGain', 'intraThalamicGain', 'EbkgThalamicGain', 'IbkgThalamicGain']

    for p in updateParams2:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad2[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad2[p]})


    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}
    
    ## Exc pops
    Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'TCM', 'HTC']  # all layers + thal + IC
    #Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'TC', 'TCM', 'HTC']  # all layers + thal + IC

    #Etune = {'target': 5, 'width': 20, 'min': 0.05}
    Etune = {'target': 5, 'width': 5, 'min': 0.5}
    
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    Ipops = ['NGF1',                            # L1
            'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
            'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
            'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
            'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
            'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
            'PV6', 'SOM6', 'VIP6', 'NGF6',       # L6
            'IRE', 'IREM', 'TI']  # Thal 
    # Ipops = [#'NGF1',  
    #         'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
    #         'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
    #         'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
    #         #'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A 
    #         #'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
    #         #'PV6', 'SOM6', 'VIP6', 'NGF6',  # L6
    #         'IRE', 'IREM', 'TI']  # Thal 

    #Itune = {'target': 10, 'width': 30, 'min': 0.05}
    Itune = {'target': 10, 'width': 15, 'min': 0.5}

    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 1000
    fitnessFuncArgs['tranges'] = initCfg['printPopAvgRates']


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        tranges = kwargs['tranges']

        popFitnessAll = []

        for trange in tranges:
            popFitnessAll.append([min(np.exp(abs(v['target'] - simData['popRates'][k]['%d_%d'%(trange[0], trange[1])])/v['width']), maxFitness) 
                if simData['popRates'][k]['%d_%d'%(trange[0], trange[1])] >= v['min'] else maxFitness for k, v in pops.items()])
        
        popFitness = np.mean(np.array(popFitnessAll), axis=0)
        
        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f' % (p, np.mean(list(simData['popRates'][p].values())), popFitness[i]) for i,p in enumerate(pops)])
        print('  ' + popInfo)

        return fitness
        
    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.method = 'optuna'

    b.optimCfg = {
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'maxFitness': fitnessFuncArgs['maxFitness'],
        'maxiters':     1e6,    #    Maximum number of iterations (1 iteration = 1 function evaluation)
        'maxtime':      None,    #    Maximum time allowed, in seconds
        'maxiter_wait': 45,
        'time_sleep': 120,
        'popsize': 1  # unused - run with mpi 
    }

    return b



# ----------------------------------------------------------------------------------------------
# Optuna optimization
# ----------------------------------------------------------------------------------------------
def optunaRatesLayersThalL12345A5B6():

    # from prev
    import json
    with open('data/v34_batch15/trial_5955/trial_5955_cfg.json', 'rb') as f:
        cfgLoad = json.load(f)['simConfig']


    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    rangeV1 = 0.2
    rangeV2 = 0.25
    rangeV3 = 2.0

    scaleLow = 0.9
    scaleHigh = 1.1

    scaleLow2 = 0.5
    scaleHigh2 = 1.5

    minV = 0.1
    maxV = 5.0

    ''' potentially:
    layer-specific EE,EI,IE,II gains for L2, with [prev-0.25:prev+0.25]
    layer-specific EE,EI,IE,II gains for L3, with [prev-0.25:prev+025]
    layer-specific EE,EI,IE,II gains for L4, with [prev-0.1:prev+0.1]
    layer-specific EE,EI,IE,II gains for L5A, with param range [prev-0.1:prev+0.1]
    layer-specific EE,EI,IE,II gains for L5B, with broad param range [0.1:4.0]
    E-PV, E-SOM, E-NGF, E-VIP (common to all layers) with range [prev-0.25:prev+0.25]
    PV-E, SOM-E, NGF-E, VIP-E (common to all layers) with range [0.1:4.0]

    '''

    # E->I cell-type-specific
    params[('EICellTypeGain', 'PV')] = [max(cfgLoad['EICellTypeGain']['PV']*scaleLow2, minV), min(cfgLoad['EICellTypeGain']['PV']*scaleHigh2, maxV)]
    params[('EICellTypeGain', 'SOM')] = [max(cfgLoad['EICellTypeGain']['SOM']*scaleLow2, minV), min(cfgLoad['EICellTypeGain']['SOM']*scaleHigh2, maxV)]
    params[('EICellTypeGain', 'VIP')] = [max(cfgLoad['EICellTypeGain']['VIP']*scaleLow2, minV), min(cfgLoad['EICellTypeGain']['VIP']*scaleHigh2, maxV)]
    params[('EICellTypeGain', 'NGF')] = [max(cfgLoad['EICellTypeGain']['NGF']*scaleLow2, minV), min(cfgLoad['EICellTypeGain']['NGF']*scaleHigh2, maxV)]

    # I->E cell-type-specific
    params[('IECellTypeGain', 'PV')] = [max(cfgLoad['IECellTypeGain']['PV']*scaleLow, minV), min(cfgLoad['IECellTypeGain']['PV']*scaleHigh, maxV)]
    params[('IECellTypeGain', 'SOM')] = [max(cfgLoad['IECellTypeGain']['SOM']*scaleLow, minV), min(cfgLoad['IECellTypeGain']['SOM']*scaleHigh, maxV)]
    params[('IECellTypeGain', 'VIP')] = [max(cfgLoad['IECellTypeGain']['VIP']*scaleLow, minV), min(cfgLoad['IECellTypeGain']['VIP']*scaleHigh, maxV)]
    params[('IECellTypeGain', 'NGF')] = [max(cfgLoad['IECellTypeGain']['NGF']*scaleLow, minV), min(cfgLoad['IECellTypeGain']['NGF']*scaleHigh, maxV)]

    # L1
    params[('EILayerGain', '1')] = [max(cfgLoad['EILayerGain']['1']*scaleLow, minV), min(cfgLoad['EILayerGain']['1']*scaleHigh, maxV)]
    params[('IILayerGain', '1')] = [max(cfgLoad['IILayerGain']['1']*scaleLow, minV), min(cfgLoad['IILayerGain']['1']*scaleHigh, maxV)]

    # L2
    params[('EELayerGain', '2')] = [max(cfgLoad['EELayerGain']['2']*scaleLow2, minV), min(cfgLoad['EELayerGain']['2']*scaleHigh2, maxV)]
    params[('EILayerGain', '2')] = [max(cfgLoad['EILayerGain']['2']*scaleLow2, minV), min(cfgLoad['EILayerGain']['2']*scaleHigh2, maxV)]
    params[('IELayerGain', '2')] = [max(cfgLoad['IELayerGain']['2']*scaleLow2, minV), min(cfgLoad['IELayerGain']['2']*scaleHigh2, maxV)]
    params[('IILayerGain', '2')] = [max(cfgLoad['IILayerGain']['2']*scaleLow2, minV), min(cfgLoad['IILayerGain']['2']*scaleHigh2, maxV)]

    # L3
    params[('EELayerGain', '3')] = [max(cfgLoad['EELayerGain']['3']*scaleLow2, minV), min(cfgLoad['EELayerGain']['3']*scaleHigh2, maxV)]
    params[('EILayerGain', '3')] = [max(cfgLoad['EILayerGain']['3']*scaleLow2, minV), min(cfgLoad['EILayerGain']['3']*scaleHigh2, maxV)]
    params[('IELayerGain', '3')] = [max(cfgLoad['IELayerGain']['3']*scaleLow2, minV), min(cfgLoad['IELayerGain']['3']*scaleHigh2, maxV)]
    params[('IILayerGain', '3')] = [max(cfgLoad['IILayerGain']['3']*scaleLow2, minV), min(cfgLoad['IILayerGain']['3']*scaleHigh2, maxV)]

    # L4
    params[('EELayerGain', '4')] = [max(cfgLoad['EELayerGain']['4']*scaleLow, minV), min(cfgLoad['EELayerGain']['4']*scaleHigh, maxV)]
    params[('EILayerGain', '4')] = [max(cfgLoad['EILayerGain']['4']*scaleLow, minV), min(cfgLoad['EILayerGain']['4']*scaleHigh, maxV)]
    params[('IELayerGain', '4')] = [max(cfgLoad['IELayerGain']['4']*scaleLow, minV), min(cfgLoad['IELayerGain']['4']*scaleHigh, maxV)]
    params[('IILayerGain', '4')] = [max(cfgLoad['IILayerGain']['4']*scaleLow, minV), min(cfgLoad['IILayerGain']['4']*scaleHigh, maxV)]
    
    # L5A
    params[('EELayerGain', '5A')] = [max(cfgLoad['EELayerGain']['5A']*scaleLow, minV), min(cfgLoad['EELayerGain']['5A']*scaleHigh, maxV)]
    params[('EILayerGain', '5A')] = [max(cfgLoad['EILayerGain']['5A']*scaleLow, minV), min(cfgLoad['EILayerGain']['5A']*scaleHigh, maxV)]
    params[('IELayerGain', '5A')] = [max(cfgLoad['IELayerGain']['5A']*scaleLow, minV), min(cfgLoad['IELayerGain']['5A']*scaleHigh, maxV)]
    params[('IILayerGain', '5A')] = [max(cfgLoad['IILayerGain']['5A']*scaleLow, minV), min(cfgLoad['IILayerGain']['5A']*scaleHigh, maxV)]

    # L5B
    params[('EELayerGain', '5B')] = [max(cfgLoad['EELayerGain']['5B']*scaleLow, minV), min(cfgLoad['EELayerGain']['5B']*scaleHigh, maxV)]
    params[('EILayerGain', '5B')] = [max(cfgLoad['EILayerGain']['5B']*scaleLow, minV), min(cfgLoad['EILayerGain']['5B']*scaleHigh, maxV)]
    params[('IELayerGain', '5B')] = [max(cfgLoad['IELayerGain']['5B']*scaleLow, minV), min(cfgLoad['IELayerGain']['5B']*scaleHigh, maxV)]
    params[('IILayerGain', '5B')] = [max(cfgLoad['IILayerGain']['5B']*scaleLow, minV), min(cfgLoad['IILayerGain']['5B']*scaleHigh, maxV)]

    # L6
    params[('EELayerGain', '6')] = [max(cfgLoad['EELayerGain']['6']*scaleLow, minV), min(cfgLoad['EELayerGain']['6']*scaleHigh, maxV)]
    params[('EILayerGain', '6')] = [max(cfgLoad['EILayerGain']['6']*scaleLow, minV), min(cfgLoad['EILayerGain']['6']*scaleHigh, maxV)]
    params[('IELayerGain', '6')] = [max(cfgLoad['IELayerGain']['6']*scaleLow, minV), min(cfgLoad['IELayerGain']['6']*scaleHigh, maxV)]
    params[('IILayerGain', '6')] = [max(cfgLoad['IILayerGain']['6']*scaleLow, minV), min(cfgLoad['IILayerGain']['6']*scaleHigh, maxV)]


    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg = {}
    initCfg['duration'] = 2500
    initCfg['printPopAvgRates'] = [[1500, 1750], [1750, 2000], [2000, 2250], [2250, 2500]]
    initCfg['dt'] = 0.05

    initCfg['scaleDensity'] = 1.0

    # plotting and saving params
    initCfg[('analysis','plotRaster','markerSize')] = 10

    initCfg[('analysis','plotRaster','timeRange')] = [1500, 2500]
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [1500, 2500]
    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'
    initCfg['recordLFP'] = None
    initCfg[('analysis', 'plotLFP')] = False

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    # from prev - best of 50% cell density
    updateParams = ['EEGain', 'EIGain', 'IEGain', 'IIGain',
                    ('EICellTypeGain', 'PV'), ('EICellTypeGain', 'SOM'), ('EICellTypeGain', 'VIP'), ('EICellTypeGain', 'NGF'),
                    ('IECellTypeGain', 'PV'), ('IECellTypeGain', 'SOM'), ('IECellTypeGain', 'VIP'), ('IECellTypeGain', 'NGF'),
                    ('EILayerGain', '1'), ('IILayerGain', '1'),
                    ('EELayerGain', '2'), ('EILayerGain', '2'),  ('IELayerGain', '2'), ('IILayerGain', '2'), 
                    ('EELayerGain', '3'), ('EILayerGain', '3'), ('IELayerGain', '3'), ('IILayerGain', '3'), 
                    ('EELayerGain', '4'), ('EILayerGain', '4'), ('IELayerGain', '4'), ('IILayerGain', '4'), 
                    ('EELayerGain', '5A'), ('EILayerGain', '5A'), ('IELayerGain', '5A'), ('IILayerGain', '5A'), 
                    ('EELayerGain', '5B'), ('EILayerGain', '5B'), ('IELayerGain', '5B'), ('IILayerGain', '5B'), 
                    ('EELayerGain', '6'), ('EILayerGain', '6'), ('IELayerGain', '6'), ('IILayerGain', '6')] 

    for p in updateParams:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad[p]})

    # good thal params for 100% cell density 
    updateParams2 = ['thalamoCorticalGain', 'intraThalamicGain', 'EbkgThalamicGain', 'IbkgThalamicGain']

    for p in updateParams2:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad[p]})

    print(initCfg)


    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}
    
    ## Exc pops
    Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'TCM', 'HTC']  # all layers + thal + IC
    # Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'TC', 'TCM', 'HTC']  # all layers + thal + IC

    Etune = {'target': 5, 'width': 20, 'min': 0.05}
    #Etune = {'target': 5, 'width': 5, 'min': 0.5}
    
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    # Ipops = ['NGF1',                            # L1
    #         'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
    #         'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
    #         'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
    #         'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
    #         'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
    #         'PV6', 'SOM6', 'VIP6', 'NGF6',       # L6
    #         'IRE', 'IREM', 'TI']  # Thal 
    Ipops = ['NGF1',  
            'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
            'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
            'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
            'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A 
            'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
            'PV6', 'SOM6', 'VIP6', 'NGF6',  # L6
            'IRE', 'IREM', 'TI']  # Thal 

    Itune = {'target': 10, 'width': 30, 'min': 0.05}
    #Itune = {'target': 10, 'width': 15, 'min': 0.5}

    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 2000
    fitnessFuncArgs['tranges'] = initCfg['printPopAvgRates']


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        tranges = kwargs['tranges']

        popFitnessAll = []

        for trange in tranges:
            popFitnessAll.append([min(np.exp(abs(v['target'] - simData['popRates'][k]['%d_%d'%(trange[0], trange[1])])/v['width']), maxFitness) 
                if simData['popRates'][k]['%d_%d'%(trange[0], trange[1])] >= v['min'] else maxFitness for k, v in pops.items()])
        
        popFitness = np.mean(np.array(popFitnessAll), axis=0)
        
        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f' % (p, np.mean(list(simData['popRates'][p].values())), popFitness[i]) for i,p in enumerate(pops)])
        print('  ' + popInfo)

        return fitness
        
    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.method = 'optuna'

    b.optimCfg = {
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'maxFitness': fitnessFuncArgs['maxFitness'],
        'maxiters':     1e6,    #    Maximum number of iterations (1 iteration = 1 function evaluation)
        'maxtime':      None,    #    Maximum time allowed, in seconds
        'maxiter_wait': 60,
        'time_sleep': 150,
        'popsize': 1  # unused - run with mpi 
    }

    return b



# ----------------------------------------------------------------------------------------------
# Optuna optimization
# ----------------------------------------------------------------------------------------------
def optunaRatesLayersWmat():

    # from prev
    import json
    #with open('data/v34_batch25/trial_2142/trial_2142_cfg.json', 'rb') as f:
    #    cfgLoad = json.load(f)['simConfig']

    with open('data/v34_batch23/trial_1937/trial_1937_cfg.json', 'rb') as f:
        cfgLoad = json.load(f)['simConfig']


    # --------------------------------------------------------
    # parameters
    params = specs.ODict()

    scaleLow = 0.75
    scaleHigh = 1.25

    #scaleLow2 = 0.1
    #scaleHigh2 = 10.0

    scaleLow2 = 0.5
    scaleHigh2 = 2.0

    # import pickle
    # with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    # wmat = connData['wmat']
    wmat = cfgLoad['wmat']    

    weightsScale = [['IT2', 'PV2'],
                    ['IT2', 'SOM2'], 
                    ['IT3', 'PV2'],
                    ['IT3', 'SOM2'],
                    ['PV2', 'PV2'],
                    ['PV2', 'VIP2'],
                    ['PV3', 'PV2'],
                    ['PV3', 'VIP2'],
                    ['SOM2', 'PV2'],
                    ['SOM2', 'VIP2'],
                    ['SOM3', 'PV2'],
                    ['SOM3', 'VIP2'],
                    ['VIP2', 'SOM2'],
                    ['VIP3', 'SOM2'],
                    ['IT2', 'SOM3'], 
                    ['IT3', 'SOM3'],
                    ['VIP2', 'SOM3'],
                    ['VIP3', 'SOM3'],
                    
                    ['IT6', 'PV6'], 
                    ['IT6', 'SOM6'], 
                    ['IT6', 'VIP6'], 
                    ['IT6', 'NGF6']]
                    

    for ws in weightsScale:
        params[('wmat', ws[0], ws[1])] = [wmat[ws[0]][ws[1]] * scaleLow, wmat[ws[0]][ws[1]] * scaleHigh]

    weightsScale2 = []

    for ws in weightsScale:
        params[('wmat', ws[0], ws[1])] = [wmat[ws[0]][ws[1]] * scaleLow2, wmat[ws[0]][ws[1]] * scaleHigh2]



    groupedParams = []

    # --------------------------------------------------------
    # initial config
    initCfg = {}
    initCfg = {}
    initCfg['duration'] = 2500
    initCfg['printPopAvgRates'] = [[1500, 1750], [1750, 2000], [2000, 2250], [2250, 2500]]
    initCfg['dt'] = 0.05

    initCfg['scaleDensity'] = 1.0

    # plotting and saving params
    initCfg[('analysis','plotRaster','markerSize')] = 10

    initCfg[('analysis','plotRaster','timeRange')] = [1500, 2500]
    initCfg[('analysis', 'plotTraces', 'timeRange')] = [1500, 2500]
    initCfg[('analysis', 'plotTraces', 'oneFigPer')] = 'trace'
    initCfg['recordLFP'] = None
    initCfg[('analysis', 'plotLFP')] = False

    initCfg['saveCellSecs'] = False
    initCfg['saveCellConns'] = False


    # from prev - best of 50% cell density
    updateParams = ['EEGain', 'EIGain', 'IEGain', 'IIGain',
                    ('EICellTypeGain', 'PV'), ('EICellTypeGain', 'SOM'), ('EICellTypeGain', 'VIP'), ('EICellTypeGain', 'NGF'),
                    ('IECellTypeGain', 'PV'), ('IECellTypeGain', 'SOM'), ('IECellTypeGain', 'VIP'), ('IECellTypeGain', 'NGF'),
                    ('EILayerGain', '1'), ('IILayerGain', '1'),
                    ('EELayerGain', '2'), ('EILayerGain', '2'),  ('IELayerGain', '2'), ('IILayerGain', '2'), 
                    ('EELayerGain', '3'), ('EILayerGain', '3'), ('IELayerGain', '3'), ('IILayerGain', '3'), 
                    ('EELayerGain', '4'), ('EILayerGain', '4'), ('IELayerGain', '4'), ('IILayerGain', '4'), 
                    ('EELayerGain', '5A'), ('EILayerGain', '5A'), ('IELayerGain', '5A'), ('IILayerGain', '5A'), 
                    ('EELayerGain', '5B'), ('EILayerGain', '5B'), ('IELayerGain', '5B'), ('IILayerGain', '5B'), 
                    ('EELayerGain', '6'), ('EILayerGain', '6'), ('IELayerGain', '6'), ('IILayerGain', '6')] 

    for p in updateParams:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad[p]})

    # good thal params for 100% cell density 
    updateParams2 = ['thalamoCorticalGain', 'intraThalamicGain', 'EbkgThalamicGain', 'IbkgThalamicGain']

    for p in updateParams2:
        if isinstance(p, tuple):
            initCfg.update({p: cfgLoad[p[0]][p[1]]})
        else:
            initCfg.update({p: cfgLoad[p]})

    print(initCfg)


    # --------------------------------------------------------
    # fitness function
    fitnessFuncArgs = {}
    pops = {}
    
    ## Exc pops
    Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'TCM', 'HTC']  # all layers + thal + IC
    # Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'TC', 'TCM', 'HTC']  # all layers + thal + IC

    Etune = {'target': 5, 'width': 20, 'min': 0.05}
    #Etune = {'target': 5, 'width': 5, 'min': 0.5}
    
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    # Ipops = ['NGF1',                            # L1
    #         'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
    #         'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
    #         'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
    #         'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
    #         'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
    #         'PV6', 'SOM6', 'VIP6', 'NGF6',       # L6
    #         'IRE', 'IREM', 'TI']  # Thal 
    Ipops = ['NGF1',  
            'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
            'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
            'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
            'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A 
            'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
            'PV6', 'SOM6', 'VIP6', 'NGF6',  # L6
            'IRE', 'IREM', 'TI']  # Thal 

    Itune = {'target': 10, 'width': 30, 'min': 0.05}
    #Itune = {'target': 10, 'width': 15, 'min': 0.5}

    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 2000
    fitnessFuncArgs['tranges'] = initCfg['printPopAvgRates']


    def fitnessFunc(simData, **kwargs):
        import numpy as np
        pops = kwargs['pops']
        maxFitness = kwargs['maxFitness']
        tranges = kwargs['tranges']

        popFitnessAll = []

        for trange in tranges:
            popFitnessAll.append([min(np.exp(abs(v['target'] - simData['popRates'][k]['%d_%d'%(trange[0], trange[1])])/v['width']), maxFitness) 
                if simData['popRates'][k]['%d_%d'%(trange[0], trange[1])] >= v['min'] else maxFitness for k, v in pops.items()])
        
        popFitness = np.mean(np.array(popFitnessAll), axis=0)
        
        fitness = np.mean(popFitness)

        popInfo = '; '.join(['%s rate=%.1f fit=%1.f' % (p, np.mean(list(simData['popRates'][p].values())), popFitness[i]) for i,p in enumerate(pops)])
        print('  ' + popInfo)

        return fitness
        
    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    b.method = 'optuna'

    b.optimCfg = {
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'maxFitness': fitnessFuncArgs['maxFitness'],
        'maxiters':     1e6,    #    Maximum number of iterations (1 iteration = 1 function evaluation)
        'maxtime':      None,    #    Maximum time allowed, in seconds
        'maxiter_wait': 60,
        'time_sleep': 150,
        'popsize': 1  # unused - run with mpi 
    }

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin'):
    if type=='mpi_bulletin':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'init_cell.py', 
            'skip': True}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'nodes': 1,
            'coresPerNode': 96,
            'script': 'init.py',
            'mpiCommand': 'mpirun',
            'skip': True}

    elif type=='hpc_torque':
        b.runCfg = {'type': 'hpc_torque',
             'script': 'init.py',
             'nodes': 3,
             'ppn': 8,
             'walltime': "12:00:00",
             'queueName': 'longerq',
             'sleepInterval': 5,
             'skip': True}

    elif type=='hpc_slurm_comet':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'shs100', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403'
            #'reservation': 'salva1',
            'walltime': '6:00:00',
            'nodes': 4,
            'coresPerNode': 24,  # comet=24, bridges=28
            'email': 'salvadordura@gmail.com',
            'folder': '/home/salvadord/m1/sim/',  # comet='/salvadord', bridges='/salvi82'
            'script': 'init.py', 
            'mpiCommand': 'ibrun', # comet='ibrun', bridges='mpirun'
            'skipCustom': '_raster.png'}

    elif type=='hpc_slurm_gcp':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'default', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403', gcp='default'
            'walltime': '24:00:00', #'48:00:00',
            'nodes': 1,
            'coresPerNode': 80,  # comet=24, bridges=28, gcp=32
            'email': 'salvadordura@gmail.com',
            'folder': '/home/ext_salvadordura_gmail_com/A1_layers/',  # comet,gcp='/salvadord', bridges='/salvi82'
            'script': 'init.py',
            'mpiCommand': 'mpirun', # comet='ibrun', bridges,gcp='mpirun' 
            'nrnCommand': 'nrniv -mpi -python', #'python3',
            'skipCustom': '_raster.png'}
            #'custom': '#SBATCH --exclude=compute[17-64000]'} # only use first 16 nodes (non-preemptible for long runs )
            # --nodelist=compute1


    elif type=='hpc_slurm_bridges':
        b.runCfg = {'type': 'hpc_slurm', 
            'allocation': 'ib4iflp', # bridges='ib4iflp', comet m1='shs100', comet nsg='csd403'
            'walltime': '06:00:00',
            'nodes': 2,
            'coresPerNode': 28,  # comet=24, bridges=28
            'email': 'salvadordura@gmail.com',
            'folder': '/home/salvi82/m1/sim/',  # comet='/salvadord', bridges='/salvi82'
            'script': 'init.py', 
            'mpiCommand': 'mpirun', # comet='ibrun', bridges='mpirun'
            'skip': True}



# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

if __name__ == '__main__':

    cellTypes = ['IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'HTC', 'IRE', 'TI']

    b = custom_spont('data/v34_batch25/trial_2142/trial_2142_cfg.json')
    # b = optunaRatesLayersThalL12345A5B6()
    # b = optunaRatesLayersWmat()

    # b = bkgWeights(pops = cellTypes, weights = list(np.arange(1,100)))
    # b = fIcurve(pops=['ITS4']) 

    b.batchLabel = 'v34_batch68' 
    b.saveFolder = 'data/'+b.batchLabel

    setRunCfg(b, 'hpc_slurm_gcp') #'hpc_slurm_gcp') #'mpi_bulletin') #'hpc_slurm_gcp')
    b.run() # run batch

