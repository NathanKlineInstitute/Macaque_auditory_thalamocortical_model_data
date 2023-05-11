"""
netParams.py 

High-level specifications for A1 network model using NetPyNE

Contributors: ericaygriffith@gmail.com, salvadordura@gmail.com
"""

from netpyne import specs
import pickle, json

netParams = specs.NetParams()   # object of class NetParams to store the network parameters

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg


#------------------------------------------------------------------------------
# VERSION 
#------------------------------------------------------------------------------
netParams.version = 34

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General network parameters
#------------------------------------------------------------------------------

netParams.scale = cfg.scale # Scale factor for number of cells # NOT DEFINED YET! 3/11/19 # How is this different than scaleDensity? 
netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder' # cylindrical (column-like) volume

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.scaleConnWeight = 1.0 # Connection weight scale factor (default if no model specified)
netParams.scaleConnWeightModels = { 'HH_reduced': 1.0, 'HH_full': 1.0} #scale conn weight factor for each cell model
netParams.scaleConnWeightNetStims = 1.0 #0.5  # scale conn weight factor for NetStims
netParams.defaultThreshold = 0.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 500.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)

#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------

Etypes = ['IT', 'ITS4', 'PT', 'CT']
Itypes = ['PV', 'SOM', 'VIP', 'NGF']
cellModels = ['HH_reduced', 'HH_full'] # List of cell models

# II: 100-950, IV: 950-1250, V: 1250-1550, VI: 1550-2000 
layer = {'1': [0.00, 0.05], '2': [0.05, 0.08], '3': [0.08, 0.475], '4': [0.475, 0.625], '5A': [0.625, 0.667], '5B': [0.667, 0.775], '6': [0.775, 1], 'thal': [1.2, 1.4], 'cochlear': [1.6, 1.8]}  # normalized layer boundaries  

layerGroups = { '1-3': [layer['1'][0], layer['3'][1]],  # L1-3
                '4': layer['4'],                      # L4
                '5': [layer['5A'][0], layer['5B'][1]],  # L5A-5B
                '6': layer['6']}                        # L6

# add layer border correction ??
#netParams.correctBorder = {'threshold': [cfg.correctBorderThreshold, cfg.correctBorderThreshold, cfg.correctBorderThreshold], 
#                        'yborders': [layer['2'][0], layer['5A'][0], layer['6'][0], layer['6'][1]]}  # correct conn border effect


#------------------------------------------------------------------------------
## Load cell rules previously saved using netpyne format (DOES NOT INCLUDE VIP, NGF and spiny stellate)
## include conditions ('conds') for each cellRule
cellParamLabels = ['IT2_reduced', 'IT3_reduced', 'ITP4_reduced', 'ITS4_reduced',
                    'IT5A_reduced', 'CT5A_reduced', 'IT5B_reduced',
                    'PT5B_reduced', 'CT5B_reduced', 'IT6_reduced', 'CT6_reduced',
                    'PV_reduced', 'SOM_reduced', 'VIP_reduced', 'NGF_reduced',
                    'RE_reduced', 'TC_reduced', 'HTC_reduced', 'TI_reduced']

for ruleLabel in cellParamLabels:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='cells/' + ruleLabel + '_cellParams.json')  # Load cellParams for each of the above cell subtype

# change weightNorm 
for k in cfg.weightNormScaling:
    for sec in netParams.cellParams[k]['secs'].values():
        for i in range(len(sec['weightNorm'])):
            sec['weightNorm'][i] *= cfg.weightNormScaling[k]



#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------

## load densities
with open('cells/cellDensity.pkl', 'rb') as fileObj: density = pickle.load(fileObj)['density']
density = {k: [x * cfg.scaleDensity for x in v] for k,v in density.items()} # Scale densities 

# ### LAYER 1:
netParams.popParams['NGF1'] = {'cellType': 'NGF', 'cellModel': 'HH_reduced','ynormRange': layer['1'],   'density': density[('A1','nonVIP')][0]}

### LAYER 2:
netParams.popParams['IT2'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',  'ynormRange': layer['2'],   'density': density[('A1','E')][1]}     #
netParams.popParams['SOM2'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','SOM')][1]}   
netParams.popParams['PV2'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','PV')][1]}    
netParams.popParams['VIP2'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','VIP')][1]}
netParams.popParams['NGF2'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['2'],   'density': density[('A1','nonVIP')][1]}

### LAYER 3:
netParams.popParams['IT3'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',  'ynormRange': layer['3'],   'density': density[('A1','E')][1]} 
netParams.popParams['SOM3'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','SOM')][1]} 
netParams.popParams['PV3'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','PV')][1]} 
netParams.popParams['VIP3'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','VIP')][1]} 
netParams.popParams['NGF3'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['3'],   'density': density[('A1','nonVIP')][1]}


### LAYER 4: 
netParams.popParams['ITP4'] =	 {'cellType': 'IT', 'cellModel': 'HH_reduced',  'ynormRange': layer['4'],   'density': 0.5*density[('A1','E')][2]}      ## CHANGE DENSITY #
netParams.popParams['ITS4'] =	 {'cellType': 'IT', 'cellModel': 'HH_reduced', 'ynormRange': layer['4'],  'density': 0.5*density[('A1','E')][2]}  
netParams.popParams['SOM4'] = 	 {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],  'density': density[('A1','SOM')][2]}
netParams.popParams['PV4'] = 	 {'cellType': 'PV', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],   'density': density[('A1','PV')][2]}
netParams.popParams['VIP4'] =	 {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],  'density': density[('A1','VIP')][2]}
netParams.popParams['NGF4'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['4'],  'density': density[('A1','nonVIP')][2]}

# # ### LAYER 5A: 
netParams.popParams['IT5A'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5A'], 	'density': 0.5*density[('A1','E')][3]}      
netParams.popParams['CT5A'] =     {'cellType': 'CT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5A'],   'density': 0.5*density[('A1','E')][3]}   
netParams.popParams['SOM5A'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],	'density': density[('A1','SOM')][3]}          
netParams.popParams['PV5A'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],	'density': density[('A1','PV')][3]}         
netParams.popParams['VIP5A'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],   'density': density[('A1','VIP')][3]}
netParams.popParams['NGF5A'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',    'ynormRange': layer['5A'],   'density': density[('A1','nonVIP')][3]}

### LAYER 5B: 
netParams.popParams['IT5B'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5B'], 	'density': (1/3)*density[('A1','E')][4]}  
netParams.popParams['CT5B'] =     {'cellType': 'CT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5B'],   'density': (1/3)*density[('A1','E')][4]}  
netParams.popParams['PT5B'] =     {'cellType': 'PT',  'cellModel': 'HH_reduced',   'ynormRange': layer['5B'], 	'density': (1/3)*density[('A1','E')][4]}  
netParams.popParams['SOM5B'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],   'density': density[('A1', 'SOM')][4]}
netParams.popParams['PV5B'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],	'density': density[('A1','PV')][4]}     
netParams.popParams['VIP5B'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],   'density': density[('A1','VIP')][4]}
netParams.popParams['NGF5B'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',    'ynormRange': layer['5B'],   'density': density[('A1','nonVIP')][4]}

# # ### LAYER 6:
netParams.popParams['IT6'] =     {'cellType': 'IT',  'cellModel': 'HH_reduced',  'ynormRange': layer['6'],   'density': 0.5*density[('A1','E')][5]}  
netParams.popParams['CT6'] =     {'cellType': 'CT',  'cellModel': 'HH_reduced',  'ynormRange': layer['6'],   'density': 0.5*density[('A1','E')][5]} 
netParams.popParams['SOM6'] =    {'cellType': 'SOM', 'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','SOM')][5]}   
netParams.popParams['PV6'] =     {'cellType': 'PV',  'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','PV')][5]}     
netParams.popParams['VIP6'] =    {'cellType': 'VIP', 'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','VIP')][5]}
netParams.popParams['NGF6'] =    {'cellType': 'NGF', 'cellModel': 'HH_reduced',   'ynormRange': layer['6'],   'density': density[('A1','nonVIP')][5]}


### THALAMIC POPULATIONS (from prev model)
thalDensity = density[('A1','PV')][2] * 1.25  # temporary estimate (from prev model)

netParams.popParams['TC'] =     {'cellType': 'TC',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.75*thalDensity}  
netParams.popParams['TCM'] =    {'cellType': 'TC',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': thalDensity} 
netParams.popParams['HTC'] =    {'cellType': 'HTC', 'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.25*thalDensity}   
netParams.popParams['IRE'] =    {'cellType': 'RE',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': thalDensity}     
netParams.popParams['IREM'] =   {'cellType': 'RE', 'cellModel': 'HH_reduced',   'ynormRange': layer['thal'],   'density': thalDensity}
netParams.popParams['TI'] =     {'cellType': 'TI',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.33 * thalDensity} ## Winer & Larue 1996; Huang et al 1999 
netParams.popParams['TIM'] =    {'cellType': 'TI',  'cellModel': 'HH_reduced',  'ynormRange': layer['thal'],   'density': 0.33 * thalDensity} ## Winer & Larue 1996; Huang et al 1999 


if cfg.singleCellPops:
    for pop in netParams.popParams.values(): pop['numCells'] = 1

## List of E and I pops to use later on
Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B' , 'PT5B', 'IT6', 'CT6']  # all layers

Ipops = ['NGF1',                            # L1
        'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
        'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
        'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
        'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
        'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
        'PV6', 'SOM6', 'VIP6', 'NGF6']      # L6 



#------------------------------------------------------------------------------
# Synaptic mechanism parameters
#------------------------------------------------------------------------------

### From M1 detailed netParams.py 
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3*cfg.AMPATau2Factor, 'e': 0}
netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93} 
netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.07, 'tau2': 18.2, 'e': -80}
netParams.synMechParams['GABAA_VIP'] = {'mod':'MyExp2SynBB', 'tau1': 0.3, 'tau2': 6.4, 'e': -80}  # Pi et al 2013
netParams.synMechParams['GABAASlow'] = {'mod': 'MyExp2SynBB','tau1': 2, 'tau2': 100, 'e': -80}
netParams.synMechParams['GABAASlowSlow'] = {'mod': 'MyExp2SynBB', 'tau1': 200, 'tau2': 400, 'e': -80}

ESynMech = ['AMPA', 'NMDA']
SOMESynMech = ['GABAASlow','GABAB']
SOMISynMech = ['GABAASlow']
PVSynMech = ['GABAA']
VIPSynMech = ['GABAA_VIP']
NGFSynMech = ['GABAA', 'GABAB']


#------------------------------------------------------------------------------
# Local connectivity parameters
#------------------------------------------------------------------------------

## load data from conn pre-processing file
with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
pmat = connData['pmat']
lmat = connData['lmat']
wmat = connData['wmat']
bins = connData['bins']
connDataSource = connData['connDataSource']

wmat = cfg.wmat

layerGainLabels = ['1', '2', '3', '4', '5A', '5B', '6']

#------------------------------------------------------------------------------
## E -> E
if cfg.addConn and cfg.EEGain > 0.0:
    for pre in Epops:
        for post in Epops:
            for l in layerGainLabels:  # used to tune each layer group independently
                if connDataSource['E->E/I'] in ['Allen_V1', 'Allen_custom']:
                    prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])
                else:
                    prob = pmat[pre][post]
                netParams.connParams['EE_'+pre+'_'+post+'_'+l] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post, 'ynorm': layer[l]},
                    'synMech': ESynMech,
                    'probability': prob,
                    'weight': wmat[pre][post] * cfg.EEGain * cfg.EELayerGain[l], 
                    'synMechWeightFactor': cfg.synWeightFractionEE,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': 1,
                    'sec': 'dend_all'}
                    

#------------------------------------------------------------------------------
## E -> I
if cfg.addConn and cfg.EIGain > 0.0:
    for pre in Epops:
        for post in Ipops:
            for postType in Itypes:
                if postType in post: # only create rule if celltype matches pop
                    for l in layerGainLabels:  # used to tune each layer group independently
                        if connDataSource['E->E/I'] in ['Allen_V1', 'Allen_custom']:
                            prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])
                        else:
                            prob = pmat[pre][post]
                        
                        if 'NGF' in post:
                            synWeightFactor = cfg.synWeightFractionENGF
                        else:
                            synWeightFactor = cfg.synWeightFractionEI
                        netParams.connParams['EI_'+pre+'_'+post+'_'+postType+'_'+l] = { 
                            'preConds': {'pop': pre}, 
                            'postConds': {'pop': post, 'cellType': postType, 'ynorm': layer[l]},
                            'synMech': ESynMech,
                            'probability': prob,
                            'weight': wmat[pre][post] * cfg.EIGain * cfg.EICellTypeGain[postType] * cfg.EILayerGain[l], 
                            'synMechWeightFactor': synWeightFactor,
                            'delay': 'defaultDelay+dist_3D/propVelocity',
                            'synsPerConn': 1,
                            'sec': 'proximal'}
                

#------------------------------------------------------------------------------
## I -> E
if cfg.addConn and cfg.IEGain > 0.0:

    if connDataSource['I->E/I'] == 'Allen_custom':

        ESynMech = ['AMPA', 'NMDA']
        SOMESynMech = ['GABAASlow','GABAB']
        SOMISynMech = ['GABAASlow']
        PVSynMech = ['GABAA']
        VIPSynMech = ['GABAA_VIP']
        NGFSynMech = ['GABAA', 'GABAB']

        for pre in Ipops:
            for preType in Itypes:
                if preType in pre:  # only create rule if celltype matches pop
                    for post in Epops:
                        for l in layerGainLabels:  # used to tune each layer group independently
                            
                            prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])
                            
                            if 'SOM' in pre:
                                synMech = SOMESynMech
                            elif 'PV' in pre:
                                synMech = PVSynMech
                            elif 'VIP' in pre:
                                synMech = VIPSynMech
                            elif 'NGF' in pre:
                                synMech = NGFSynMech

                            netParams.connParams['IE_'+pre+'_'+preType+'_'+post+'_'+l] = { 
                                'preConds': {'pop': pre}, 
                                'postConds': {'pop': post, 'ynorm': layer[l]},
                                'synMech': synMech,
                                'probability': prob,
                                'weight': wmat[pre][post] * cfg.IEGain * cfg.IECellTypeGain[preType] * cfg.IELayerGain[l], 
                                'synMechWeightFactor': cfg.synWeightFractionEI,
                                'delay': 'defaultDelay+dist_3D/propVelocity',
                                'synsPerConn': 1,
                                'sec': 'proximal'}
                    

#------------------------------------------------------------------------------
## I -> I
if cfg.addConn and cfg.IIGain > 0.0:

    if connDataSource['I->E/I'] == 'Allen_custom':

        for pre in Ipops:
            for post in Ipops:
                for l in layerGainLabels: 
                    
                    prob = '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post])

                    if 'SOM' in pre:
                        synMech = SOMISynMech
                    elif 'PV' in pre:
                        synMech = PVSynMech
                    elif 'VIP' in pre:
                        synMech = VIPSynMech
                    elif 'NGF' in pre:
                        synMech = NGFSynMech

                    netParams.connParams['II_'+pre+'_'+post+'_'+l] = { 
                        'preConds': {'pop': pre}, 
                        'postConds': {'pop': post,  'ynorm': layer[l]},
                        'synMech': synMech,
                        'probability': prob,
                        'weight': wmat[pre][post] * cfg.IIGain * cfg.IILayerGain[l], 
                        'synMechWeightFactor': cfg.synWeightFractionII,
                        'delay': 'defaultDelay+dist_3D/propVelocity',
                        'synsPerConn': 1,
                        'sec': 'proximal'}
                        

#------------------------------------------------------------------------------
# Thalamic connectivity parameters
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
## Intrathalamic 

TEpops = ['TC', 'TCM', 'HTC']
TIpops = ['IRE', 'IREM', 'TI', 'TIM']

if cfg.addConn and cfg.addIntraThalamicConn:
    for pre in TEpops+TIpops:
        for post in TEpops+TIpops:
            if post in pmat[pre]:
                # for syns use ESynMech, SOMESynMech and SOMISynMech 
                if pre in TEpops:     # E->E/I
                    syn = ESynMech
                    synWeightFactor = cfg.synWeightFractionEE
                elif post in TEpops:  # I->E
                    syn = SOMESynMech
                    synWeightFactor = cfg.synWeightFractionIE
                else:                  # I->I
                    syn = SOMISynMech
                    synWeightFactor = [1.0]
                    
                netParams.connParams['ITh_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'synMech': syn,
                    'probability': pmat[pre][post],
                    'weight': wmat[pre][post] * cfg.intraThalamicGain, 
                    'synMechWeightFactor': synWeightFactor,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': 1,
                    'sec': 'soma'}  


#------------------------------------------------------------------------------
## Corticothalamic 
if cfg.addConn and cfg.addCorticoThalamicConn:
    for pre in Epops:
        for post in TEpops+TIpops:
            if post in pmat[pre]:
                netParams.connParams['CxTh_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'synMech': ESynMech,
                    'probability': pmat[pre][post],
                    'weight': wmat[pre][post] * cfg.corticoThalamicGain, 
                    'synMechWeightFactor': cfg.synWeightFractionEE,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': 1,
                    'sec': 'soma'}  

#------------------------------------------------------------------------------
## Thalamocortical 
if cfg.addConn and cfg.addThalamoCorticalConn:
    for pre in TEpops+TIpops:
        for post in Epops+Ipops:
            if post in pmat[pre]:
                # for syns use ESynMech, SOMESynMech and SOMISynMech 
                if pre in TEpops:     # E->E/I
                    syn = ESynMech
                    synWeightFactor = cfg.synWeightFractionEE
                elif post in Epops:  # I->E
                    syn = SOMESynMech
                    synWeightFactor = cfg.synWeightFractionIE
                else:                  # I->I
                    syn = SOMISynMech
                    synWeightFactor = [1.0]

                netParams.connParams['ThCx_'+pre+'_'+post] = { 
                    'preConds': {'pop': pre}, 
                    'postConds': {'pop': post},
                    'synMech': syn,
                    'probability': '%f * exp(-dist_2D/%f)' % (pmat[pre][post], lmat[pre][post]),
                    'weight': wmat[pre][post] * cfg.thalamoCorticalGain, 
                    'synMechWeightFactor': synWeightFactor,
                    'delay': 'defaultDelay+dist_3D/propVelocity',
                    'synsPerConn': 1,
                    'sec': 'soma'}  


#------------------------------------------------------------------------------
# Subcellular connectivity (synaptic distributions)
#------------------------------------------------------------------------------  

# Set target sections (somatodendritic distribution of synapses)
# From Billeh 2019 (Allen V1) (fig 4F) and Tremblay 2016 (fig 3)

if cfg.addSubConn:
    #------------------------------------------------------------------------------
    # E -> E2/3,4: soma,dendrites <200um
    netParams.subConnParams['E->E2,3,4'] = {
        'preConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']}, 
        'postConds': {'pops': ['IT2', 'IT3', 'ITP4', 'ITS4']},
        'sec': 'proximal',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    # E -> E5,6: soma,dendrites (all)
    netParams.subConnParams['E->E5,6'] = {
        'preConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']}, 
        'postConds': {'pops': ['IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']},
        'sec': 'all',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'}
        
    #------------------------------------------------------------------------------
    # E -> I: soma, dendrite (all)
    netParams.subConnParams['E->I'] = {
        'preConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']}, 
        'postConds': {'cellType': ['PV','SOM','NGF', 'VIP']},
        'sec': 'all',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    # NGF1 -> E: apic_tuft
    netParams.subConnParams['NGF1->E'] = {
        'preConds': {'pops': ['NGF1']}, 
        'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
        'sec': 'apic_tuft',
        'groupSynMechs': NGFSynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    # NGF2,3,4 -> E2,3,4: apic_trunk
    netParams.subConnParams['NGF2,3,4->E2,3,4'] = {
        'preConds': {'pops': ['NGF2', 'NGF3', 'NGF4']}, 
        'postConds': {'pops': ['IT2', 'IT3', 'ITP4', 'ITS4']},
        'sec': 'apic_trunk',
        'groupSynMechs': NGFSynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    # NGF2,3,4 -> E5,6: apic_uppertrunk
    netParams.subConnParams['NGF2,3,4->E5,6'] = {
        'preConds': {'pops': ['NGF2', 'NGF3', 'NGF4']}, 
        'postConds': {'pops': ['IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']},
        'sec': 'apic_uppertrunk',
        'groupSynMechs': NGFSynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    # NGF5,6 -> E5,6: apic_lowerrunk
    netParams.subConnParams['NGF5,6->E5,6'] = {
        'preConds': {'pops': ['NGF5A', 'NGF5B', 'NGF6']}, 
        'postConds': {'pops': ['IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']},
        'sec': 'apic_lowertrunk',
        'groupSynMechs': NGFSynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    #  SOM -> E: all_dend (not close to soma)
    netParams.subConnParams['SOM->E'] = {
        'preConds': {'cellType': ['SOM']}, 
        'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
        'sec': 'dend_all',
        'groupSynMechs': SOMESynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    #  PV -> E: proximal
    netParams.subConnParams['PV->E'] = {
        'preConds': {'cellType': ['PV']}, 
        'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
        'sec': 'proximal',
        'groupSynMechs': PVSynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    #  TC -> E: proximal
    netParams.subConnParams['TC->E'] = {
        'preConds': {'cellType': ['TC', 'HTC']}, 
        'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
        'sec': 'proximal',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'} 

    #------------------------------------------------------------------------------
    #  TCM -> E: apical
    netParams.subConnParams['TCM->E'] = {
        'preConds': {'cellType': ['TCM']}, 
        'postConds': {'cellType': ['IT', 'ITS4', 'PT', 'CT']},
        'sec': 'apic',
        'groupSynMechs': ESynMech, 
        'density': 'uniform'}
        

#------------------------------------------------------------------------------
# Background inputs 
#------------------------------------------------------------------------------  
if cfg.addBkgConn:
    # add bkg sources for E and I cells
    netParams.stimSourceParams['excBkg'] = {'type': 'NetStim', 'start': cfg.startBkg, 'rate': cfg.rateBkg['exc'], 'noise': cfg.noiseBkg, 'number': 1e9}
    netParams.stimSourceParams['inhBkg'] = {'type': 'NetStim', 'start': cfg.startBkg, 'rate': cfg.rateBkg['inh'], 'noise': cfg.noiseBkg, 'number': 1e9}
    
    if cfg.cochlearThalInput:
        from input import cochlearInputSpikes
        numCochlearCells = cfg.cochlearThalInput['numCells']
        cochlearSpkTimes = cochlearInputSpikes(numCells = numCochlearCells,
                                               duration = cfg.duration,
                                               freqRange = cfg.cochlearThalInput['freqRange'],
                                               toneFreq=cfg.cochlearThalInput['toneFreq'],
                                               loudnessDBs=cfg.cochlearThalInput['loudnessDBs'])
                                              
        netParams.popParams['cochlea'] = {'cellModel': 'VecStim', 'numCells': numCochlearCells, 'spkTimes': cochlearSpkTimes, 'ynormRange': layer['cochlear']}

    if cfg.ICThalInput:
        # load file with IC output rates
        from scipy.io import loadmat
        import numpy as np

        data = loadmat(cfg.ICThalInput['file'])
        fs = data['RsFs'][0][0]
        ICrates = data['BE_sout_population'].tolist()
        ICtimes = list(np.arange(0, cfg.duration, 1000./fs))  # list with times to set each time-dep rate
        
        ICrates = ICrates * 4 # 200 cells
        
        numCells = len(ICrates)

        # Option 1: create population of DynamicNetStims with time-varying rates
        #netParams.popParams['IC'] = {'cellModel': 'DynamicNetStim', 'numCells': numCells, 'ynormRange': layer['cochlear'],
        #    'dynamicRates': {'rates': ICrates, 'times': ICtimes}}

        # Option 2:
        from input import inh_poisson_generator
        
        maxLen = min(len(ICrates[0]), len(ICtimes))
        spkTimes = [[x+cfg.ICThalInput['startTime'] for x in inh_poisson_generator(ICrates[i][:maxLen], ICtimes[:maxLen], cfg.duration, cfg.ICThalInput['seed']+i)] for i in range(len(ICrates))]
        netParams.popParams['IC'] = {'cellModel': 'VecStim', 'numCells': numCells, 'ynormRange': layer['cochlear'],
            'spkTimes': spkTimes}


    # excBkg/I -> thalamus + cortex
    with open('cells/bkgWeightPops.json', 'r') as f:
        weightBkg = json.load(f)
    pops = list(cfg.allpops)
    pops.remove('IC')

    for pop in ['TC', 'TCM', 'HTC']:
        weightBkg[pop] *= cfg.EbkgThalamicGain 

    for pop in ['IRE', 'IREM', 'TI', 'TIM']:
        weightBkg[pop] *= cfg.IbkgThalamicGain 


    for pop in pops:
        netParams.stimTargetParams['excBkg->'+pop] =  {
            'source': 'excBkg', 
            'conds': {'pop': pop},
            'sec': 'apic', 
            'loc': 0.5,
            'synMech': ESynMech,
            'weight': weightBkg[pop],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}

        netParams.stimTargetParams['inhBkg->'+pop] =  {
            'source': 'inhBkg', 
            'conds': {'pop': pop},
            'sec': 'proximal',
            'loc': 0.5,
            'synMech': 'GABAA',
            'weight': weightBkg[pop],
            'delay': cfg.delayBkg}

    # cochlea -> thal
    if cfg.cochlearThalInput:
        netParams.connParams['cochlea->ThalE'] = { 
            'preConds': {'pop': 'cochlea'}, 
            'postConds': {'cellType': ['TC', 'HTC']},
            'sec': 'soma', 
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': cfg.probInput['ThalE'], 
            'weight': cfg.weightInput['ThalE'],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}
        
        netParams.connParams['cochlea->ThalI'] = { 
            'preConds': {'pop': 'cochlea'}, 
            'postConds': {'cellType': ['RE']},
            'sec': 'soma', 
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': cfg.probInput['ThalI'], 
            'weight': cfg.weightInput['ThalI'],
            'synMechWeightFactor': cfg.synWeightFractionEI,
            'delay': cfg.delayBkg}  

    # cochlea/IC -> thal
    if cfg.ICThalInput:
        netParams.connParams['IC->ThalE'] = { 
            'preConds': {'pop': 'IC'}, 
            'postConds': {'cellType': ['TC', 'HTC']},
            'sec': 'soma', 
            'loc': 0.5,
            'synMech': ESynMech,
            'probability': cfg.ICThalInput['probE'],
            'weight': cfg.ICThalInput['weightE'],
            'synMechWeightFactor': cfg.synWeightFractionEE,
            'delay': cfg.delayBkg}
        
        netParams.connParams['IC->ThalI'] = { 
            'preConds': {'pop': 'IC'}, 
            'postConds': {'cellType': ['RE', 'TI']},
            'sec': 'soma', 
            'loc': 0.5,
            'synMech': 'GABAA',
            'probability': cfg.ICThalInput['probI'],
            'weight': cfg.ICThalInput['weightI'],
            'delay': cfg.delayBkg}  


#------------------------------------------------------------------------------
# NetStim inputs (to simulate short external stimuli; not bkg)
#------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith('NetStim')]:
        params = getattr(cfg, key, None)
        [pop, ynorm, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay] = \
        [params[s] for s in ['pop', 'ynorm', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay']] 

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}
        
        if not isinstance(pop, list):
            pop = [pop]

        for eachPop in pop:
            # connect stim source to target 
            print(key, eachPop)
            netParams.stimTargetParams[key+'_'+eachPop] =  {
                'source': key, 
                'conds': {'pop': eachPop, 'ynorm': ynorm},
                'sec': sec, 
                'loc': loc,
                'synMech': synMech,
                'weight': weight,
                'synMechWeightFactor': synMechWeightFactor,
                'delay': delay}

#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------

netParams.description = """
v7 - Added template for connectivity
v8 - Added cell types
v9 - Added local connectivity
v10 - Added thalamic populations from prev model
v11 - Added thalamic conn from prev model
v12 - Added CT cells to L5B
v13 - Added CT cells to L5A
v14 - Fixed L5A & L5B E cell densities + added CT5A & CT5B to 'Epops'
v15 - Added cortical and thalamic conn to CT5A and CT5B 
v16 - Updated multiple cell types
v17 - Changed NGF -> I prob from strong (1.0) to weak (0.35)
v18 - Fixed bug in VIP cell morphology
v19 - Added in 2-compartment thalamic interneuron model 
v20 - Added TI conn and updated thal pop
v21 - Added exc+inh bkg inputs specific to each cell type
v22 - Made exc+inh bkg inputs specific to each pop; automated calculation
v23 - IE/II specific layer gains and simplified code (assume 'Allen_custom')
v24 - Fixed bug in IE/II specific layer gains
v25 - Fixed subconnparams TC->E and NGF1->E; made IC input deterministic
v26 - Changed NGF AMPA:NMDA ratio 
v27 - Split thalamic interneurons into core and matrix (TI and TIM)
v28 - Set recurrent TC->TC conn to 0
v29 - Added EI specific layer gains
v30 - Added EE specific layer gains; and split combined L1-3 gains into L1,L2,L3
v31 - Added EI postsyn-cell-type specific gains; update ITS4 and NGF
v32 - Added IE presyn-cell-type specific gains
v33 - Fixed bug in matrix thalamocortical conn (were very low)
v34 - Added missing conn from cortex to matrix thalamus IREM and TIM
"""