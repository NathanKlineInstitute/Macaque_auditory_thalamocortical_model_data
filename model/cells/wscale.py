"""
wscale.py

Code to analyze and calculate the scaling of weights as a function of input dendritic location

Contributors: salvadordura@gmail.com
"""

import utils
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from cycler import cycler
import sys, os

def calculateEPSPs(params, data, somaLabel='soma', stimRange=[3000,4000], syn='exc'):
    out = {}
    secs = [s for s in params[0]['values']]
    locs = [s for s in params[1]['values']]

    for key, d in data.items():
        cellLabel = d['simData']['V_soma'].keys()[0]
        vsoma = d['simData']['V_'+somaLabel][cellLabel]
        if syn == 'exc':
            epsp = max(vsoma[stimRange[0]:stimRange[1]]) - vsoma[stimRange[0]-1] # max voltage between stim time - baseline
        elif syn == 'inh':
            epsp = min(vsoma[stimRange[0]:stimRange[1]]) - vsoma[stimRange[0]-1] # min voltage between stim time - baseline
        seg = (d['paramValues'][0], d['paramValues'][1])
        weight = d['paramValues'][1]
        print(seg, weight, epsp, len(vsoma))
        out[tuple(seg)].append([weight, epsp])

    return out


def calculateEPSPsPops(params, data, somaLabel='soma', stimRange=[3000,4000], syn='exc'):
    out = {}

    pops = [p for p in params[2]['values']]
    secs = [s for s in params[0]['values']]
    locs = [s for s in params[1]['values']]

    for pop in pops: 
        out[pop] = {}
        for sec,loc in zip(secs,locs): out[pop][(sec,loc)] = []

    for key, d in data.items():
        
        cellLabel = list(d['V_soma'].keys())[0] # d['simData']['V_soma'].keys()[0]
        vsoma = d['V_'+somaLabel][cellLabel]  #d['simData']['V_'+somaLabel][cellLabel]
        if syn == 'exc':
            epsp = max(vsoma[stimRange[0]:stimRange[1]]) - vsoma[stimRange[0]-1] # max voltage between stim time - baseline
        elif syn == 'inh':
            epsp = min(vsoma[stimRange[0]:stimRange[1]]) - vsoma[stimRange[0]-1] # min voltage between stim time - baseline
        pop = d['paramValues'][2]
        seg = (d['paramValues'][0], d['paramValues'][1])
        weight = d['paramValues'][3]
        print(pop, seg, weight, epsp, len(vsoma))
        out[pop][tuple(seg)].append([weight, epsp])

    return out


def calculateWeightNorm(params, data, epspNorm=0.5, somaLabel='soma', stimRange=[3000,4000], savePath=None):
    epsp = calculateEPSPs(params, data, somaLabel=somaLabel, stimRange=stimRange)

    
    segs = [s for s in params[1]['values']]
    segs.sort()
    weightNorm = {}
    for seg in segs: weightNorm[seg[0]] = []  # empty list for each section
    for seg in segs:
        epspSeg = epsp[tuple(seg)]
        epspSeg.sort()
        x,y = zip(*epspSeg)
        f = interp1d(y,x,fill_value="extrapolate")
        w = f(epspNorm)
        wnorm = w / epspNorm
        weightNorm[seg[0]].append(wnorm)
        print('\n%s wscale = %.6f' % (str(seg), wnorm))

        if savePath:
            import pickle
            with open(savePath+'_weightNorm.pkl', 'wb') as fileObj:
                pickle.dump(weightNorm, fileObj)


def calculateWeightNormPops(params, data, epspNorm=0.5, somaLabel='soma', stimRange=[3000, 4000], savePath=None,
    popSaveLabels ={'IT2': 'IT2_reduced', 'IT4': 'IT4_reduced', 'IT5A': 'IT5A_full', 'IT5B': 'IT5B_reduced', 
                    'PT5B': 'PT_full', 'IT6': 'IT6_reduced', 'CT6': 'CT6_reduced', 'PV2': 'PV_simple', 'SOM2': 'SOM_simple'}):
    epsp = calculateEPSPsPops(params, data, somaLabel=somaLabel, stimRange=stimRange)


    pops = [p for p in params[2]['values']]
    secs = [s for s in params[0]['values']]
    locs = [s for s in params[1]['values']]
    segs=[]
    for sec,loc in zip(secs,locs): segs.append((sec,loc))
    segs.sort()
    weightNorm = {}

    for pop in pops:
        print(pop)
        weightNorm[pop] = {}
        for seg in segs: weightNorm[pop][seg[0]] = []  # empty list for each section
        for seg in segs:
            epspSeg = epsp[pop][tuple(seg)]
            epspSeg.sort()
            x,y = zip(*epspSeg)
            print(x,y)
            f = interp1d(y,x,fill_value="extrapolate")
            w = f(epspNorm)
            print(w)
            wnorm = w / epspNorm
            print(wnorm)
            weightNorm[pop][seg[0]].append(wnorm)
            print('\n%s %s wscale = %.6f' % (pop, str(seg), wnorm))

        if savePath:
            import pickle
            with open(savePath+popSaveLabels[pop]+'_weightNorm.pkl','wb') as fileObj:
                pickle.dump(weightNorm[pop], fileObj)


    return weightNorm



def plotEPSPs(epsp, dataFolder, batchLabel, addLegend=True, includeSegs = None):
    utils.setPlotFormat(numColors = 8)

    if len(params) == 3:
        pops = ['_']
        segs = includeSegs if includeSegs else [s for s in params[0]['values']]
        weights = params[2]['values']
        epspPops = {'_': epsp}
    elif len(params) == 4:
        pops = params[2]['values']
        segs = includeSegs if includeSegs else epsp[pops[0]].keys() # params[1]['values']
        weights = params[3]['values']
        epspPops = epsp

    for pop in pops:
        plt.figure(figsize=((12,8)))

        for seg in segs:
            if not seg[0].startswith('axon') and seg[1] not in [0.0,1.0]:
                epspSeg = epspPops[pop][tuple(seg)]
                if epspSeg:
                    epspSeg.sort()
                    x,y = zip(*epspSeg)
                    handles = plt.plot(x, y, marker='o', markersize=10, label=str(seg))

        plt.xlabel('Weight (of NetStim connection)')
        #xtick = np.arange(0.0, 0.0022, 0.0002)
        #plt.xticks(xtick, xtick)
        plt.ylabel('Somatic EPSP amplitude (mV) in response to 1 NetStim spike')
        if addLegend: plt.legend(title = 'Section', loc=2)
        if includeSegs:
            plt.savefig('%s/%s/%s_%s_epsp_subset.png' % (dataFolder, batchLabel, batchLabel, pop))
        else:
            plt.savefig('%s/%s/%s_%s_epsp.png' % (dataFolder, batchLabel, batchLabel, pop))
    
    #plt.show()



# main code
if __name__ == '__main__':

    # run batch E cells
    
    dataFolder = '../data/'
    batchLabels = ['v22_batch13'] #, 'v22_batch19'] #, 'v22_batch1', 'v22_batch16', 'v22_batch17']
    loadFromFile = 1

    ''' run via batch.py
    b = batch.weightNormE(pops=['IT2', 'IT4'], rule='IT2_reduced', weight=[0.0001])
    b.batchLabel = batchLabel  
    b.saveFolder = dataFolder+b.batchLabel
    b.method = 'grid'
    setRunCfg(b, 'mpi')
    b.run() # run batch
    '''

    popSaveLabels = {'IT2': 'IT2_reduced', 'IT3': 'IT3_reduced', 'ITP4': 'ITP4_reduced', 'ITS4': 'ITS4_reduced',
                     'IT5A': 'IT5A_reduced', 'CT5A': 'CT5A_reduced', 'IT5B': 'IT5B_reduced', 'CT5B': 'CT5B_reduced', 'PT5B': 'PT5B_reduced',
                     'IT6': 'IT6_reduced', 'CT6': 'CT6_reduced',
                     'PV2': 'PV_reduced', 'SOM2': 'SOM_reduced', 'VIP2': 'VIP_reduced', 'NGF2': 'NGF_reduced',
                     'IRE': 'RE_reduced', 'TC': 'TC_reduced', 'HTC': 'HTC_reduced', 'TI': 'TI_reduced'}
    

    # analyze batch E cells    
    for batchLabel in batchLabels:
        params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadFromFile, saveAll=1-loadFromFile, vars=[('simData','V_soma')], maxCombs=None) 
        epsp = calculateEPSPsPops(params, data, somaLabel='soma', stimRange=[10*700,10*800], syn='exc')
        plotEPSPs(epsp, dataFolder, batchLabel, addLegend=1)
        #plotEPSPs(epsp, dataFolder, batchLabel, addLegend=1, includeSegs=[('apic_28',0.5), ('apic_36',0.5), ('apic_49',0.5), ('apic_56',0.5)])
        weightNorm = calculateWeightNormPops(params, data,  somaLabel='soma', stimRange=[10*700,10*800], popSaveLabels=popSaveLabels, savePath=dataFolder+'/'+batchLabel+'/')