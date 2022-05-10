"""
wscale.py

Code to calculate the background input weights for each cell type

Contributors: salvadordura@gmail.com
"""

import utils
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from cycler import cycler
import sys, os
import json, pickle


def axisFontSize(ax, fontsize):
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 

def extractPopRates(params, data, popSubset=[]):
    if popSubset:
        Lvals = [pop for pop in params[0]['values'] if pop in pops]
    else:
        Lvals = params[0]['values']
    Ivals = params[1]['values']

    Lvalsdic = {val: i for i,val in enumerate(Lvals)}
    Ivalsdic = {val: i for i,val in enumerate(Ivals)}

    rates = [[0 for x in range(len(Ivals))] for y in range(len(Lvals))] 
    for key, d in data.items():
        try:
            rate = d['popRates'][d['paramValues'][0]]
            
            Lindex = Lvalsdic[d['paramValues'][0]]
            Iindex = Ivalsdic[d['paramValues'][1]]
            #print(Lindex, Iindex)
            rates[Lindex][Iindex] = rate
            #print(d['paramValues'])
            #print(rate)
        except:
            print('missing',Lindex,Iindex)
            pass

    filename = '%s/%s/%s_rateVsWeight.json' % (dataFolder, batchLabel, batchLabel)
    with open(filename, 'w') as fileObj:
        json.dump(rates, fileObj)
    
    return Lvals, Ivals, rates

def plotRateVsWeight(pops, weights, rates, showLegend = True, legendLabel ='', maxRate=50):
    fontsiz = 5
    plt.rcParams.update({'font.size': fontsiz})

    plt.figure(figsize=(25, 18))

    fig, ax  = plt.subplots(4, 5, sharey='row', sharex='col')

    for i, pop in enumerate(pops):
        x,y = np.unravel_index(i, (4,5))
        ax[x][y].plot(rates[i], marker='o', markersize=2)

        if i==0:
            ax[x][y].set_xlabel('Bkg weight (EPSP mV)')
            ax[x][y].set_ylabel('Rate (Hz)')
        ax[x][y].set_title(pop)
        ax[x][y].set_xlim(0,len(weights))
        ax[x][y].set_ylim(0,maxRate)
        axisFontSize(ax[x][y], fontsiz)
    
    plt.tight_layout()
    plt.savefig('%s/%s/%s_rateVsWeight.png' % (dataFolder, batchLabel, batchLabel), dpi=600)


def animateRateVsWeight(dataFolder, batchLabel, params):
    import imageio
    from pathlib import Path

    Lvals = params[0]['values']
    Ivals = params[1]['values']

    for ipop, pop in enumerate(Lvals):
        print('Generating traces gif for pop %s ...' % (pop))
        #v22_batch1_18_98_traces.png
        images = ['%s/%s/%s_%d_%d_traces.png' % (dataFolder, batchLabel, batchLabel, ipop, iweight) for iweight in range(len(Ivals))]
        #images = list(image_path.glob())
        image_list = []
        for file_name in images:
            image_list.append(imageio.imread(file_name))
        pass
        imageio.mimwrite('%s/%s/%s_%s_traces.gif' % (dataFolder, batchLabel, batchLabel, pop), image_list)


def calculateBkgWeightPops(pops, weights, rates, targetRates = {}, manualScaling = {}, savePath=None, saveToCells=False):
    bkgWeights = {}

    print ('\nCalculating bkg weights for each pop ...')
    for ipop, pop in enumerate(pops):
        success = False
        numPoints = 100 #20
        numPointsIncrease = 10
        while not success:
            x = weights[:numPoints]
            y = rates[ipop][:numPoints] 
            f = interp1d(y, x, fill_value="extrapolate")
            w = f(targetRates[pop])
            if (0.0 < w < 150) or numPoints == 110:
                success = True
            else:
                numPoints += numPointsIncrease

        if w < 0.0: w = 0.1
            

        bkgWeights[pop] = float(w)

        print('Bkg weight to get %.1f Hz in pop %s = %.1f' %(targetRates[pop], pop, w))

    if savePath:

        # fill in bkgWeights so all pops have the corresponding value for its cell type
        # note we only calculated bkg weights for each cell type, since several pops use same cell type
        # e.g. PV2, PV4, PV5A, PV5B, PV6 will all use the value calculated for PV2 

        allpops = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'PV4', 'SOM4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']

        import difflib
        for pop in allpops:
            if pop not in bkgWeights:
                popClose = difflib.get_close_matches(pop[:3], bkgWeights.keys())
                #print(pop,popClose)
                bkgWeights[pop] = bkgWeights[popClose[0]]
        
        for pop in allpops:
            if pop in manualScaling:
                bkgWeights[pop] *= manualScaling[pop]
        
        print ('\nFinal saved weights after completing and manual scaling: \n', bkgWeights)

        with open('%s/bkgWeightPops.json' % savePath,'w') as f:
            json.dump(bkgWeights, f)

        if saveToCells:
            with open('../cells/bkgWeightPops.json', 'w') as f:
                json.dump(bkgWeights, f)


    return bkgWeights


# main code
if __name__ == '__main__':

    # run batch E cells
    
    dataFolder = '../data/'
    batchLabel = 'v22_batch28'
    loadFromFile = 1

    ''' run via batch.py
    cellTypes = ['IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'HTC', 'IRE', 'TI']
    b = bkgWeights(pops = cellTypes, weights = list(range(1,100)))
    b.batchLabel = 'v22_batch21' 
    b.saveFolder = 'data/'+b.batchLabel
    b.method = 'grid'  # evol
    setRunCfg(b, 'mpi_bulletin') 
    b.run() # run batch
    '''

    # analyze batch E cells    
    params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadFromFile, saveAll=1 - loadFromFile, vars=[('simData', 'popRates')], maxCombs=None)
    pops, weights, rates = extractPopRates(params, data)
    plotRateVsWeight(pops, weights, rates)
    #animateRateVsWeight(dataFolder, batchLabel, params)

    Erate = 0.1  # Hz
    Irate = 0.1 # Hz
    targetRates = {p: Erate for p in ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'HTC']}
    targetRates.update({p: Irate for p in ['PV2', 'SOM2', 'VIP2', 'NGF2', 'IRE', 'TI', 'TIM']})

    # manual scaling adjustments based on simulation with automatically calculated bkg weights (finetuning)
    NGFfactor = 2.5
    
    ITS4factor = 25 # make similar to ITP4 since updated cell but did not update bkg weight estimates (v22_batch28)

    manualScaling = {'NGF1': 1.0 * NGFfactor, 'SOM2': 0.75, 'VIP2': 0.75, 'NGF2': 1.0 * NGFfactor, 
                    'SOM3': 1.0, 'VIP3': 1.25, 'NGF3': 1.0 * NGFfactor, 
                    'ITP4': 1.1, 'ITS4': ITS4factor, 'SOM4': 1.0, 'PV4': 0.9, 'VIP4': 1.0, 'NGF4': 1.0 * NGFfactor, 
                    'IT5A': 0.075, 'CT5A': 0.75, 'SOM5A': 1.25, 'PV5A': 1.25, 'VIP5A': 1.1, 'NGF5A': 1.0 * NGFfactor, 
                    'PT5B': 2.0, 'IT5B': 0.075, 'CT5B': 0.75, 'SOM5B': 1.25, 'PV5B': 1.25, 'VIP5B': 1.1, 'NGF5B': 1.0 * NGFfactor, 
                    'IT6': 0.1, 'CT6': 0.75, 'SOM6': 1.1, 'PV6': 0.75, 'NGF6': 1.0 * NGFfactor, 
                    'TC': 1.25, 'TCM': 1.25, 'HTC': 1.25, 'TI': 1.25, 'TIM': 2.0}  
    # run calculation
    bkgWeights = calculateBkgWeightPops(pops, weights, rates, targetRates, manualScaling,
                                        savePath=dataFolder + '/' + batchLabel + '/', saveToCells=True)
