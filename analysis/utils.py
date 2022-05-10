"""
utils.py 

General functions to analyse simulation data

Contributors: salvadordura@gmail.com
"""
import json
import pickle
import numpy as np
from pylab import *
from itertools import product
import pandas as pd
from pprint import pprint
import os

from netpyne import specs
from collections import OrderedDict

def toPandas(params, data):
    if 'simData' in data[list(data.keys())[0]]:
        rows = [list(d['paramValues'])+[s for s in d['simData'].values()] for d in data.values()]
        cols = [str(d['label']) for d in params]+[s for s in data[list(data.keys())[0]]['simData'].keys()]
    else:
        rows = [list(d['paramValues'])+[s for s in d.values()] for d in data.values()]
        cols = [str(d['label']) for d in params]+[s for s in data[list(data.keys())[0]].keys()]
    
    df = pd.DataFrame(rows, columns=cols) 
    df['simLabel'] = data.keys()

    colRename=[]
    for col in list(df.columns):
        if col.startswith("[u'"):
            colName = col.replace(", u'","_'").replace("[u","").replace("'","").replace("]","").replace(", ","_")
            colRename.append(colName)
        else: 
            colRename.append(col)
    print(colRename)
    df.columns = colRename

    return df


def setPlotFormat(numColors=8):
    # plt.style.use('ggplot')
    #plt.style.use(['dark_background', 'presentation'])
    plt.style.use('seaborn-whitegrid')

    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 'large'

    NUM_COLORS = numColors
    colormap = plt.get_cmap('nipy_spectral')
    colorlist = [colormap(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
    
    # add red and blue at the beginning
    # colorlist.insert(0, [1, 0, 0])
    # colorlist.insert(0, [0, 0, 1])

    plt.rc('axes', prop_cycle=(cycler('color', colorlist)))


def compare(source_file, target_file, source_key=None, target_key=None):
    from deepdiff import DeepDiff 
    with open(source_file, 'r') as fileObj:
        if source_file.endswith('.json'):
            source = json.load(fileObj, object_pairs_hook=OrderedDict)
        elif source_file.endswith('.pkl'):
            source = pickle.load(fileObj)
    if source_key: source = source[source_key]

    with open(target_file, 'r') as fileObj:
        if target_file.endswith('.json'):
            target = json.load(fileObj, object_pairs_hook=OrderedDict)
        elif source_file.endswith('.pkl'):
            target = pickle.load(fileObj)
    if target_key: target = target[target_key]
    
    ddiff = DeepDiff(source, target)
    pprint(ddiff)
    return ddiff


def readBatchData(dataFolder, batchLabel, loadAll=False, saveAll=True, vars=None, maxCombs=None, listCombs=None):
    # load from previously saved file with all data
    if loadAll:
        print('\nLoading single file with all data...')
        filename = '%s/%s/%s_allData.json' % (dataFolder, batchLabel, batchLabel)
        with open(filename, 'r') as fileObj:
            dataLoad = json.load(fileObj, object_pairs_hook=OrderedDict)
        params = dataLoad['params']
        data = dataLoad['data']
        return params, data

    if isinstance(listCombs, str):
        filename = str(listCombs)
        with open(filename, 'r') as fileObj:
            dataLoad = json.load(fileObj)
        listCombs = dataLoad['paramsMatch']

    # read the batch file and cfg
    batchFile = '%s/%s/%s_batch.json' % (dataFolder, batchLabel, batchLabel)
    with open(batchFile, 'r') as fileObj:
        b = json.load(fileObj)['batch']

    # read params labels and ranges
    params = b['params']

    # reorder so grouped params come first
    preorder = [p for p in params if 'group' in p and p['group']]
    for p in params:
        if p not in preorder: preorder.append(p)
    params = preorder

    # read vars from all files - store in dict 
    if b['method'] == 'grid':
        labelList, valuesList = zip(*[(p['label'], p['values']) for p in params])
        valueCombinations = product(*(valuesList))
        indexCombinations = product(*[range(len(x)) for x in valuesList])
        data = {}
        print('Reading data...')
        missing = 0
        for i,(iComb, pComb) in enumerate(zip(indexCombinations, valueCombinations)):
            if (not maxCombs or i<= maxCombs) and (not listCombs or list(pComb) in listCombs):
                print(i, iComb)
                # read output file
                iCombStr = ''.join([''.join('_'+str(i)) for i in iComb])
                simLabel = b['batchLabel']+iCombStr
                outFile = dataFolder+'/'+batchLabel+'/'+simLabel
                if os.path.isfile(outFile+'.json'):
                    outFile = outFile + '.json'
                    with open(outFile, 'rb') as fileObj:
                        output = json.load(fileObj, object_pairs_hook=OrderedDict)
                elif os.path.isfile(outFile+'.pkl'):
                    outFile = outFile + '.pkl'
                    with open(outFile, 'rb') as fileObj:
                        output = pickle.load(fileObj)
                else:
                    print('... file missing')
                    missing = missing + 1
                    output = {}
                    continue

                try:
                    # save output file in data dict
                    data[iCombStr] = {}  
                    data[iCombStr]['paramValues'] = pComb  # store param values
                    if not vars: vars = output.keys()

                    for key in vars:
                        if isinstance(key, tuple):
                            container = output
                            for ikey in range(len(key)-1):
                                container = container[key[ikey]]
                            data[iCombStr][key[1]] = container[key[-1]]

                        elif isinstance(key, str): 
                            data[iCombStr][key] = output[key]
                except:
                    print('... file missing')
                    missing = missing + 1
                    output = {}

                    #import IPython; IPython.embed()

            else:
                missing = missing + 1

        print('%d files missing' % (missing))

        # save
        if saveAll:
            print('Saving to single file with all data')
            filename = '%s/%s/%s_allData.json' % (dataFolder, batchLabel, batchLabel)
            dataSave = {'params': params, 'data': data}
            with open(filename, 'w') as fileObj:
                json.dump(dataSave, fileObj)
        
        return params, data

def loadFromFile(filename, pops=None):
    from netpyne import specs,sim
    import collections

    with open(filename, 'rb') as fileObj:
        data = json.load(fileObj, object_pairs_hook=OrderedDict)

    cfg = specs.SimConfig(data['simConfig'])
    cfg.createNEURONObj = False

    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data, instantiate=False)
    sim.setSimCfg(cfg)

    if pops:    
        temp = collections.OrderedDict()

        # order pops by keys in popColors
        for k in pops:
            if k in sim.net.params.popParams:
                temp[k] = sim.net.params.popParams.pop(k)

        # add remaining pops at the end
        for k in list(sim.net.params.popParams.keys()):
            temp[k] = sim.net.params.popParams.pop(k)
        
        sim.net.params.popParams = temp
    
    try:
        print('Cells created: ',len(sim.net.allCells))
    except:
        sim.net.createPops()     
        sim.net.createCells()
        sim.setupRecording()
        sim.gatherData() 

    sim.allSimData = data['simData']

    return sim



def plotsFromFile(filename, raster=True, stats=True, rates=False, syncs=False, hist=True, psd=True, grang=True, traces=True, plotAll=False, 
                timeRange=None, textTop='', include={'raster':[]}, popColors=None, orderBy='gid'):
    from netpyne import specs
    from collections import OrderedDict

    # try:
    with open(filename, 'rb') as fileObj:
        if filename.endswith('.json'):
            data = json.load(fileObj, object_pairs_hook=OrderedDict)
        elif filename.endswith('.pkl'):
            data = pickle.load(fileObj)
    # except:
    #     print('Error opening file %s' % (filename))
    #     return 0
    sim,data, out = plotsFromData(data, raster=raster, stats=stats, rates=rates, syncs=syncs, hist=hist, psd=psd, traces=traces, grang=grang,
                             plotAll=plotAll, timeRange=timeRange, textTop=textTop, include=include, popColors=popColors, orderBy=orderBy)
    return sim,data, out


def plotsFromData(data, textTop = '', raster=True, stats=False, rates=False, syncs=False, hist=True, psd=True, traces=True, grang=True, 
                plotAll=0, timeRange=None, include=None, popColors=None, orderBy='gid'):
    import matplotlib.pyplot as plt
    from netpyne import specs,sim
    cfg = specs.SimConfig(data['simConfig'])
    cfg.createNEURONObj = False

    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data, instantiate=False)
    sim.setSimCfg(cfg)
    
    import collections
    temp = collections.OrderedDict()
    # order pops by keys in popColors
    for k in include['raster']:
        if k in sim.net.params.popParams:
            temp[k] = sim.net.params.popParams.pop(k)

    # add remaining pops at the end
    for k in list(sim.net.params.popParams.keys()):
        temp[k] = sim.net.params.popParams.pop(k)
    
    sim.net.params.popParams = temp
    
    try:
        print('Cells created: ',len(sim.net.allCells))
    except:
        #import IPython; IPython.embed()
        sim.net.createPops()     
        sim.net.createCells()
        sim.setupRecording()
        sim.gatherData() 

    sim.allSimData = data['simData']

    plt.style.use('seaborn-ticks') # seaborn-whitegrid 'seaborn-deep', 'seaborn-whitegrid', 'seaborn-poster'

    out = None

    if raster or plotAll:
        fig1 = sim.analysis.plotRaster(include=include['raster'], timeRange=timeRange, popRates=True, orderInverse=True, lw=0, markerSize=10, 
            marker='.', popColors=popColors, showFig=0, saveFig=0, figSize=(10,11), orderBy=orderBy)
        ax = plt.gca()
        plt.text(0.05, 1.07, textTop, transform=ax.transAxes, fontsize=12, color='k')
        plt.ylabel('Neuron number')
        plt.title('')
        plt.savefig('../'+cfg.filename+'_raster_%d_%d_%s.png'%(timeRange[0], timeRange[1], orderBy), dpi=600)

    if stats or plotAll:
        filename = cfg.filename+'_%d_%d'%(timeRange[0], timeRange[1])
        fig1 = sim.analysis.plotSpikeStats(include=include['stats'], figSize=(4,8), timeRange=timeRange, stats = ['rate', 'isicv'], popColors=popColors, showFig=0, saveFig=filename)
        #fig1 = sim.analysis.plotSpikeStats(include=include['stats'], xlim=[0,1], timeRange=timeRange, stats = ['isicv'], popColors=popColors, showFig=0, saveFig=filename)
        #, 'pairsync', 'sync'] # xlim[]

    if rates or plotAll:
        midpoint = (timeRange[1]+timeRange[0]) / 2.0
        timeRanges = [[timeRange[0], midpoint], [midpoint, timeRange[1]]]
        colors = [popColors[inc] if inc in popColors else (0,0,0) for inc in include['rates']]
        out = sim.analysis.plotRates(include=include['rates'], timeRanges=timeRanges, figSize=(4,2), timeRangeLabels=['Pre stim', 'Post stim'], colors=colors, showFig=0, saveFig=1)

    if syncs or plotAll:
        midpoint = (timeRange[1]+timeRange[0]) / 2.0
        timeRanges = [[timeRange[0], midpoint], [midpoint, timeRange[1]]]
        colors = [popColors[inc] if inc in popColors else (0,0,0) for inc in include['syncs']]
        fig1 = sim.analysis.plotSyncs(include=include['syncs'], timeRanges=timeRanges, figSize=(5,4), timeRangeLabels=['Pre stim', 'Post stim'], colors=colors, showFig=0, saveFig=1)

    if hist or plotAll:
        fig2 = sim.analysis.plotSpikeHist(include=include['hist'], yaxis='rate', binSize=5, graphType='bar', timeRange=timeRange,  popColors=popColors, showFig=False, saveFig=0, figSize=(12,8))
        ax = plt.gca()
        plt.text(0.05, 1.07, textTop, transform=ax.transAxes, fontsize=12, color='k')
        plt.savefig(cfg.filename+'_spikeHist_replot_%d_%d.png'%(timeRange[0], timeRange[1]))

    if psd or plotAll:
        popColors['allCells'] = 'k'
        fig3 = sim.analysis.plotRatePSD(include=include['psd'], timeRange=timeRange, Fs=160, smooth=16 , showFig=0, saveFig=0, popColors=popColors, figSize=(11,5))
        # ylim=[-55, 0]
        ax = plt.gca()
        plt.text(0.05, 1.07, textTop, transform=ax.transAxes, fontsize=12, color='k')
        plt.savefig(cfg.filename+'_spikePSD_%d_%d.png'%(timeRange[0], timeRange[1]), dpi=300)

    if traces or plotAll:
        colors = popColors; colors=['r','b']
        fig4 = sim.analysis.plotTraces(include=include['traces'], timeRange=timeRange, overlay = True, oneFigPer = 'trace', rerun = False, ylim=[-90, 30], axis='on', figSize = (12,4), saveData = None, saveFig = '../'+cfg.filename+'_traces.png', showFig = 0)


    if grang or plotAll:
        print('Calculating and plotting Granger...')
        sim.allSimData = data['simData']
        spkts = sim.allSimData['spkt']
        spkids = sim.allSimData['spkid']
        spkids,spkts = zip(*[(spkid,spkt-timeRange[0]) for spkid,spkt in zip(spkids,spkts) if timeRange[0] <= spkt <= timeRange[1]])

        popLabels = sim.net.allPops.keys()
        gidPops = [cell['tags']['pop'] for cell in sim.net.allCells]
        popNumCells = [float(gidPops.count(pop)) for pop in popLabels]
        spktPop = {}
        for pop, popNum in zip(popLabels, popNumCells):
            spktPop[pop] = [spkt for spkt,spkid in zip(spkts,spkids) if sim.net.allCells[int(spkid)]['tags']['pop']==pop]


        F, Fx2y, Fy2x, Fxy = granger(spktPop['IT2'], spktPop['PT5B'])
        F_b, Fx2y_b, Fy2x_b, Fxy_b = granger(spktPop['IT2'], spktPop['IT5A'])
        F_b2, Fx2y_b2, Fy2x_b2, Fxy_b2 = granger(spktPop['IT5A'], spktPop['PT5B'])
        # F_c, Fx2y_c, Fy2x_c, Fxy_c = granger(spktPop['S2'], spktPop['IT5A'])
        # F_d, Fx2y_d, Fy2x_d, Fxy_d = granger(spktPop['M2'], spktPop['PT5B'])
        # F_e, Fx2y_e, Fy2x_e, Fxy_e = granger(spktPop['S2'], spktPop['PT5B'])
        # F_f, Fx2y_f, Fy2x_f, Fxy_f = granger(spktPop['M2'], spktPop['IT5A'])
        

        figh = figure(figsize=(11,5))
        plot(F_b, Fy2x_b, 'r-', label = 'IT2 -> IT5A')
        plot(F_b, Fx2y_b, 'r:', label = 'IT5A -> IT2')
        plot(F, Fy2x, 'b-', label = 'IT2 -> PT5B')
        plot(F, Fx2y, 'b:', label = 'PT5B -> IT2')
        plot(F, Fy2x_b2, 'm-', label = 'IT5A -> PT5B')
        plot(F, Fx2y_b2, 'm:', label = 'PT5B -> IT5A')
        # plot(F_c, Fy2x_c, 'r', label = 'S2 -> IT5A')
        # plot(F_e, Fy2x_e, 'c', label = 'S2 -> PT5B')
        # plot(F_f, Fy2x_f, 'tab:orange', label = 'M2 -> IT5A')
        # plot(F_d, Fy2x_d, 'b', label = 'M2 -> PT5B')
        
        #plot(F_c, Fx2y_c, 'r', label = 'IT5A -> S2')

        xlabel('Frequency (Hz)')
        ylabel('Granger Causality')
        legend()
        plt.ylim(0,5)
        plt.xlim(0,80)
        ax = plt.gca()
        plt.tight_layout()

        #plt.text(0.05, 1.07, textTop, transform=ax.transAxes, fontsize=12, color='k')
        plt.savefig(cfg.filename+'_Granger_%d_%d.png'%(timeRange[0], timeRange[1]), dpi=400)

    return sim, data, out


def getSpksPop(data, timeRange=None):
    from netpyne import specs,sim
    cfg = specs.SimConfig(data['simConfig'])
    cfg.createNEURONObj = False

    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data)
    sim.setSimCfg(cfg)
    sim.net.createPops()     
    sim.net.createCells()
    sim.gatherData() 

    sim.allSimData = data['simData']
    spkts = sim.allSimData['spkt']
    spkids = sim.allSimData['spkid']

    popLabels = sim.net.allPops.keys()
    gidPops = [cell['tags']['pop'] for cell in sim.net.allCells]
    popNumCells = [float(gidPops.count(pop)) for pop in popLabels]
    numCellsPop = {pop:float(gidPops.count(pop)) for pop in popLabels}
    spktPop = {}
    for pop, popNum in zip(popLabels, popNumCells):
        if timeRange:
            spktPop[pop] = [spkt for spkt,spkid in zip(spkts,spkids) if sim.net.allCells[int(spkid)]['tags']['pop']==pop if timeRange[0]<=spkt<=timeRange[1]]
        else:
            spktPop[pop] = [spkt for spkt,spkid in zip(spkts,spkids) if sim.net.allCells[int(spkid)]['tags']['pop']==pop]

    return spktPop, numCellsPop


def readBatchUpdatePopRates(dataFolder, batchLabel, timeRange=None, saveAll=True, loadAll=False):
    # load from previously saved file with all data
    if loadAll:
        from netpyne import specs
        print('\nLoading single file with all data...')
        filename = '%s/%s/%s_allData_popRates.json' % (dataFolder, batchLabel, batchLabel)
        with open(filename, 'r') as fileObj:
            dataLoad = json.load(fileObj, object_pairs_hook=OrderedDict)
        params = dataLoad['params']
        data = dataLoad['data']
        return params, data

    # load single sim to get netParams and simConfig
    params, data = readBatchData(dataFolder, batchLabel, maxCombs=1)  

    # recreate net to get cell+pop params
    from netpyne import specs,sim
    cfg = specs.SimConfig(data.values()[0]['simConfig'])
    cfg.createNEURONObj = False
    sim.initialize()  # create network object and set cfg and net params
    sim.loadAll('', data=data.values()[0])
    sim.setSimCfg(cfg)
    sim.net.createPops()     
    sim.net.createCells()
    sim.gatherData() 

    # caclulate pop stats
    if not timeRange: timeRange = [0, cfg.duration]
    tsecs = (timeRange[1]-timeRange[0])/1e3 
    popLabels = sim.net.allPops.keys()
    gidPops = [cell['tags']['pop'] for cell in sim.net.allCells]
    popNumCells = [float(gidPops.count(pop)) for pop in popLabels]
    #numCellsPop = {pop:float(gidPops.count(pop)) for pop in popLabels}
    popLabels = sim.net.allPops.keys()

    # load spkt and spkid
    params, data = readBatchData(dataFolder, batchLabel, vars=[('simData','spkt'),('simData','spkid'),('simData','popRates')], saveAll=False)
    for d in data.values():
        spkts = d['spkt']
        spkids = d['spkid']

        spktPop = {}
        for pop, popNum in zip(popLabels, popNumCells):
            spktPop[pop] = [spkt for spkt,spkid in zip(spkts,spkids) if sim.net.allCells[int(spkid)]['tags']['pop']==pop if timeRange[0]<=spkt<=timeRange[1]]
            d['popRates'][pop] = len(spktPop[pop]) / popNum / tsecs  # calcualte new popRate
            d['spkt'] = None  # remove to save space
            d['spkid'] = None  

    # save
    if saveAll:
        print('Saving to single file with all data')
        filename = '%s/%s/%s_allData_popRates.json' % (dataFolder, batchLabel, batchLabel)
        dataSave = {'params': params, 'data': data}
        with open(filename, 'w') as fileObj:
            json.dump(dataSave, fileObj)

    return params, data



def nTE(spk1, spk2, trange = [0,2000]):
    from neuron import h
    h.load_file('/usr/site/nrniv/local/hoc/nqs.hoc')
    h.load_file('/usr/site/nrniv/local/hoc/infot.hoc') # for nTE code
    inputVec = h.Vector()
    outputVec = h.Vector()
    binSize=20
    histo1 = histogram(spk1, bins = np.arange(trange[0], trange[1], binSize))
    histoCount1 = histo1[0] 
    histo2 = histogram(spk2, bins = np.arange(trange[0], trange[1], binSize))
    histoCount2 = histo2[0] 

    inputVec.from_python(histoCount1)
    outputVec.from_python(histoCount2)
    #te = inputVec.tentropspks(outputVec, vtmp, nshuf)
    nshuf = 30
    out = h.normte(inputVec, outputVec, nshuf)
    TE, H, nTE, _, _ = out.to_python()
    return nTE

def granger(spk1, spk2, binSize=5, trange=[0,2000]):
    """
    Typical usage is as follows:
    from bsmart import pwcausalr
    F,pp,cohe,Fx2y,Fy2x,Fxy=pwcausalr(x,ntrls,npts,p,fs,freq);

    Outputs:
        F is the frequency vector for the remaining quantities
    pp is the spectral power
    cohe is the coherence
    Fx2y is the causality of channel X to channel Y
    Fy2x is the causality of channel Y to channel X
    Fxy is the "instantaneous" causality (cohe-Fx2y-Fy2x I think)
    Inputs:
    x is the data for at least two channels, e.g. a 2x8000 array consisting of two LFP time series
    ntrls is the number of trials (whatever that means -- just leave it at 1)
    npts is the number of points in the data (in this example, 8000)
    p is the order of the polynomial fit (e.g. 10 for a smooth fit, 20 for a less smooth fit)
    fs is the sampling rate (e.g. 200 Hz)
    freq is the maximum frequency to calculate (e.g. fs/2=100, which will return 0:100 Hz)
    """
    
    from pylab import histogram, plot, show
    from netpyne.support.bsmart import pwcausalr

    histo1 = histogram(spk1, bins = np.arange(trange[0], trange[1], binSize))
    histoCount1 = histo1[0] 

    histo2 = histogram(spk2, bins = np.arange(trange[0], trange[1], binSize))
    histoCount2 = histo2[0] 

    fs = 1000/binSize
    F,pp,cohe,Fx2y,Fy2x,Fxy = pwcausalr(np.array([histoCount1, histoCount2]), 1, len(histoCount1), 10, fs, fs/2)

    return F, Fx2y[0],Fy2x[0], Fxy[0]



def stars(p):
    if p < 0.0001:
        return "**** (p<0.0001)"
    elif (p < 0.001):
        return "*** (p<0.001)"
    elif (p < 0.01):
        return "** (p<0.01)"
    elif (p < 0.05):
        return "* (p<0.1)"
    else:
        return "-"
            
def boxplot(data, labels, saveFile=None):
    import glob
    import scipy.stats

    import brewer2mpl
    bmap = brewer2mpl.get_map('Set2', 'qualitative', 7)
    colors = bmap.mpl_colors

    params = {
        'axes.labelsize': 14,
        'text.fontsize': 14,
        'legend.fontsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'text.usetex': False,
        'figure.figsize': [5, 9]
    }
    rcParams.update(params)

    # def load(dir):
    #     f_list = glob.glob(dir + '/*/*/bestfit.dat')
    #     num_lines = sum(1 for line in open(f_list[0]))
    #     i = 0;
    #     data = np.zeros((len(f_list), num_lines)) 
    #     for f in f_list:
    #         data[i, :] = np.loadtxt(f)[:,1]
    #         i += 1
    #     return data


    # data_low_mut = load('data/low_mut')
    # data_high_mut = load('data/high_mut')
    # low_mut_100 = data_low_mut[:, 100]
    # high_mut_100 =  data_high_mut[:, 100]

    fig = figure()
    ax = fig.add_subplot(111)

    bp = ax.boxplot(data, notch=0, sym='b+', vert=1, whis=1.5, 
                 positions=None, widths=0.6)


    # for i in range(len(bp['boxes'])):
    #     box = bp['boxes'][i]
    #     box.set_linewidth(0)
    #     boxX = []
    #     boxY = []
    #     for j in range(5):
    #         boxX.append(box.get_xdata()[j])
    #         boxY.append(box.get_ydata()[j])
    #         boxCoords = zip(boxX,boxY)
    #         boxPolygon = Polygon(boxCoords, facecolor = colors[i % len(colors)], linewidth=0)
    #         ax.add_patch(boxPolygon)

    for i in range(0, len(bp['boxes'])):
        bp['boxes'][i].set_color(colors[i])
        bp['boxes'][i].set_linewidth(2)
        # we have two whiskers!
        bp['whiskers'][i*2].set_color(colors[i])
        bp['whiskers'][i*2 + 1].set_color(colors[i])
        bp['whiskers'][i*2].set_linewidth(0)
        bp['whiskers'][i*2 + 1].set_linewidth(0)
        # top and bottom fliers
        # if len(bp['fliers']) > i*2:
           #  bp['fliers'][i * 2].set(markerfacecolor=colors[i],
           #                  marker='o', alpha=0.75, markersize=6,
           #                  markeredgecolor='none')
           #  bp['fliers'][i * 2 + 1].set(markerfacecolor=colors[i],
           #                  marker='o', alpha=0.75, markersize=6,
           #                  markeredgecolor='none')
        bp['medians'][i].set_color(colors[i])
        bp['medians'][i].set_linewidth(4)
        # and 4 caps to remove
        for c in bp['caps']:
            c.set_linewidth(0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', direction='out')
    ax.tick_params(axis='y', length=0)

    # scatter data
    for i,d in enumerate(data):
        scatter([i+1]*len(d), d, color=colors[i], s=20, marker='o', linewidth=1, facecolors='none')

    ax.grid(axis='y', color="0.9", linestyle='-', linewidth=1)
    ax.set_axisbelow(True)

    ax.set_xticklabels(labels)

    # the stars
    if len(data)==2:
        z, p = scipy.stats.mannwhitneyu(data[0], data[1])
        p_value = p * 2
        s = stars(p)

        y_max = np.max(np.concatenate((data[0], data[1])))
        y_min = np.min(np.concatenate((data[0], data[1])))
        ax.annotate("", xy=(1, y_max), xycoords='data',
                    xytext=(2, y_max), textcoords='data',
                    arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                    connectionstyle="bar,fraction=0.2"))
        ax.text(1.5, y_max + abs(y_max - y_min)*0.1, stars(p_value),
                horizontalalignment='center',
                verticalalignment='center') 


        fig.subplots_adjust(left=0.2)

    if len(data)==3:

        # 0-1
        z, p = scipy.stats.mannwhitneyu(data[0], data[1])
        p_value = p * 2
        s = stars(p_value)

        if s != '-':
            y_max = np.max(np.concatenate((data[0], data[1])))
            y_min = np.min(np.concatenate((data[0], data[1])))
            ax.annotate("", xy=(1, y_max), xycoords='data',
                        xytext=(2, y_max), textcoords='data',
                        arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                        connectionstyle="bar,fraction=0.2"))
            ax.text(1.5, y_max + abs(y_max - y_min)*0.11, stars(p_value),
                    horizontalalignment='center',
                    verticalalignment='center') 

        # 0-2
        z, p = scipy.stats.mannwhitneyu(data[0], data[2])
        p_value = p * 2
        s = stars(p)

        if s != '-':
            y_max = np.max(np.concatenate((data[0], data[2])))
            y_min = np.min(np.concatenate((data[0], data[2])))
            ax.annotate("", xy=(1, y_max), xycoords='data',
                        xytext=(3, y_max), textcoords='data',
                        arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                        connectionstyle="bar,fraction=0.2"))
            ax.text(2, y_max + abs(y_max - y_min)*0.11, stars(p_value),
                    horizontalalignment='center',
                    verticalalignment='center') 

        # 1-2
        z, p = scipy.stats.mannwhitneyu(data[1], data[2])
        p_value = p * 2
        s = stars(p)

        if s != '-':
            y_max = np.max(np.concatenate((data[1], data[2])))
            y_min = np.min(np.concatenate((data[1], data[2])))
            ax.annotate("", xy=(2, y_max), xycoords='data',
                        xytext=(3, y_max), textcoords='data',
                        arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
                                        connectionstyle="bar,fraction=0.2"))
            ax.text(2.5, y_max + abs(y_max - y_min)*0.06, stars(p_value),
                    horizontalalignment='center',
                    verticalalignment='center') 



        fig.subplots_adjust(left=0.2)


    if saveFile:
        savefig(saveFile, dpi=600)


def move_element(odict, thekey, newpos):
    odict[thekey] = odict.pop(thekey)
    i = 0
    for key, value in odict.items():
        if key != thekey and i >= newpos:
            odict[key] = odict.pop(key)
        i += 1
    return odict


def addPopRatesToJson(dataFolder, batchLabel):
    from os import listdir
    from os.path import isfile, join
    from netpyne import sim
    
    path = dataFolder+batchLabel 
    outFiles = [f for f in listdir(path) if isfile(join(path, f)) and f.endswith('.pkl')]

    for outfile in outFiles:
        try:
            filepath = path+outfile
            sim.load(filepath, instantiate=False)
            sim.allSimData['popRates']=sim.analysis.popAvgRates(tranges=[[1500, 1750], [1750,2000], [2000,2250], [2250,2500]])
            os.rename(filepath, filepath[:-4]+'_orig.pkl')
            sim.saveData(filename=filepath[:-4])
        except:
            print('Error in file %s' % (outfile))

        