"""
batchAnalysis.py 

Code to anlayse batch sim results

Contributors: salvadordura@gmail.com
"""

import utils
import json, pickle
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sb
import os
from batchAnalysisPlotCombined import dfPopRates
import IPython as ipy
plt.style.use('seaborn-whitegrid')

# ---------------------------------------------------------------------------------
# Support funcs
# ---------------------------------------------------------------------------------

def filterRates(df, condlist=['rates', 'I>E', 'E5>E6>E2', 'PV>SOM'], copyFolder=None, dataFolder=None, batchLabel=None, skipDepol=True):
    from os.path import isfile, join
    from glob import glob

    rangesE = {}
    rangesI = {}

    Erange = [0.01,200]
    Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B', 'PT5B', 'IT6','CT6'] # , 'IT5A', 'CT5A', 'IT5B', 'CT5B', 'PT5B', 'IT6', 'CT6']  # all layers
#'ITS4',
    for pop in Epops:
        rangesE[pop] = Erange

    Irange = [0.01,200]
    #Ipops = ['PV2']

    Ipops = ['NGF1',                        # L1
        'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
        'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
        'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
        'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
        'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
        'PV6', 'SOM6', 'VIP6', 'NGF6']      # L6 

    for pop in Ipops:
        rangesI[pop] = Irange

    conds = []

    # check pop rate ranges
    if 'rates' in condlist:
        for k,v in rangesE.items(): conds.append(str(v[0]) + '<=' + k + '<=' + str(v[1]))
    
    # check I > E in each layer
    if 'I>E' in condlist:
        conds.append('PV2 > IT2 and SOM2 > IT2')
        conds.append('PV5A > IT5A and SOM5A > IT5A')
        conds.append('PV5B > IT5B and SOM5B > IT5B')
        conds.append('PV6 > IT6 and SOM6 > IT6')

    # check E L5 > L6 > L2
    if 'E5>E6>E2' in condlist:
        #conds.append('(IT5A+IT5B+PT5B)/3 > (IT6+CT6)/2 > IT2')
        conds.append('(IT5A+IT5B+PT5B)/3 > (IT6+CT6)/2')
        conds.append('(IT6+CT6)/2 > IT2')
        conds.append('(IT5A+IT5B+PT5B)/3 > IT2')
    
    # check PV > SOM in each layer
    if 'PV>SOM' in condlist:
        conds.append('PV2 > IT2')
        conds.append('PV5A > SOM5A')
        conds.append('PV5B > SOM5B')
        conds.append('PV6 > SOM6')

    # construct query and apply
    condStr = ''.join([''.join(str(cond) + ' and ') for cond in conds])[:-4]
    

    #import IPython; IPython.embed()

    dfcond = df.query(condStr)

    conds = []
    if 'rates' in condlist:
        for k, v in rangesI.items(): conds.append(str(v[0]) + '<=' + k + '<=' + str(v[1]))
    condStr = ''.join([''.join(str(cond) + ' and ') for cond in conds])[:-4]

    dfcond = dfcond.query(condStr)
            
    print('\n Filtering based on: ' + str(condlist) + '\n' + condStr)
    print(dfcond)
    print(len(dfcond))



    # 
    if copyFolder:
        targetFolder = dataFolder+batchLabel+'/'+copyFolder
        try: 
            os.mkdir(targetFolder)
        except:
            pass
        
        for i,row in dfcond.iterrows():     
            if skipDepol:
                sourceFile1 = dataFolder+batchLabel+'/noDepol/'+batchLabel+row['simLabel']+'*.png'  
            else:
                sourceFile1 = dataFolder+batchLabel+'/'+batchLabel+row['simLabel']+'*.png'   
            #sourceFile2 = dataFolder+batchLabel+'/'+batchLabel+row['simLabel']+'.json'
            if len(glob(sourceFile1))>0:
                cpcmd = 'cp ' + sourceFile1 + ' ' + targetFolder + '/.'
                #cpcmd = cpcmd + '; cp ' + sourceFile2 + ' ' + targetFolder + '/.'
                os.system(cpcmd) 
                print(cpcmd)


    return dfcond

def detectDepolBlock(volt, minSamples=200, vRange=[-50, -10]):
    cumm = 0
    for v in volt:
        if v >= vRange[0] and v<= vRange[1]:
            cumm += 1
            if cumm > minSamples:
                return 1
        else:
            cumm = 0

    return 0 


# ---------------------------------------------------------------------------------
# Wrapper funcs
# ---------------------------------------------------------------------------------
               

def applyFilterRates(dataFolder, batchLabel, loadAll, skipDepol=True):
    # load data 
    var = [('simData','popRates')]
    params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, vars=var, maxCombs=None)
 
    #pipy.embed()

    # convert to pandas and add pop Rates
    df1 = utils.toPandas(params, data)
    dfpop = dfPopRates(df1)

    # filter based on pop rates
    condLists = [['rates']]#,
                # ['rates','I>E'],
                # ['rates', 'E5>E6>E2'],
                # ['rates', 'PV>SOM']]

                # [['I>E', 'E5>E6>E2', 'PV>SOM'],
                # ['I>E'],
                # ['E5>E6>E2'],
                # ['PV>SOM'],
                # ['rates']]
                # ['E5>E6>E2'],
                # ['PV>SOM']]


    for i,condList in enumerate(condLists):
        print(len(filterRates(dfpop, condlist=condList, copyFolder='selected_'+str(i), dataFolder=dataFolder, batchLabel=batchLabel, skipDepol=skipDepol)))
    
    return dfpop


def filterStimRates(dataFolder, batchLabel, simLabel=None, load=False, subset=''):
    from os import listdir
    from os.path import isfile, join

    # from mpi4py import MPI
    # comm = MPI.COMM_WORLD
    # rank = comm.Get_rank()
    # mode = MPI.MODE_CREATE | MPI.MODE_WRONLY 

    ihlevels = [0,1,2]
    ihmax = 3

    for ihlevel in ihlevels:
        jsonFolder = batchLabel 
        path = dataFolder+batchLabel #+'/noDepol'
        #onlyFiles = [f for f in listdir(path) if isfile(join(path, f)) and not f.endswith('batch.json') and not f.endswith('cfg.json')]
        onlyFiles = [f.replace('_raster.png', '.json') for f in listdir(path) if isfile(join(path, f)) and f.endswith('raster.png')]

        if type(simLabel) is list:
            outfiles = [f for f in onlyFiles if f in simLabel] #any([f.endswith(sl+'.json') for sl in simLabel])] 
        elif type(simLabel) is '':
            outfiles = [f for f in onlyFiles if f.endswith(simLabel+'.json')]
        else:
            outfiles = [f for f in onlyFiles if f.endswith('_%d.json'%(ihlevel))] 

        include = {}
        include['rates'] = ['IT5A', 'PT5B', 'IT2']

        targetFolder = dataFolder+batchLabel+'/stimRates'
        try: 
            os.mkdir(targetFolder)
        except:
            pass

        if load:
            filename = dataFolder+jsonFolder+'/'+'stimRatesData'+subset+'%d.json'%(ihlevel)
            with open(filename, 'r') as fileObj: loadData = json.load(fileObj)

        with open('../sim/cells/popColors.pkl', 'r') as fileObj: popColors = pickle.load(fileObj)['popColors']
        saveData = []
        counter=0
        copied = 0
        for ifile, outfile in enumerate(outfiles):
            if 1: #rank == ifile:    
                try: 
                    filename = dataFolder+jsonFolder+'/'+outfile
                    #print filename     

                    if load:
                        filename = filename.replace('_%d.json'%(ihlevel), '_%d.json'%(ihmax))
                        found=0
                        for d in loadData:
                            if str(d[0]) == str(filename):
                                filename, avgsPT5Bzd, peaksPT5Bzd, avgsPT5Bih, peaksPT5Bih, avgsIT5Azd, peaksIT5Azd, avgsIT5Aih, peaksIT5Aih, \
                                    avgsIT2zd, peaksIT2zd, avgsIT2ih, peaksIT2ih = d
                                found=1
                                counter=counter+1
                        if not found: continue

                    else:
                        sim,data,out=utils.plotsFromFile(filename, raster=0, stats=0, rates=1, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0, 
                            timeRange=[500,1500], include=include, textTop='', popColors=popColors, orderBy='gid')

                        avgsIT5Azd, avgsPT5Bzd, avgsIT2zd = [out[2][0][0],out[2][1][0]], [out[2][0][1],out[2][1][1]], [out[2][0][2],out[2][1][2]]
                        peaksIT5Azd, peaksPT5Bzd, peaksIT2zd = [out[3][0][0],out[3][1][0]], [out[3][0][1],out[3][1][1]], [out[3][0][2],out[3][1][2]]

                        # ih = 1.0 (ends with 1)
                        filename = dataFolder+jsonFolder+'/'+outfile.replace('_%d.json'%(ihlevel), '_%d.json'%(ihmax))
                        print(filename)     

                        sim,data,out=utils.plotsFromFile(filename, raster=0, stats=0, rates=1, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0, 
                            timeRange=[500,1500], include=include, textTop='', popColors=popColors, orderBy='gid')

                        avgsIT5Aih, avgsPT5Bih, avgsIT2ih = [out[2][0][0],out[2][1][0]], [out[2][0][1],out[2][1][1]], [out[2][0][2],out[2][1][2]]
                        peaksIT5Aih, peaksPT5Bih, peaksIT2ih = [out[3][0][0],out[3][1][0]], [out[3][0][1],out[3][1][1]], [out[3][0][2],out[3][1][2]]

                        saveData.append([filename, avgsPT5Bzd, peaksPT5Bzd, avgsPT5Bih, peaksPT5Bih, avgsIT5Azd, peaksIT5Azd, avgsIT5Aih, peaksIT5Aih, 
                                        avgsIT2zd, peaksIT2zd, avgsIT2ih, peaksIT2ih])

                    # conditions to select 
                    avgsPT5BzdInc = (avgsPT5Bzd[1] - avgsPT5Bzd[0]) / avgsPT5Bzd[0] if avgsPT5Bzd[0]>0 else 0 # % inc
                    avgsPT5BihInc = (avgsPT5Bih[1] - avgsPT5Bih[0]) / avgsPT5Bih[0] if avgsPT5Bih[0]>0 else 0# % inc

                    peaksPT5BzdInc = (peaksPT5Bzd[1] - peaksPT5Bzd[0]) / peaksPT5Bzd[0]  if peaksPT5Bzd[0]>0 else 0# % inc
                    peaksPT5BihInc = (peaksPT5Bih[1] - peaksPT5Bih[0]) / peaksPT5Bih[0]  if peaksPT5Bih[0]>0 else 0 #% inc


                    # peaks PT5B
                    if (peaksPT5Bzd[1] > 1.0*peaksPT5Bih[1] and
                    # peaksPT5BzdInc > 1.0*peaksPT5BihInc and
                    peaksPT5BzdInc > 0 and
                    
                    # peaks IT2
                    # peaksIT2zd[0] < 250.0 and
                    # peaksIT2ih[0] < 250 and
                    # peaksIT2zd[1] < 250.0 and
                    # peaksIT2ih[1] < 250  and
                    abs(peaksPT5BihInc) < 1.5 and 

                    #avgs PT5B
                    avgsPT5Bzd[1] > 1.0*avgsPT5Bih[1] and  # 1.5
                    # avgsPT5BzdInc > 1.0*avgsPT5BihInc and  # 1.25
                    avgsPT5BzdInc > 0 and 
                    avgsPT5Bih[0] > 0. and
                    # abs(avgsPT5BihInc) < 1.5 and
                    
                    # avgs IT2
                    avgsIT2zd[0] > 0.05 and
                    avgsIT2ih[0] > 0.05 and
                    avgsIT2zd[1] > 0.05 and
                    avgsIT2ih[1] > 0.05): #and

                        copied += 1
                        # print avgsPT5Bzd, avgsPT5BzdInc
                        # print peaksPT5Bzd, peaksPT5BzdInc

                        # print avgsPT5Bih, avgsPT5BihInc
                        # print peaksPT5Bih, peaksPT5BihInc

                        sourceFiles = []
                        sourceFiles.append(filename)
                        sourceFiles.append(filename.split('.json')[0]+'_raster.png')
                        sourceFiles.append(filename.split('.json')[0]+'_spikeHist.png')
                        sourceFiles.append(filename.split('.json')[0]+'_avgRates.png')
                        sourceFiles.append(filename.split('.json')[0]+'_peakRates.png')
                        filename2 = filename.replace('_%d.json'%(ihmax), '_%d.json'%(ihlevel))
                        sourceFiles.append(filename2)
                        sourceFiles.append(filename2.split('.json')[0]+'_raster.png')
                        sourceFiles.append(filename2.split('.json')[0]+'_spikeHist.png')
                        sourceFiles.append(filename2.split('.json')[0]+'_avgRates.png')
                        sourceFiles.append(filename2.split('.json')[0]+'_peakRates.png')

                        for sf in sourceFiles:
                            cpcmd = 'cp ' + sf + ' ' + targetFolder + '/.'
                            try:
                                os.system(cpcmd) 
                                #print cpcmd
                            except:
                                pass
                    else: 
                        pass
                        #print 'not copied'

                    if not load:
                        filename = dataFolder+jsonFolder+'/'+'stimRatesData'+subset+'%d.json'%(ihlevel)
                        with open(filename, 'w') as fileObj: json.dump(saveData, fileObj)

                        # filename2 = dataFolder+jsonFolder+'/'+'stimRatesData_mpi'+subset+'.json'
                        # fh = MPI.File.Open(comm, filename2, mode) 
                        # line1 = str(saveData) + '\n' 
                        # fh.Write_ordered(line1) 
                        # fh.Close() 

                except:
                   pass

        print('Total: ', counter)
        print('Copied: ', copied)


def filterDepolBlock(dataFolder, batchLabel, loadAll, gids=None): 
    if gids:  #  e.g. gids = [3136, 5198])
        var = [('simData')]
        params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, vars=var, maxCombs=None)
        df = utils.toPandas(params, data)

        copyFolder = 'noDepol'
        targetFolder = dataFolder+batchLabel+'/'+copyFolder
        try: 
            os.mkdir(targetFolder)
        except:
            pass
        
        for i,row in df.iterrows():
            block = False
            for gid in gids:
                if detectDepolBlock(row.V_soma['cell_'+str(gid)]): block = True

            if not block:
                sourceFile1 = dataFolder+batchLabel+'/'+batchLabel+row['simLabel']+'_*.png'
                #sourceFile2 = dataFolder+batchLabel+'/'+batchLabel+row['simLabel']+'.json'
                cpcmd = 'cp ' + sourceFile1 + ' ' + targetFolder + '/.'
                #cpcmd = cpcmd + '; cp ' + sourceFile2 + ' ' + targetFolder + '/.'
                os.system(cpcmd) 
                print(cpcmd)

        return df

    else:
        var = [('simData', 'V_soma', 'cell_5198'), ('simData', 'simLabel')]
        params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, vars=var, maxCombs=None)
        df = utils.toPandas(params, data)

        copyFolder = 'noDepol'
        targetFolder = dataFolder+batchLabel+'/'+copyFolder
        try: 
            os.mkdir(targetFolder)
        except:
            pass
        
        for i,row in df.iterrows():
            block = False
            if detectDepolBlock(row.V_soma): block = True

            if not block:
                sourceFile1 = dataFolder+batchLabel+'/'+batchLabel+row.simLabel+'_*.png'
                sourceFile2 = dataFolder+batchLabel+'/'+batchLabel+row['simLabel']+'.json'
                cpcmd = 'cp ' + sourceFile1 + ' ' + targetFolder + '/.'
                cpcmd = cpcmd + '; cp ' + sourceFile2 + ' ' + targetFolder + '/.'
                os.system(cpcmd) 
                print(cpcmd)
                
    return df



def addFitness(dataFolder, batchLabel, loadAll, tranges=[[1500, 1750], [1750,2000], [2000,2250], [2250,2500]], Epops=None, Ipops=None):
        
    params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll)
    df = utils.toPandas(params, data)

    # calculate fitness value for each row
    fitnessFuncArgs = {}
    pops = {}

    ## Exc pops
    if not Epops:
        Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'TCM', 'HTC']  # all layers + thal + IC

    Etune = {'target': 5, 'width': 5, 'min': 0.5}
    
    for pop in Epops:
        pops[pop] = Etune
    
    ## Inh pops 
    if not Ipops:
        Ipops = ['NGF1',                            # L1
                'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
                'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
                'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
                'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
                'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
                'PV6', 'SOM6', 'VIP6', 'NGF6',       # L6
                'IRE', 'IREM', 'TI']  # Thal 


    Itune = {'target': 10, 'width': 15, 'min': 0.5}

    for pop in Ipops:
        pops[pop] = Itune
    
    fitnessFuncArgs['pops'] = pops
    fitnessFuncArgs['maxFitness'] = 1000
    fitnessFuncArgs['tranges'] = tranges

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
        return fitness

    df['fitness'] = [fitnessFunc(simData=row, **fitnessFuncArgs) for index, row in df.iterrows()]

    # add also avg rates 
    for k in df.iloc[0]['popRates'].keys():
        df[k] = [np.mean([float(x) for x in r['popRates'][k].values()]) for i,r in df.iterrows()]

    return df