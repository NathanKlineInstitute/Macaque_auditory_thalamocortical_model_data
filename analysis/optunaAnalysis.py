# analyzeEvol.py 

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv
import json
import seaborn as sns
import IPython as ipy
import os
import utils
from pprint import pprint
import optuna
import pickle

def getParamLabels(dataFolder, batchSim):
    # get param labels
    with open('%s/%s/%s_batch.json' % (dataFolder, batchSim, batchSim), 'r') as f: 
        paramLabels = [str(x['label'][0])+str(x['label'][1]) if isinstance(x['label'], list) else str(x['label']) for x in json.load(f)['batch']['params']]
    return paramLabels

def loadData(dataFolder, batchSim, pops, rateTimeRanges = [], loadStudyFromFile=False, loadDataFromFile=False, maxNumber=None):
 
    if loadStudyFromFile:
        with open('%s/%s/%s_df_study.pkl' % (dataFolder, batchSim, batchSim), 'rb') as f:
            df = pickle.load(f)
    else:
        study = optuna.create_study(study_name=batchSim, storage='sqlite:///%s/%s/%s_storage.db' % (dataFolder, batchSim, batchSim), load_if_exists=True) 
        df = study.trials_dataframe(attrs=('number', 'value', 'params'))

        # replace column labels
        for col in df.columns:
            if col.startswith('params_'):
                newName = col.replace('params_', '').replace('(', '').replace(')', '').replace("'", "").replace(', ','')
                df = df.rename(columns={col: newName})

        with open('%s/%s/%s_df_study.pkl' % (dataFolder, batchSim, batchSim), 'wb') as f:
            pickle.dump(df, f)

    if loadDataFromFile:
        try:
            with open('%s/%s/%s_df.pkl' % (dataFolder, batchSim, batchSim), 'rb') as f:
                df = pickle.load(f)
        except:
            print('Could not find _df.pkl file...')
    else:
        # load json for each trial with pop rates and add to df
        popRates = {p: [] for p in pops}
        for p in list(popRates.keys()):
            for t in rateTimeRanges:
                try:
                    popRates[p+'_'+t] = []
                except:
                    print('%s not found ...' % (p+'_'+t))

        if not maxNumber:
            maxNumber = len(df.number)
        for i in range(maxNumber):
            # try:
                filename = '%s/%s/trial_%d/trial_%d' % (dataFolder, batchSim, int(i), int(i))
                if os.path.exists(filename+'.json'):
                    with open(filename+'.json', 'r') as f:
                        popRatesLoad = json.load(f)['simData']['popRates']
                elif os.path.exists(filename + '.pkl'):
                    with open(filename+'.pkl', 'rb') as f:
                        popRatesLoad = pickle.load(f)['simData']['popRates']          
                for p in popRatesLoad:
                    popRates[p].append(np.mean(list(popRatesLoad[p].values())))
                    for t in popRatesLoad[p].keys():
                        popRates[p+'_'+t].append(popRatesLoad[p][t])

                print('Added trial %d' % (i))
            # except:
            #     for p in popRates:
            #         popRates[p].append(0.0)
            #     print('Skipped trial %d' % (i))
            
        for p, rates in popRates.items():
            df.insert(len(df.columns), p, rates+[0]*(len(df)-len(rates)))

        with open('%s/%s/%s_df.pkl' % (dataFolder, batchSim, batchSim), 'wb') as f:
            pickle.dump(df, f)


    return df

def plotScatterPopVsParams(dataFolder, batchsim, df, pops, skipCols=[], ylim=None):

    dfcorr=df.corr('pearson')
    
    for pop in pops:

        for param in df.columns:
            if not any([skipCol in param for skipCol in skipCols]): 
                try:
                    print('Plotting scatter of %s vs %s param (R=%.2f) ...' %(pop, param, dfcorr[pop][param]))
                    df.plot.scatter(param, pop, s=4, c='number', colormap='viridis', alpha=0.5, figsize=(8, 8), colorbar=False)
                    if ylim: plt.ylim(0, ylim)
                    plt.title('%s vs %s R=%.2f' % (pop, param, dfcorr[pop][param]))
                    plt.savefig('%s/%s/%s_scatter_%s_%s.png' %(dataFolder, batchSim, batchSim, pop, param), dpi=300)
                except:
                    print('Error plotting %s vs %s' % (pop,param))


def plotScatterFitnessVsParams(dataFolder, batchsim, df, excludeAbove=None, skipCols=[]):

    if excludeAbove:
        df = df[df.value < excludeAbove]

    dfcorr=df.corr('pearson')

    for param in df.columns:
        if not any([skipCol in param for skipCol in skipCols]): 
            try:
                print('Plotting scatter of %s vs %s param (R=%.2f) ...' %('fitness', param, dfcorr['value'][param]))
                df.plot.scatter(param, 'value', s=4, c='number', colormap='viridis', alpha=0.5, figsize=(8, 8), colorbar=False)
                plt.ylabel('fitness error')
                plt.title('%s vs %s R=%.2f' % ('fitness', param.replace('tune', ''), dfcorr['value'][param]))
                plt.savefig('%s/%s/%s_scatterEvol_%s_%s.png' % (dataFolder, batchSim, batchSim, 'fitness', param.replace('tune', '')), dpi=300)
            except:
                print('Error plotting %s vs %s' % ('fitness',param))


def plotScatterTrialVsParams(dataFolder, batchsim, df, excludeAbove=None, skipCols=[]):

    if excludeAbove:
        df = df[df.value < excludeAbove]

    dfcorr=df.corr('pearson')

    for param in df.columns:
        if not any([skipCol in param for skipCol in skipCols]): 
            try:
                print('Plotting scatter of %s vs %s param (R=%.2f) ...' %('fitness', param, dfcorr['value'][param]))
                #df.plot.scatter(param, 'value', s=4, c='number', colormap='viridis', alpha=0.5, figsize=(8, 8), colorbar=False)
                df.plot.scatter('number', param, s=4, c='value', colormap='viridis', alpha=0.5, figsize=(8, 8), colorbar=True, vmin=100, vmax=400)
                plt.ylabel('param value')
                plt.xlabel('trial')
                plt.title('%s vs %s R=%.2f' % ('trial', param.replace('tune', ''), dfcorr['value'][param]))
                plt.savefig('%s/%s/%s_scatter_%s_%s.png' % (dataFolder, batchSim, batchSim, 'trial', param.replace('tune', '')), dpi=300)
            except:
                print('Error plotting %s vs %s' % ('fitness',param))

                
def plotJointplotFitnessVsParams(dataFolder, batchsim, df, excludeAbove=None):

    if excludeAbove:
        df = df[df.value < excludeAbove]

    dfcorr=df.corr('pearson')

    for param in df.columns:
        
        try:
            g = sns.jointplot(x=list(df[param]), y=list(df['value']), scatter=False, kind='kde', fill=True, thresh=0, cmap='Blues', height=8)
            g.set_axis_labels(param, 'fitness', fontsize=16)
            #plt.title('%s vs %s R=%.2f' % ('fitness', param.replace('tune', ''), dfcorr['value'][param]))
            plt.subplots_adjust(left=0.1, bottom=0.1)#,right=0.95, top=0.9, bottom=0.1)
            plt.savefig('%s/%s/%s_jointplot_%s_%s.png' % (dataFolder, batchSim, batchSim, 'fitness', param.replace('tune', '')), dpi=300)
        
        except:
             print('Error plotting %s vs %s' % ('fitness',param))


def plotParamsVsFitness(dataFolder, batchSim, df, paramLabels, excludeAbove=None, ylim=None):

    if excludeAbove:
        df = df[df.value < excludeAbove]

    df2 = df.drop(['value', 'number'], axis=1)
    fits = list(df['value'])
    plt.figure(figsize=(16, 12))

    for i, (k,v) in enumerate(list(df2[paramLabels].items())):
        y = v #(np.array(v)-min(v))/(max(v)-min(v)) # normalize
        x = np.random.normal(i, 0.04, size=len(y))         # Add some random "jitter" to the x-axis
        s = plt.scatter(x, y, alpha=0.3, c=[int(f-1) for f in fits], cmap='jet_r')
    plt.colorbar(label = 'fitness')
    plt.ylabel('Parameter value')
    plt.xlabel('Parameter')
    plt.subplots_adjust(top=0.95, bottom=0.2, right=0.95)
    if ylim: plt.ylim(0, ylim)
    plt.savefig('%s/%s/%s_scatter_params_%s.png' % (dataFolder, batchSim, batchSim, 'excludeAbove-'+str(excludeAbove) if excludeAbove else ''))
    #plt.show()



def filterRates(df, condlist=['rates', 'I>E', 'E5>E6>E2', 'PV>SOM'], Epops=[], Ipops=[], rateTimeRanges=[], 
                copyFolder=None, dataFolder=None, batchLabel=None, skipDepol=False, verbose=True):
    from os.path import isfile, join
    from glob import glob

    # allpops = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'PV4', 'SOM4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI']

    ranges = {}
    Erange = [0.01,200]
    #Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B',  'IT6', 'CT6', 'TC',  'HTC', 'TCM']
    L23 = ['IT2', 'IT3' ]
    L4 = [ 'ITP4', 'ITS4']
    L5A = ['IT5A', 'CT5A']
    L5B = ['IT5B', 'PT5B', 'CT5B']
    L6 = ['IT6', 'CT6']
    thal = ['TC',  'HTC', 'TCM']

    #Epops = L23+L4

    for pop in Epops:
        ranges[pop] = Erange
    
    # check pop rate ranges
    if rateTimeRanges:
        dfcond = df
        if 'rates' in condlist:
            for t in rateTimeRanges:
                conds = []
                for k,v in ranges.items(): conds.append(str(v[0]) + '<=' + k+'_'+t + '<=' + str(v[1]))
                condStr = ''.join([''.join(str(cond) + ' and ') for cond in conds])[:-4]
                if condStr: dfcond = dfcond.query(condStr)
    else:
        conds = []
        if 'rates' in condlist:
            for k,v in ranges.items(): conds.append(str(v[0]) + '<=' + k + '<=' + str(v[1]))
        condStr = ''.join([''.join(str(cond) + ' and ') for cond in conds])[:-4]
        dfcond = df.query(condStr)

    if verbose:
        print('\n Filtering based on E pops: ' + str(condlist) + '\n' + condStr)
        print(dfcond)
        print(len(dfcond))

    ranges = {}
    Irange = [0.01,200]
    L1 = ['NGF1']
    L2 = ['PV2', 'SOM2', 'VIP2']#, 'NGF2']      # L2
    L3 = ['PV3', 'SOM3', 'VIP3']#, 'NGF3']      # L3
    L4 = ['PV4', 'SOM4', 'VIP4']#, 'NGF4']      # L4
    L5A = ['PV5A', 'SOM5A', 'VIP5A']#, 'NGF5A']  # L5A  
    L5B = ['PV5B', 'SOM5B', 'VIP5B']#, 'NGF5B'] #,  # L5B
    L6 = ['PV6', 'SOM6', 'VIP6']#, 'NGF6'] #,
    thal = ['IRE', 'IREM', 'TI', 'TIM']      # L6 PV6
 
    #Ipops = L1+L2+L3

    for pop in Ipops:
        ranges[pop] = Irange

    # check pop rate ranges

    if rateTimeRanges:
        if 'rates' in condlist:
            for t in rateTimeRanges:
                conds = []
                for k,v in ranges.items(): conds.append(str(v[0]) + '<=' + k+'_'+t + '<=' + str(v[1]))
                condStr = ''.join([''.join(str(cond) + ' and ') for cond in conds])[:-4]
                if condStr: dfcond = dfcond.query(condStr)
    else:
        conds = []
        if 'rates' in condlist:
            for k,v in ranges.items(): conds.append(str(v[0]) + '<=' + k + '<=' + str(v[1]))
        condStr = ''.join([''.join(str(cond) + ' and ') for cond in conds])[:-4]
        dfcond = dfcond.query(condStr)


    # # check I > E in each layer
    # if 'I>E' in condlist:
    #     conds.append('PV2 > IT2 and SOM2 > IT2')
    #     conds.append('PV5A > IT5A and SOM5A > IT5A')
    #     conds.append('PV5B > IT5B and SOM5B > IT5B')
    #     conds.append('PV6 > IT6 and SOM6 > IT6')

    # # check E L5 > L6 > L2
    # if 'E5>E6>E2' in condlist:
    #     #conds.append('(IT5A+IT5B+PT5B)/3 > (IT6+CT6)/2 > IT2')
    #     conds.append('(IT5A+IT5B+PT5B)/3 > (IT6+CT6)/2')
    #     conds.append('(IT6+CT6)/2 > IT2')
    #     conds.append('(IT5A+IT5B+PT5B)/3 > IT2')
    
    # # check PV > SOM in each layer
    # if 'PV>SOM' in condlist:
    #     conds.append('PV2 > IT2')
    #     conds.append('PV5A > SOM5A')
    #     conds.append('PV5B > SOM5B')
    #     conds.append('PV6 > SOM6')


    # construct query and apply
    # condStr = ''.join([''.join(str(cond) + ' and ') for cond in conds])[:-4]
    # dfcond = df.query(condStr)

    if verbose:
        print('\n Filtering based on E + I pops: ' + str(condlist) + '\n' + condStr)
        print(dfcond)
        print(len(dfcond))

    # copy files
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
                sourceFile1 = dataFolder+batchLabel+'/gen_'+i.split('_')[0]+'/gen_'+i.split('_')[0]+'_cand_'+i.split('_')[1]+'_*raster*.png'   
            #sourceFile2 = dataFolder+batchLabel+'/'+batchLabel+row['simLabel']+'.json'
            if len(glob(sourceFile1))>0:
                cpcmd = 'cp ' + sourceFile1 + ' ' + targetFolder + '/.'
                #cpcmd = cpcmd + '; cp ' + sourceFile2 + ' ' + targetFolder + '/.'
                os.system(cpcmd) 
                print(cpcmd)


    return dfcond


def calculateParamImportance(dataFolder, batchSim):
    study = optuna.create_study(study_name=batchSim, storage='sqlite:///%s/%s/%s_storage.db' % (dataFolder, batchSim, batchSim), load_if_exists=True) 
    importance = optuna.importance.get_param_importances(study=study)
    print(importance)
    return importance




# -----------------------------------------------------------------------------
# Main code
# -----------------------------------------------------------------------------
if __name__ == '__main__': 
    dataFolder = '../data/'
    batchSim = 'v34_batch31'
    loadFromFile = 1
    
    allpops = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'PV4', 'SOM4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']  #, 'IC']
    #allpops = ['IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'PV4', 'SOM4', 'VIP4', 'NGF4','IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM'] # , 'IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']  #, 'IC']

    #rateTimeRanges = ['1000_1250', '1250_1500', '1500_1750', '1750_2000']
    rateTimeRanges = ['1500_1750', '1750_2000', '2000_2250', '2250_2500']

    #pd.options.display.max_columns = None  
    #pd.options.display.max_rows = None  

    # set font size
    plt.rcParams.update({'font.size': 18})

    # get param labelsc
    paramLabels = getParamLabels(dataFolder, batchSim)

    # load evol data from files
    df = loadData(dataFolder, batchSim, pops=allpops, rateTimeRanges=rateTimeRanges, loadStudyFromFile=loadFromFile, loadDataFromFile=loadFromFile)#, maxNumber=3038)

    #plotParamsVsFitness(dataFolder, batchSim, df, paramLabels, excludeAbove=200, ylim=None)

    #plotScatterFitnessVsParams(dataFolder, batchSim, df, excludeAbove=None, skipCols=rateTimeRanges)

    #plotScatterTrialVsParams(dataFolder, batchSim, df, excludeAbove=None, skipCols=rateTimeRanges)

    #plotJointplotFitnessVsParams(dataFolder, batchSim, df, excludeAbove=500)

    #plotScatterPopVsParams(dataFolder, batchSim, df, pops = ['NGF6','SOM6', 'SOM2', 'SOM3'], skipCols=rateTimeRanges)
