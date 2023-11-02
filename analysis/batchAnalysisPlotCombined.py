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
plt.style.use('seaborn-whitegrid')

# ---------------------------------------------------------------------------------
# Support funcs
# ---------------------------------------------------------------------------------

def axisFontSize(ax, fontsize):
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 


def plotfI(dataFolder, batchLabel, params, data, saveLabel='', maxRate=150, pops=None):
    utils.setPlotFormat(numColors = 8)
    
    if pops:
        Lvals = [pop for pop in params[0]['values'] if pop in pops]
    else:
        Lvals = params[0]['values']
    Ivals = params[1]['values']

    Lvalsdic = {val: i for i,val in enumerate(Lvals)}
    Ivalsdic = {val: i for i,val in enumerate(Ivals)}

    rates = [[0 for x in range(len(Ivals))] for y in range(len(Lvals))] 
    for key, d in data.items():
        rate = len(d['simData']['spkt'])
        try:
            Lindex = Lvalsdic[d['paramValues'][0]]
            Iindex = Ivalsdic[d['paramValues'][1]]

            rates[Lindex][Iindex] = rate
            print(d['paramValues'])
            print(rate)
        except:
            pass

    filename = '%s/%s/%s_fI%s.json' % (dataFolder, batchLabel, batchLabel, saveLabel)
    with open(filename, 'w') as fileObj:
        json.dump(rates, fileObj)

    fontsiz = 5
    plt.rcParams.update({'font.size': fontsiz})

    # PLotting format options
    plt.figure(figsize=(25, 18))
    nrows = 4
    ncols = 5
    fig, ax=plt.subplots(nrows, ncols, sharey='row', sharex='col')
    maxRate = 150

    #import IPython; IPython.embed()

    for i, pop in enumerate(Lvals):
        x, y = np.unravel_index(i, (nrows, ncols))
        ax[x][y].plot(rates[i], marker='o', markersize=2)
        ax[x][y].set_xticks(list(range(len(Ivals)))[::2])
        ax[x][y].set_xticklabels(Ivals[::2])

        if i==0:
            ax[x][y].set_xlabel('Current amplitude (nA)', fontsize=fontsiz)
            ax[x][y].set_ylabel('Rate (Hz)', fontsize=fontsiz)
        ax[x][y].set_title(pop, fontsize=fontsiz*2)
        ax[x][y].set_ylim(0,maxRate)
        axisFontSize(ax[x][y], fontsiz)

    #if showLegend: plt.legend(handles, Lvals, title = legendLabel, loc=2)
    plt.tight_layout()
    plt.savefig('%s/%s/%s_fI%s.png' % (dataFolder, batchLabel, batchLabel, saveLabel), dpi=600)
    #plt.show()


def plot1Dparams(df, par1, val, valLabel=None, varLabel=None, groupStat='mean', normRange=False, saveFile=None, line=False):
    if isinstance(val, list):
        if not varLabel: varLabel = 'var'
        if not valLabel: valLabel = 'value'
        dfmelt = pd.melt(df, id_vars=[par1], value_vars=val, var_name=varLabel, value_name=valLabel)
        print(dfmelt)
        plot2Dparams(dfmelt, par1, varLabel, valLabel, groupStat=groupStat, cluster=False, normRange=normRange, line=True)
        return

    if not valLabel: valLabel = val
    dfsubset = df[[par1,val]] 
    dfgroup = dfsubset.groupby(by=[par1])
    if groupStat=='first':
        dfgroup2 = dfgroup.first()
    elif groupStat=='last':
        dfgroup2 = dfgroup.last()
    elif groupStat=='mean':
        dfgroup2 = dfgroup.mean()
    elif groupStat=='sum':
        dfgroup2 = dfgroup.sum()
    dfgroup3 = pd.DataFrame(dfgroup2).reset_index()
    print(dfgroup3)
    dfpiv = dfgroup3.pivot(index=par1, values=val)
    #dfpiv=dfgroup3
    utils.setPlotFormat(numColors = len(dfpiv.columns))
    dfpiv.plot(marker='o')
    
    normStr = '_norm' if normRange else ''
    if saveFile:
        plt.savefig(saveFile)
    else:
        plt.savefig(dataFolder+batchLabel+'/heatmap_'+par1+'_'+val+normStr+'.png')
    plt.show()
    return dfpiv


def plot2Dparams(df, par1, par2, val, valLabel=None, varLabel=None, groupStat='mean', cluster=False, line=False, normRange=False, saveFile=None):
    if isinstance(val, list):
        if len(val) == 1: 
            val = val[0]
        else:
            if not varLabel: varLabel = 'var'
            if not valLabel: valLabel = 'value'
            dfmelt = pd.melt(df, id_vars=[par1,par2], value_vars=val, var_name=varLabel, value_name=valLabel)
            plot3Dparams(dfmelt, par1, par2, varLabel, valLabel, groupStat=groupStat, cluster=False, normRange=normRange)
            return
    elif val == par1 or val == par2: # output value depends on param - eg. pop
        df['rate'] = [df.iloc[i][df.iloc[i][val]] for i in range(len(df))] # create combined column with rate of each pop
        val = 'rate'

    if isinstance(par2, list):
        dfsubset = df[[par1]+par2+[val]] 
        pd.options.mode.chained_assignment = None
        groupLabel = ''.join([', '+str(par) for par in par2])[1:]
        dfsubset[groupLabel] = [''.join([', '+str(dfsubset.iloc[i][par]) for par in par2])[1:] for i in range(len(dfsubset))]
        dffinal = dfsubset[[par1,groupLabel,val]] 
        par2 = groupLabel
    else:
        if not valLabel: valLabel = val
        dfsubset = df[[par1,par2,val]] 
        dfgroup = dfsubset.groupby(by=[par1,par2])
        if groupStat=='first':
            dfgroup2 = dfgroup.first()
        elif groupStat=='last':
            dfgroup2 = dfgroup.last()
        elif groupStat=='mean':
            dfgroup2 = dfgroup.mean()
        elif groupStat=='sum':
            dfgroup2 = dfgroup.sum()
        dffinal = pd.DataFrame(dfgroup2).reset_index()

    dfpiv = pd.pivot_table(dffinal, index=par1, columns=par2, values=val)
#    pandas.pivot_table(df,values='count',index='site_id',columns='week')
    if cluster:
        sb.clustermap(dfpiv, square=True, cbar_kws={'label': valLabel})
    elif line:
        utils.setPlotFormat(numColors = len(dfpiv.columns))

        #dfpiv = dfpiv[['IT2','IT4','IT5A','IT5B','PT5B','IT6','CT6']]
        dfpiv.plot(marker='o')  
    else:
        sb.heatmap(dfpiv, square=True, cbar_kws={'label': valLabel})
    normStr = '_norm' if normRange else ''
    try:
        if saveFile:
            plt.savefig(saveFile)
        else:
            plt.savefig(dataFolder+batchLabel+'/heatmap_'+par1+'_'+par2+'_'+val+normStr+'.png')
    except:
        print('Error saving figure...')

    plt.show()
    return dfpiv


def plot3Dparams(df, par1, par2, par3, val, valLabel=None, groupStat='mean', cluster=False, normRange=False, saveFile=None):
    if isinstance(val, list):
        valGroup = list(val)
        for ival in valGroup:
            plot3Dparams(df, par1, par2, par3, ival, ival+' '+valLabel, groupStat, cluster, normRange)
        return

    if not valLabel: valLabel = val
    dfsubset = df[[par1,par2,par3,val]] 
    dfgroup = dfsubset.groupby(by=[par1,par2,par3])
    if groupStat=='first':
        dfgroup2 = dfgroup.first()
    elif groupStat=='last':
        dfgroup2 = dfgroup.last()
    elif groupStat=='mean':
        dfgroup2 = dfgroup.mean()
    elif groupStat=='sum':
        dfgroup2 = dfgroup.sum()
    dfgroup3 = pd.DataFrame(dfgroup2).reset_index()

    def draw_heatmap(*args, **kwargs):
        data = kwargs.pop('data')
        d = data.pivot(index=par2, columns=par1, values=val)
        if cluster:
            sb.clustermap(d, **kwargs)
        else:
            sb.heatmap(d, **kwargs)

    fg = sb.FacetGrid(dfgroup3, col=par3)
    vmin,vmax = (dfgroup3[val].min(), dfgroup3[val].max()) if normRange else (None,None)
    fg.map_dataframe(draw_heatmap, par1, par2, val, cbar_kws={'label': valLabel}, square=True, vmin=vmin, vmax=vmax)#, cmap='jet')

    # get figure background color
    facecolor=plt.gcf().get_facecolor()
    for ax in fg.axes.flat:
        # set aspect of all axis
        #ax.set_aspect('equal','box-forced')
        # set background color of axis instance
        ax.set_axis_bgcolor(facecolor)
    normStr = '_norm' if normRange else ''
    if saveFile:
        plt.savefig(saveFile)
    else:
        plt.savefig(dataFolder+batchLabel+'/heatmap_'+par1+'_'+par2+'_'+val+normStr+'.png')    
    plt.show()
    return dfpiv


def plot4Dparams(dataFolder, batchLabel, df, par1, par2, par3, par4, val, valLabel=None, groupStat='mean', cluster=False, normRange=False, saveFile=None):
    if isinstance(val, list):
        valGroup = list(val)
        for ival in valGroup:
            plot4Dparams(dataFolder, batchLabel, df, par1, par2, par3, par4, ival, ival+' '+valLabel, groupStat, cluster, normRange)
        return

    if not valLabel: valLabel = val

    dfsubset = df[[par1, par2, par3, par4, val]]    
    dfgroup = dfsubset.groupby(by=[par1,par2,par3,par4])
    if groupStat=='first':
        dfgroup2 = dfgroup.first()
    elif groupStat=='last':
        dfgroup2 = dfgroup.last()
    elif groupStat=='mean':
        dfgroup2 = dfgroup.mean()
    elif groupStat=='sum':
        dfgroup2 = dfgroup.sum()
    dfgroup3 = pd.DataFrame(dfgroup2).reset_index()

    def draw_heatmap(*args, **kwargs):
        data = kwargs.pop('data')
        d = data.pivot(index=par2, columns=par1, values=val)
        if cluster:
            sb.clustermap(d, **kwargs)
        else:
            sb.heatmap(d, **kwargs)

    fg = sb.FacetGrid(dfgroup3, col=par3, row=par4)
    vmin,vmax = (dfgroup3[val].min(), dfgroup3[val].max()) if normRange else (None,None)
    fg.map_dataframe(draw_heatmap, par1, par2, val, cbar_kws={'label': valLabel}, square=True, vmin=vmin, vmax=vmax, cmap='viridis')

    # get figure background color
    facecolor=plt.gcf().get_facecolor()
    for ax in fg.axes.flat:
        # set aspect of all axis
        #ax.set_aspect('equal','box-forced')
        # set background color of axis instance
        ax.set_facecolor(facecolor)
    normStr = '_norm' if normRange else ''

    if saveFile:
        plt.savefig(saveFile)
    else:
        plt.savefig(dataFolder+batchLabel+'/heatmap_'+par1+'_'+par2+'_'+par3+'_'+par4+'_'+val+'_'+groupStat+normStr+'.png')
    #plt.show()
    return dfgroup3


def dfPopRates(df, numParams=6, pops=[]):
    dfpop = df.iloc[:,0:numParams] # get param columns of all rows
    dfpop['simLabel'] = df['simLabel']
    if not pops: pops = list(df.popRates[0].keys())
    for k in pops: dfpop[k] = [1.0*r[k] if k in r else 0 for r in df.popRates] 

    return dfpop


def dfEPSPratios(df):
    dfepsp = df.iloc[:,0:4]
    dfepsp['ratioDiff'] = [x[1][1]-x[1][0] for x in df.EPSPratio]
    dfepsp['ampDiff'] = [x[1][1]-x[1][0] for x in df.EPSPamp]
    dfepsp['peakDiff'] = [x[1][1]-x[1][0] for x in df.EPSPpeak]
    for k in list(df.popRates[0].keys()): dfepsp[k] = [r[k] for r in df.popRates]
    dfepsp['rateDiff'] = [x['PT5B_ZD'] - x['PT5B_1'] for i,x in dfepsp.iterrows()]
    return dfepsp

def dfEPSPamps(df, numParams=4):
    dfepsp = df.iloc[:,0:numParams]
    try: dfepsp['EPSPamp_PTih'] = [x[0][0] for x in df.EPSPamp]
    except: pass
    try: dfepsp['Vpeak_PTih'] = [x[0][0] if isinstance(x, list) else 0 for x in df.EPSPpeak]
    except: pass
    try: dfepsp['EPSPratio_PTih'] = [x[1][0]-1.0 for x in df.EPSPratio] 
    except: pass
    try: dfepsp['EPSPamp_PTzd'] = [x[0][1] for x in df.EPSPamp]
    except: pass
    try: dfepsp['Vpeak_PTzd'] = [x[0][1] for x in df.EPSPpeak]
    except: pass
    try: dfepsp['EPSPratio_PTzd'] = [x[1][1]-1.0 for x in df.EPSPratio] 
    except: pass
    return dfepsp


# ---------------------------------------------------------------------------------
# Wrapper funcs
# ---------------------------------------------------------------------------------
def fIAnalysis(dataFolder, batchLabel, loadAll):
    params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, vars=None, maxCombs=None) 
    plotfI(dataFolder, batchLabel, params, data) #, saveLabel='_I', showLegend=True, pops=['SOM2','PV2']) #['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6'])#)


def popRateAnalysis(dataFolder, batchLabel, loadAll, pars=['IEweights_0','IIweights_0','IEweights_2','IIweights_2'], 
    groupPars = [], vals=['IT2','IT4','IT5A', 'IT5B', 'PT5B','IT6','CT6'], groupStat='first', query=None, plotLine=False):
    # load data 
    var = [('simData','popRates')]
    params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, vars=var, maxCombs=None)
    # params, data = utils.readBatchUpdatePopRates(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, timeRange=[1000,2000])
    
    # hack to fix bug after py3
    for i in range(len(params)):
        if isinstance(params[i]['label'], list):
            params[i]['label'] = '%s_%s' % (params[i]['label'][0], params[i]['label'][1])    

    # convert to pandas and add pop Rates
    df1 = utils.toPandas(params, data)
    dfpop = dfPopRates(df1, 8, pops=vals)

    if query: dfpop = dfpop.query(query)

    print(dfpop)
    # Subtrat num spks for ihGbar=0.0
    zdcomp = 0
    if zdcomp:
        #pd.options.mode.chained_assignment = None
        #dfpop.is_copy = False
        dfzd = dfpop.query('ihGbar == 0.0')
        missing = 0
        for i in range(0, len(dfpop)):
            row = dfpop.iloc[i]
            #zdval = dfzd.query('groupWeight==' + str(row['groupWeight']) + ' and ihlkc==' + str(row['ihlkc']) + ' and ihlke==' + str(row['ihlke'])).iloc[0]['PT5B']
            try:
                zdval = dfzd.query('axonNa==' + str(row['axonNa']) + ' and ihlkc==' + str(row['ihlkc']) + ' and ihlke==' + str(row['ihlke'])).iloc[0]['PT5B']
            except:
                zdval = 0.0
                missing+=1
            #dfpop.iloc[i,dfpop.columns.get_loc('PT5B')] = row['PT5B'] - zdval
            dfpop.iloc[i,dfpop.columns.get_loc('PT5B')] = (zdval/row['PT5B'])-1.0
        print('MISSING: '+str(missing))
        

    # plot param grids
    if len(pars) == 1:
        if len(groupPars) > 0:
            dfpiv=plot2Dparams(dfpop, par1=pars[0], par2=groupPars, val=vals, valLabel='rate (Hz)',  groupStat='first', normRange=1, line=True)
        else:
            dfpiv=plot1Dparams(dfpop, par1=pars[0], val=vals, valLabel='rate (Hz)', groupStat='first', normRange=1, line=True)
    elif len(pars) == 2:
        dfpiv=plot2Dparams(dfpop, par1=pars[0], par2=pars[1], val=vals, valLabel='rate (Hz)',  groupStat='first', normRange=1, line=plotLine)
    elif len(pars) == 3:
        dfpiv=plot3Dparams(dfpop, par1=pars[0], par2=pars[1], par3=pars[2], val=vals, valLabel='rate (Hz)',  groupStat='first', normRange=1)
    elif len(pars) == 4:
        dfpiv=plot4Dparams(dataFolder, batchLabel, dfpop, par1=pars[0], par2=pars[1], par3=pars[2], par4=pars[3], val=vals, valLabel='rate (Hz)',  groupStat=groupStat, normRange=1)
              
    return dfpop, dfpiv           


def ihEPSPAnalysis(dataFolder, batchLabel, loadAll, query=None, 
                    vals = ['ratioDiff', 'ampDiff', 'peakDiff', 'rateDiff'], 
                    pars = ['groupWeight','excTau2Factor','groupRate','ihFactor2'],
                    groupPars = [],
                    zdcomp = False, 
                    plotLine=True):
    # load data
    varRequiredForVal = {'Vpeak_PTih': ('simData','EPSPpeak'), 
                        'EPSPamp_PTih': ('simData', 'EPSPamp'),
                        'ratioDiff': ('simData','EPSPratio'), 
                        'ampDiff': ('simData','EPSPamp'), 
                        'peakDiff': ('simData','EPSPpeak'), 
                        'rateDiff': ('simData','popRates')}
    var = [v for k,v in varRequiredForVal.items() if k in vals]
    params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, vars=var, maxCombs=None) #, listCombs=filename)

    # convert to pandas and add EPSP data
    df1 = utils.toPandas(params, data)
    if 'ratioDiff' in vals: dfepsp = dfEPSPratios(df1)
    elif 'EPSPamp_PTih' in vals or 'EPSPratio_PTih' in vals or 'Vpeak_PTih' in vals: 
        print(df1)
        dfepsp = dfEPSPamps(df1, 6)
        print(dfepsp)

    if query: 
        df2 = dfepsp.query(query)
    else:
        df2 = dfepsp

    #zdcomp = 1
    if zdcomp:
        # df2['Vpeak_PTih'] = [row['Vpeak_PTih'] - peak0ih[row['groupWeight']] for i,row in df2.iterrows()]

        # Subtrat EPSP amps of ihGbar=0.0
        # dfzd = dfepsp.query('ihGbar == 0.0')
        # for i,row in df2.iterrows():
        #     zdval = dfzd.query('groupWeight==' + str(row['groupWeight']) + ' and dendNa==' + str(row['dendNa']) + \
        #         ' and axonNa==' + str(row['axonNa']) + ' and axonRa==' + str(row['axonRa'])).iloc[0]['Vpeak_PTih']
        #     df2.iloc[i, df2.columns.get_loc('Vpeak_PTih')] = row['Vpeak_PTih'] - zdval

        # Subtrat EPSP amps of ihGbar=0.0
        dfzd = df2.query('ihGbar == 0.0')
        #for i,row in df2.iterrows():
        missing = 0
        for i in range(0, len(df2)):
            row = df2.iloc[i]
            #print row
            dfzd2 = dfzd.query('groupWeight > (' + str(row['groupWeight']) + '- 2e-5)' + 'and groupWeight < (' + str(row['groupWeight']) + '+ 2e-5)' +
                ' and axonNa==' + str(row['axonNa']) + ' and gpas==' + str(row['gpas']) + ' and epas==' + str(row['epas']))
            #print dfzd2
            try:
                zdval = dfzd2.iloc[0, dfzd2.columns.get_loc('Vpeak_PTih')]
            except:
                zdval = 0.0
                missing+=1
            df2.iloc[i, df2.columns.get_loc('Vpeak_PTih')] = row['Vpeak_PTih'] - zdval
        print('MISSING: '+str(missing))

    print(df2)
    # plot param grids
    if len(pars) == 1:
        if len(groupPars) > 0:
            plot2Dparams(df2, par1=pars[0], par2=groupPars, val=vals, valLabel='(mV)',  groupStat='first', normRange=1, line=True)
        else:
            plot1Dparams(df2, par1=pars[0], val=vals, valLabel='(mV)',  groupStat='first', normRange=1, line=plotLine)
    if len(pars) == 2:
        plot2Dparams(df2, par1=pars[0], par2=pars[1], val=vals, valLabel='(mV)',  groupStat='first', normRange=1, line=plotLine)
    elif len(pars) == 3:
        plot3Dparams(df2, par1=pars[0], par2=pars[1], par3=pars[2], val=vals, valLabel='(mV)',  groupStat='first', normRange=0)
    elif len(pars) == 4:
        plot4Dparams(dataFolder, batchLabel, df2, par1=pars[0], par2=pars[1], par3=pars[2], par4=pars[3], val=vals, valLabel='(mV)',  groupStat='first', normRange=1)
               

def extractRates(dataFolder, batchLabel, include=None):
    from os import listdir
    from os.path import isfile, join
    
    jsonFolder = batchLabel 
    path = dataFolder+batchLabel #+'/noDepol'
    onlyFiles = [f for f in listdir(path) if isfile(join(path, f)) and not f.endswith('batch.json') and not f.endswith('cfg.json')]

    #outfiles = [f for f in onlyFiles if f.startswith(batchLabel+'_0') and f.endswith('.json')] 
    outfiles = [f for f in onlyFiles if f.endswith('.json')] 

    pops = ['IT5A', 'PT5B', 'IT2']
    avgs = {}
    peaks = {}            
    include = {}
    include['rates'] = pops

    with open('../sim/cells/popColors.pkl', 'r') as fileObj: popColors = pickle.load(fileObj)['popColors']
    saveData = {}
    counter=0
    for outfile in outfiles:    
        try: 
            filename = dataFolder+jsonFolder+'/'+outfile
            print(filename)     

            sim,data,out=utils.plotsFromFile(filename, raster=0, stats=0, rates=1, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0, 
                timeRange=[600,1400], include=include, textTop='', popColors=popColors, orderBy='gid')
            
            for ipop, pop in enumerate(pops):
                avgs[pop] = [float(out[2][0][ipop]), float(out[2][1][ipop])]
                peaks[pop] = [float(out[3][0][ipop]), float(out[3][1][ipop])]

            saveData[outfile] = {'avgs': dict(avgs), 'peaks': dict(peaks)}
            filename = dataFolder+jsonFolder+'/'+'ratesData.json'
            with open(filename, 'w') as fileObj: json.dump(saveData, fileObj)
        except:
            pass


def plotSeedsRates(dataFolder, batchLabel, include=None):
    import copy
    # load
    loadData={}
    batchLabels = [batchLabel] #['v49_batch24', 'v49_batch28'] # [batchLabel] 
    for batchLabel in batchLabels:
        filename = dataFolder+batchLabel+'/'+'ratesData.json'
        with open(filename, 'r') as fileObj: loadData[batchLabel] = json.load(fileObj)

    # iterate data and calculate normalized % increase: (post-pre)/pre
    # calculate absolute increase in Hz (remove baseline rate) -- ih affects overall PT activity
    rows = []
    cols = ['longPop', 'pop', 'connSeed', 'stimSeed', 'ih',  'incAvg', 'incPeak']
    longPops = ['None', 'TPO', 'TVL', 'S2', 'M2']
    pops = ['IT2', 'IT5A', 'PT5B'] 
    ihs = [0.3, 1.0] #[0.3, 0.4, 0.5, 1.0]
    avgs = {k1: None for k1 in  pops}
    peaks = {k1: None for k1 in  pops}
    baseAvgs = None
    basePeaks = None
    fileRoot = dataFolder+batchLabel+'/'+batchLabel
    missing = []
    for iih, ih in enumerate(ihs):
        for ilong,longPop in enumerate(longPops):
            for connSeed in range(5):
                for stimSeed in range(5):
                    try:
                        # if ih in [0.5,1.0]:
                        #     batchLabel = 'v49_batch24'
                        #     ihOffset = -2
                        # else:
                        #     batchLabel = 'v49_batch28'
                        #     ihOffset = 0
                        ihOffset = 0
                        fileRoot = dataFolder+batchLabel+'/'+batchLabel
                        label = '%s_%d_%d_%d_%d.json'%(batchLabel,iih+ihOffset,connSeed,stimSeed,ilong)
                        avgs = loadData[batchLabel][label]['avgs']
                        peaks = loadData[batchLabel][label]['peaks']

                        if longPop == 'None':
                            baseAvgs, basePeaks = copy.deepcopy(avgs), copy.deepcopy(peaks)
                        for pop in pops:
                            avg, baseAvg = avgs[pop], baseAvgs[pop]
                            incAvg = avg[1]  - baseAvg[1] # / baseAvg[1] if baseAvg[1] > 0 else None # avg[1] - avg[0] if avg[0] > 0 else None
                            peak, basePeak = peaks[pop], basePeaks[pop]
                            incPeak = peak[1] - baseAvg[0] # / basePeak[1] if basePeak[1] > 0 else None  # peak[1] - peak[0] if peak[0] > 0 else None
                            rows.append([longPop, pop, connSeed, stimSeed, ih, incAvg, incPeak])
                    except Exception as e:
                        print(e)
                        missing.append(label)
    print(missing) 
    print(len(missing))
    
    #batchLabel = 'v49_batch28'
    fileRoot = dataFolder+batchLabel+'/'+batchLabel
    
    # boxplots
    df = pd.DataFrame(rows, columns=cols) 
    for longPop in longPops:
        for measure in ['incAvg', 'incPeak']:
            dflong = df.query('longPop == "%s"'%(longPop))
            plt.figure()
            sb.stripplot(x='pop', y=measure, hue='ih', data=dflong, palette='Set3', jitter=True, split=True,linewidth=1,edgecolor='gray')
            ax=sb.boxplot(x='pop', y=measure, hue='ih', data=dflong, palette='Set3', showfliers=False)

            handles, labels = ax.get_legend_handles_labels()
            l = plt.legend(handles[0:1], labels[0:1], title='PT ih', loc=2, bbox_to_anchor=(1.02, 1), borderaxespad=0.)  #remove duplicate legend; 
            #plt.ylim(-5,20)
            #plt.tight_layout()
            plt.savefig(fileRoot+'_%s_%s_boxplot_basePeak_sub.png'%(longPop, measure))


def plotSimultLongRates(dataFolder, batchLabel, include=None):
    import copy
    # load
    filename = dataFolder+batchLabel+'/'+'ratesData.json'
    with open(filename, 'r') as fileObj: loadData = json.load(fileObj)

    # iterate data and calculate normalized % increase: (post-pre)/pre
    # calculate absolute increase in Hz (remove baseline rate) -- ih affects overall PT activity
    rows = []
    cols = ['long2Pop', 'pop', 'interval', 'ih',  'incAvg', 'incPeak']
    long2Pops = ['TPO+TVL', 'TVL+TPO', 'TVL+S2', 'S2+TVL', 'S2+M2', 'M2+S2'] 
    pops = ['IT2', 'IT5A', 'PT5B'] 
    intervals = list(np.arange(1000, 1220, 20))
    ihs = ['high', 'low']

    avgs = {k1:{k2:None for k2 in ihs} for k1 in pops}
    peaks = {k1:{k2:None for k2 in ihs} for k1 in pops}

    fileRoot = dataFolder+batchLabel+'/'+batchLabel
    missing = []
    for ilong,long2Pop in enumerate(long2Pops):
        for iinterval,interval in enumerate(intervals):
            # try:
            label = '%s_%d_%d_%d_1.json'%(fileRoot,ilong,ilong,iinterval)
            print(label)
            [avgs['PT5B']['low'], peaks['PT5B']['low'],\
             avgs['PT5B']['high'], peaks['PT5B']['high'], \
             avgs['IT5A']['low'], peaks['IT5A']['low'], \
             avgs['IT5A']['high'], peaks['IT5A']['high'],\
             avgs['IT2']['low'], peaks['IT2']['low'], \
             avgs['IT2']['high'], peaks['IT2']['high']] = loadData[label]

            for pop in pops:
                for ih in ihs:
                    avg = avgs[pop][ih]
                    incAvg = (avg[1] - avg[0]) # / baseAvg[1] if baseAvg[1] > 0 else None # avg[1] - avg[0] if avg[0] > 0 else None
                    peak = peaks[pop][ih]
                    incPeak = (peak[1] - avg[0]) # / basePeak[1] if basePeak[1] > 0 else None  # peak[1] - peak[0] if peak[0] > 0 else None
                    rows.append([long2Pop, pop, interval, ih, incAvg, incPeak])
            # except Exception as e:
            #     print e
            #    missing.append(label)
    print(missing) 
    print(len(missing))
    
    # boxplots
    long2PopGroups = [['TPO+TVL', 'TVL+TPO'], ['TVL+S2', 'S2+TVL'], ['S2+M2', 'M2+S2']] 
    df = pd.DataFrame(rows, columns=cols) 
    for pop in pops:
        for long2PopGroup in long2PopGroups:
            for measure in ['incAvg', 'incPeak']:
                plt.figure()
                dflong = df.query('pop==@pop and long2Pop == @long2PopGroup')
                print(dflong)
                dfpiv = pd.pivot_table(dflong, index='interval', columns=['ih', 'long2Pop'], values=measure)
                blue = [0.32, 0.44, 0.67]
                green = [0.41, 0.64, 0.43]
                dfpiv.plot(color= [blue,green,blue,green], style=['--','--','-','-'], marker='o')
                plt.tight_layout()
                plt.savefig(fileRoot+'_%s_%s_%s_plot.png'%(pop, long2PopGroup, measure))


def extractRatePSD(dataFolder, batchLabel, simLabel=None, include=None, copyFiles=False):
    from os import listdir
    from os.path import isfile, join
    from netpyne import sim,specs
    
    jsonFolder = batchLabel 
    path = dataFolder+batchLabel #+'/noDepol'
    onlyFiles = [f for f in listdir(path) if isfile(join(path, f)) and not f.endswith('batch.json') and not f.endswith('cfg.json')]

    if type(simLabel) is list:
        outfiles = [f for f in onlyFiles if any([f.endswith(sl+'.json') for sl in simLabel])] 
    elif type(simLabel) is '':
        outfiles = [f for f in onlyFiles if f.endswith(simLabel+'.json')]
    else:
        outfiles = [f for f in onlyFiles if f.endswith('.json')] 
  
    with open('../sim/cells/popColors.pkl', 'r') as fileObj: popColors = pickle.load(fileObj)['popColors']
    
    saveData = {}
    counter=0
    
    include = ['IT2', 'IT5A', 'PT5B']

    saveFilename = dataFolder+jsonFolder+'/'+'ratesPSDData.json'

    for outfile in outfiles:
        try:
            filename = dataFolder+jsonFolder+'/'+outfile
            print(filename)   
            with open(filename, 'r') as fileObj:
                data = json.load(fileObj, object_pairs_hook=specs.OrderedDict)
            cfg = specs.SimConfig(data['simConfig'])
            cfg.createNEURONObj = False
            sim.initialize()  # create network object and set cfg and net params
            sim.loadAll('', data=data)
            sim.setSimCfg(cfg)
            sim.net.createPops()     
            sim.net.createCells()
            sim.gatherData() 

            sim.allSimData = data['simData']

            timeRange = [cfg.NetStim1['start'], cfg.NetStim1['start']+1000.0]

            fig, allSignal, allPower, allFreqs = sim.analysis.plotRatePSD(include=include, timeRange=timeRange, ylim=[-40,20], Fs=100, smooth=16, showFig=0, saveFig=1, popColors=popColors, figSize=(8,5))

            saveData[filename] = list(allSignal[0])
            
            with open(saveFilename, 'w') as fileObj: json.dump(saveData, fileObj)
        except:
            pass

    saveData['allFreqs'] = allFreqs
    with open(saveFilename, 'w') as fileObj: json.dump(data, fileObj)


def plotResonantFreq():
    freqs = [4,8,12,16,20,24,28,32,36,40]
    starts = [500, 550]
    ihGbars = [0.5, 1.0]
    percents = [10,20,30]
    include = ['IT2', 'IT5A', 'PT5B']

    rows = []
    cols = ['freq', 'start', 'ih', 'percent', 'pop', 'signalPSD']

    findices = [10, 20, 31, 41, 51, 61, 72, 82, 92, 102]  # use interpolate function

    for ifreq, freq in enumerate(freqs):
        for istart, start in enumerate(starts):
            for iih, ih in enumerate(ihGbars):
                for ipercent, percent in enumerate(percents):

            #print allSignal
                    for ipop, pop in enumerate(include):
                        signalPSD = allSignal[ipop]
                        rows.append([freq, start, ih, percent, pop, signalPSD])

    # labels = ['start=500ms; low ih', 'start=550ms; low ih', 'start=500ms; high ih', 'start=550ms; high ih']

    df = pd.DataFrame(rows, columns=cols) 




