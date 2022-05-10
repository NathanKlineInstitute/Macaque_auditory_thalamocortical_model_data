"""
paper_figs1-5.py 

Code to generate Figures 1 to 5 in the following publication:

Dura-Bernal S, Griffith EY, Barczak A, O’Connell MN, McGinnis T, Schroeder C, Lytton WW, Lakatos P, Neymotin SA. (2022) 
Data-driven multiscale model of macaque auditory thalamocortical circuits reproduces in vivo dynamics. BioRxiv 2022.02.03.479036

"""

#import seaborn as sb
import os
import pickle

import IPython as ipy
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from netpyne.analysis.utils import colorList
from netpyne.support.scalebar import add_scalebar


# ---------------------------------------------------------------------------------------------------------------
# Population params
allpops = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM', 'IC']
popColor = {}
for i,pop in enumerate(allpops):
    popColor[pop] = colorList[i]


def loadSimData(dataFolder, batchLabel, simLabel):
    ''' load sim file'''
    root = dataFolder+batchLabel+'/'
    sim,data,out = None, None, None
    if isinstance(simLabel, str): 
        filename = root+simLabel+'.pkl'
        print(filename)
        sim,data,out = utils.plotsFromFile(filename, raster=0, stats=0, rates=0, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0)
    
    return sim, data, out, root

def axisFontSize(ax, fontsize):
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 

def plot_conn():

    with open('../conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    pmat = connData['pmat']
    lmat = connData['lmat']
        
    popsPre = allpops
    popsPost = allpops  # NOTE: not sure why CT5B and PT5B order was switched
    
    connMatrix = np.zeros((len(popsPre), len(popsPost)))

    d = 50
    for ipre, pre in enumerate(popsPre):
        for ipost, post in enumerate(popsPost):
            print(pre,post)
            try:
                if pre in lmat and post in lmat[pre]:
                    connMatrix[ipre, ipost] = pmat[pre][post] * np.exp(-d / lmat[pre][post])** 2
                else:
                    connMatrix[ipre, ipost] = pmat[pre][post]
            except:
                connMatrix[ipre, ipost] = 0.0
                #connMatrix[ipre, ipost] = pmat[pre][post]

    # font
    fontsiz = 14
    plt.rcParams.update({'font.size': fontsiz})

    # ----------------------- 
    # conn matrix full

    vmin = np.nanmin(connMatrix)
    vmax = np.nanmax(connMatrix)
        
    plt.figure(figsize=(12, 12))
    h = plt.axes()
    plt.imshow(connMatrix, interpolation='nearest', cmap='viridis', vmin=vmin, vmax=vmax)  #_bicolormap(gap=0)


    ipopBoundaries = [1, 6, 11, 17, 23, 30, 40]
    for ipop, pop in enumerate(popsPre):
        if ipop in ipopBoundaries: # thicker, brighter, dotted lines for layer boundaries
            plt.plot(np.array([0,len(popsPost)])-0.5,np.array([ipop,ipop])-0.5,'-',c=(0.8,0.8,0.8), lw=3)
        else:
            plt.plot(np.array([0,len(popsPost)])-0.5,np.array([ipop,ipop])-0.5,'-',c=(0.7,0.7,0.7))
    for ipop, pop in enumerate(popsPost):
        if ipop in ipopBoundaries: # thicker, brighter, dotted lines for layer boundaries
            plt.plot(np.array([ipop,ipop])-0.5,np.array([0,len(popsPre)])-0.5,'-',c=(0.8,0.8,0.8), lw=3)
        else:
            plt.plot(np.array([ipop,ipop])-0.5,np.array([0,len(popsPre)])-0.5,'-',c=(0.7,0.7,0.7))

    ipop = 36 # thal boundary
    plt.plot(np.array([0,len(popsPost)])-0.5,np.array([ipop,ipop])-0.5,'-', c='orange', lw=3)
    plt.plot(np.array([ipop,ipop])-0.5,np.array([0,len(popsPre)])-0.5,'-', c='orange', lw=3)


    # Make pretty
    h.set_yticks(list(range(len(popsPre))))
    h.set_xticks(list(range(len(popsPost))))
    h.set_yticklabels(popsPre)
    h.set_xticklabels(popsPost, rotation=90)
    h.xaxis.set_ticks_position('top')

    plt.grid(False)

    vmax = 0.5
    clim = [vmin, vmax]
    plt.clim(clim[0], clim[1])
    plt.colorbar(label='probability', shrink=0.8) 
    plt.xlabel('Post')
    plt.ylabel('Pre')
    title = 'Connection probability matrix' # empirical
    h.xaxis.set_label_coords(0.5, 1.11)
    plt.title(title, y=1.12, fontWeight='bold')
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.00)

    filename = 'Full_Allen_custom_prob_conn_empirical.png'
    plt.savefig('../conn/'+filename, dpi=300)



def fig_raster(batchLabel, simLabel, timeRange=[500,1500]):
    dataFolder = '../data/'

    sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)
    
    #raster
    include = allpops
    orderBy = ['pop'] #, 'y']

    # use these params if short 1-sec raster 
    fig1 = sim.analysis.plotRaster(include=['allCells'], timeRange=timeRange, labels='legend', 
        popRates=False, orderInverse=True, lw=0, markerSize=12, marker='.',  
        showFig=0, saveFig=0, figSize=(9*0.95, 13*0.9), orderBy=orderBy)# 
    ax = plt.gca()

    [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner
    
    plt.xticks(timeRange, ['0', '1'])
    plt.yticks([0, 5000, 10000], [0, 5000, 10000])

    plt.ylabel('Neuron ID') #Neurons (ordered by NCD within each pop)')
    plt.xlabel('Time (s)')
    
    plt.title('')
    filename='%s%s_raster_%d_%d_%s.png'%(root, simLabel, timeRange[0], timeRange[1], orderBy)
    plt.savefig(filename, dpi=300)


def fig_stats(batchLabel, simLabel, timeRange):

    dataFolder = '../data/'
    
    popColor['SOM'] = popColor['SOM2'] 
    popColor['PV'] = popColor['PV2'] 
    popColor['VIP'] = popColor['VIP2'] 
    popColor['NGF'] = popColor['NGF2'] 
    

    sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

    #timeRange = [1000, 6000] #[2000, 4000]

    statPops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6',
                ('SOM2','SOM3','SOM4','SOM5A','SOM5B','SOM6'),
                ('PV2','PV3','PV4','PV5A','PV5B','PV6'),
                ('VIP2', 'VIP3', 'VIP4', 'VIP5A', 'VIP5B', 'VIP6'),
                ('NGF1', 'NGF2', 'NGF3', 'NGF4','NGF5A','NGF5B','NGF6'),
                'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']
    
    labels = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6',
            'SOM', 'PV', 'VIP', 'NGF', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']

    colors = [popColor[p] for p in labels]
    
    xlim = [0,40]

    fig1,modelData = sim.analysis.plotSpikeStats(include=statPops, stats=['rate'], timeRange=timeRange, includeRate0=False,
        showFig=0, saveFig=0, figSize=(8.5, 13))
    modelData = modelData['statData']

    plt.figure(figsize=(6*2, 6.5*2))
    meanpointprops = dict(marker = (5, 1, 0), markeredgecolor = 'black', markerfacecolor = 'white')
    fontsiz = 20    

    bp=plt.boxplot(modelData[::-1], labels=labels[::-1], notch=False, sym='k+', meanprops=meanpointprops, whis=1.5, widths=0.6, vert=False, showmeans=True, showfliers=False, patch_artist=True)  #labels[::-1] #positions=np.array(range(len(statData)))+0.4,
    plt.xlabel('Rate (Hz)', fontsize=fontsiz)
    plt.ylabel('Population', fontsize = fontsiz)
    plt.subplots_adjust(left=0.3,right=0.95, top=0.9, bottom=0.1)

    icolor=0
    borderColor = 'k'
    for i in range(0, len(bp['boxes'])):
        icolor = i
        bp['boxes'][i].set_facecolor(colors[::-1][icolor])
        bp['boxes'][i].set_linewidth(2)
        # we have two whiskers!
        bp['whiskers'][i*2].set_color(borderColor)
        bp['whiskers'][i*2 + 1].set_color(borderColor)
        bp['whiskers'][i*2].set_linewidth(2)
        bp['whiskers'][i*2 + 1].set_linewidth(2)
        bp['medians'][i].set_color(borderColor)
        bp['medians'][i].set_linewidth(3)
        for c in bp['caps']:
            c.set_color(borderColor)
            c.set_linewidth(1)

    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis='y', direction='out')
    ax.grid(axis='x', color="0.9", linestyle='-', linewidth=1)
    ax.set_axisbelow(True)
    if xlim: ax.set_xlim(xlim)
    axisFontSize(ax, fontsiz)

    plt.title('')
    filename='%s%s_stats_%d_%d.png'%(root, simLabel, timeRange[0], timeRange[1])
    plt.savefig(filename, dpi=300)



def fig_traces(batchLabel, simLabel, timeRange):
    dataFolder = '../data/'

    sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)
    #popParamLabels = list(data['simData']['popRates'])

    allpops = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM', 'IC']

    firingpops = ['NGF1', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3',  'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'CT5B',  'PV5B', 'VIP5B', 'NGF5B',  'VIP6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI']

    firingpops = ['NGF1', 'PV2', 'NGF2', 'IT3',  'SOM3',  'VIP3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'IT5A', 'CT5A', 'SOM5A', 'VIP5A', 'IT5B', 'CT5B',  'PV5B', 'NGF5B', 'TC', 'HTC', 'IRE', 'TI']

    timeRanges = [timeRange] #[[3000, 5000], [5000, 7000], [7000, 9000], [9000,11000]]

    cellNames = data['simData']['V_soma'].keys()
    popCells = {}
    for popName,cellName in zip(allpops,cellNames):
        popCells[popName] = cellName

    fontsiz = 20   
    for timeRange in timeRanges:

        plt.figure(figsize=(9*1.05, 1.2*2)) 
        time = np.linspace(timeRange[0], timeRange[1], 10001)
        plt.ylabel('V (mV)', fontsize=fontsiz)
        plt.xlabel('Time (ms)', fontsize=fontsiz)
        plt.xlim(500, 1500)
        # plt.ylim(-80, -30)
        plt.ylim(-120*len(firingpops),20)
        plt.yticks(np.arange(-120*len(firingpops)+60,60,120), firingpops[::-1], fontsize=fontsiz)
        number = 0
        
        for popName in firingpops: #allpops:
            cellName = popCells[popName]   
            Vt = np.array(data['simData']['V_soma'][cellName][timeRange[0]*10:(timeRange[1]*10)+1])
            plt.plot(time, (Vt-number*120.0), color=popColor[popName]) 
            number = number + 1

        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        filename='%s%s_%d_%d_firingTraces.png'%(root, simLabel, timeRange[0], timeRange[1])
        plt.savefig(filename, facecolor = 'white', bbox_inches='tight' , dpi=300)


def fig_CSD():
    from netpyne import sim
    
    targetFolder = 'data/v34_batch27'

    filenames = ['data/v34_batch27/v34_batch27_%d_%d.pkl' % (iseed, cseed) for iseed in [0] for cseed in [0]]

    # find all individual sim labels whose files need to be gathered
    #filenames = [targetFolder+f for f in os.listdir(targetFolder) if f.endswith('.pkl')]

    layer_bounds= {'L1': 100, 'L2': 160, 'L3': 950, 'L4': 1250, 'L5A': 1334, 'L5B': 1550, 'L6': 2000}

    for filename in filenames:
        sim.load(filename, instantiate=False)
        
        tranges = [[x, x+1000] for x in range(500, 600, 100)]
        for t in tranges:# (2100, 2200,100):    
            sim.analysis.plotCSD(**{
                'spacing_um': 100, 
                'layer_lines': 1, 
                'layer_bounds': layer_bounds, 
                'overlay': 'LFP',
                'timeRange': [t[0], t[1]], 
                'smooth': 30,
                'saveFig': filename[:-4]+'_CSD_LFP_%d-%d' % (t[0], t[1]), 
                'figSize': (4.1,8.2), 
                'dpi': 300, 
                'showFig': 0})
            


def fig_EEG():

    from netpyne import sim
    sim.load('../data/v34_batch53/v34_batch53_0_0_data.pkl', instantiate=False)

    import matplotlib; matplotlib.rcParams.update({'font.size': 16})

    sim.analysis.plotDipole(saveFig='../data/v34_batch53/dipole_all', timeRange=[500,1500], figSize=(9,4.5))

    sim.analysis.plotEEG(sim, saveFig='../data/v34_batch53/EEG_all', timeRange=[500,1500])

def fig_LFP_PSD_matrix():
    import random

    import seaborn as sns
    from scipy.spatial.distance import cdist, pdist
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler as Sc

    dataFiles = ['../data/NHPdata/spont/spont_LFP_PSD/1-bu031032017@os_eye06_20_10sec_allData.pkl', 
                '../data/NHPdata/spont/spont_LFP_PSD/2-ma031032023@os_eye06_20_10sec_allData.pkl', 
                '../data/NHPdata/spont/spont_LFP_PSD/2-rb031032016@os_eye06_20_10sec_allData.pkl', 
                '../data/NHPdata/spont/spont_LFP_PSD/2-rb045046026@os_eye06_20_10sec_allData.pkl',
                '../data/v34_batch57/v34_batch57_10sec_allData.pkl']

    psd = []
    clusters = []
    start = 0
    elec = 0

    for dataFile in dataFiles:            
        NHP = dataFile[:-4]
    
        with open(dataFile, 'rb') as f:
            loadedData = pickle.load(f)
            allData = loadedData['allData']
            freqs = loadedData['allFreqs'][0]

        for itime in range(len(allData)):
            psd.append(allData[itime][elec])

        clusters.append(range(start, start+len(allData)))
        start = start + len(allData)


    labels = ['NHP 1', 'NHP 2', 'NHP 3', 'NHP 4', 'Model', 'Model shuffled']
    df = pd.DataFrame(psd)
    
    # add model shuffled
    model_shuffled = []
    for row in df.iloc[-25:].iterrows():
        rowshuff = list(row[1]) 
        random.shuffle(rowshuff)
        model_shuffled.append(rowshuff)
    clusters.append(range(start, start+len(model_shuffled)))
    
    model_shuffled = pd.DataFrame(model_shuffled)
    df=pd.concat([df, model_shuffled])#.reset_index()
    
    # corr matrix
    dfT = df.T
    corrMatrix = dfT.corr()
    plt.figure()
    sns.heatmap(corrMatrix, cmap='viridis')
    plt.savefig('../data/NHPdata/spont/spont_LFP_PSD/PSD_corr_matrix.png')

    # mean correlation
    meanCorr = {}
    for cluster, label in zip(clusters, labels):
        meanCorr[label] = corrMatrix.iloc[cluster, range(clusters[0][0], clusters[3][-1])].mean().mean()

    # stats
    import scipy
    modelCorr=corrMatrix.iloc[clusters[-2], range(clusters[0][0], clusters[3][-1])]
    modelShuffCorr=corrMatrix.iloc[clusters[-2], range(clusters[0][0], clusters[3][-1])]    
    scipy.stats.mannwhitneyu(modelCorr, modelShuffCorr)    

    # PCA
    dfnorm = df.div(df.max(axis=1), axis=0)   
    #dfnorm = Sc().fit_transform(df) 
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(dfnorm)

    print('Explained variance ratio: ', pca.explained_variance_ratio_)

    fig = plt.figure(figsize=(7, 5))
    colors = ['navy',  'blue', 'purple', 'indigo',  'red', 'limegreen']

    colors = 'bcmygr' 
    colors = ['blue', 'r', 'm', 'y',  'limegreen', 'k']

    for cluster, label, c in zip(clusters, labels, colors): # 'rgbcmykw'):
        plt.scatter(X_pca[cluster, 0], X_pca[cluster, 1], c=c, label=label)

    cluster_points = {}
    for cluster, label, c in zip(clusters, labels, colors): # 'rgbcmykw'):
        cluster_points[label] = np.array(list(zip(X_pca[cluster, 0], X_pca[cluster, 1])))

    # distance between NHP 2 and NHP 3 vs NHP2 vs model
    print(np.mean(cdist(cluster_points['NHP 2'], cluster_points['NHP 3'])))
    print(np.mean(cdist(cluster_points['NHP 2'], cluster_points['Model'])))

    nhps = ['NHP 1', 'NHP 2', 'NHP 3', 'NHP 4']

    # distnace between macaques and model vs macaques and shuffled
    distToModel = 0
    distToModelShuffled = 0
    
    for nhp in nhps:
        distToModel += np.mean(cdist(cluster_points[nhp], cluster_points['Model']))
        distToModelShuffled += np.mean(cdist(cluster_points[nhp], cluster_points['Model shuffled']))

    distToModel /= 4
    distToModelShuffled /= 4

    plt.xlabel('PC 1')
    plt.ylabel('PC 2')
    #plt.xlim(-5,5)
    #plt.ylim(-5,5)
    plt.subplots_adjust(left=0.1,right=0.7, top=0.9, bottom=0.1)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig('../data/NHPdata/spont/spont_LFP_PSD/PSD_norm_PCA.png')


    #fig = plt.figure()
    kmeans = KMeans(n_clusters=6).fit(dfnorm)
    pca_data = pd.DataFrame(X_pca, columns=['PC1','PC2']) 
    pca_data['cluster'] = pd.Categorical(kmeans.labels_)
    sns.scatterplot(x="PC1", y="PC2", hue="cluster", data=pca_data, sizes=2, marker='o')#, legend=False)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.savefig('../data/NHPdata/spont/spont_LFP_PSD/PSD_norm_PCA_kmeans.png')




def fig_CSD_comparison():

    from netpyne import sim

    exp = 1
    model = 1

    #-------------------------
    # Experiment
    if exp:

        from expDataAnalysis import getbandpass, getflayers, loadfile

        # spont
        testFiles = ['1-bu031032017@os_eye06_20.mat', '2-ma031032023@os_eye06_20.mat', '2-rb031032016@os_eye06_20.mat', '2-rb045046026@os_eye06_20.mat']   # CHANGE FILE HERE IF LOOKING AT SPECIFIC MONKEY

        # speech
        testFiles = [sim]

        # spont
        datafile = '../data/NHPdata/spont/2-rb031032016@os_eye06_20.mat'   # LOCAL DIR 
        
        # speech
        datafile = '../data/NHPdata/speech/2-bu037038046@os.mat'  #2-ma033034023@os.mat'# 2-bu037038046@os.mat'  
        datafile_pkl = '../data/NHPdata/speech/2-bu037038046@os_small.pkl'

        vaknin = 0    

        # setup netpyne
        samprate = 11*1e3  # in Hz
        spacing = 100
        sim.initialize()
        sim.cfg.recordStep = 1000./samprate # in ms


        # load reduced file size
        with open(datafile_pkl, 'rb') as f:
            dataLoad = pickle.load(f)

        [sampr, LFP_data, dt, tt, CSD_data, trigtimes] = dataLoad['data_small']

        if datafile_pkl == '../data/NHPdata/speech/2-bu037038046@os_small.pkl':
            LFP_data = getbandpass(LFP_data.T, sampr=samprate, minf=0.5, maxf=300)
        
        dbpath = '../data/NHPdata/spont/21feb02_A1_spont_layers.csv'  # GCP # CHANGE ACCORDING TO MACHINE USED TO RUN
        dbpath = '../data/NHPdata//speech/A1_speech_layers.csv'   
            
        ##### GET LAYERS FOR OVERLAY #####
        s1low,s1high,s2low,s2high,glow,ghigh,i1low,i1high,i2low,i2high = getflayers(datafile,dbpath=dbpath,getmid=False,abbrev=False) # fullPath is to data file, dbpath is to .csv layers file 
        lchan = {}
        lchan['S'] = s2high
        lchan['G'] = ghigh
        lchan['I'] = CSD_data.shape[0]-1 #i2high
        print('s2high: ' + str(s2high))

        depth = LFP_data.shape[0]*100 
        layer_bounds= {'S': (glow*spacing)+spacing, 'G': (i1low*spacing)+spacing, 'I': depth-spacing} #(i2high*spacing)+spacing}  #list(range(s1low, glow)), list(range(glow, i1low)), list(range(i1low, i2high))

        # plot CSD
        smooth = 30

        trigtimes_ms = [tt[idx]*1000 for idx in trigtimes if idx < len(tt)]
        tranges_speech = [[x+0, x+200] for x in trigtimes_ms] #int(CSD_data.shape[1]/samprate*1000)

        for i,trange in enumerate(tranges_speech[30:]):# 
            print(i)
            sim.analysis.plotCSD(**{
                'CSD_data': CSD_data,
                'LFP_input_data': LFP_data.T,
                'spacing_um': 100, 
                'dt': sim.cfg.recordStep,
                'ymax': depth,
                'layer_lines': 1, 
                'layer_bounds': layer_bounds, 
                'overlay': 'LFP',
                'timeRange': [trange[0], trange[1]], 
                'smooth': smooth,
                'vaknin': vaknin,
                'saveFig': datafile[:-4]+'_CSD_LFP_smooth-%d_%d-%d_vaknin-%d' % (smooth, trange[0], trange[1], vaknin), 
                'figSize': (4.1,8.2), 
                'dpi': 300, 
                'showFig': 0})

    # ------------------------------------
    # Model
    if model:
        
        #layer_bounds= {'L1': 100, 'L2': 160, 'L3': 950, 'L4': 1250, 'L5A': 1334, 'L5B': 1550, 'L6': 2000}
        layer_bounds= {'S': 950, 'G': 1250, 'I': 1900}#, 'L6': 2000}
        
        filename = '../data/v34_batch56/v34_batch56_0_0_data.pkl'
        
        sim.load(filename, instantiate=False)

        tranges = [[x, x+200] for x in range(2400, 35000, 100)]
        vaknin = 0
        smooth = 30
        for trange in tranges:
            sim.analysis.plotCSD(**{
                'spacing_um': 100, 
                'layer_lines': 1, 
                'layer_bounds': layer_bounds, 
                'overlay': 'LFP',
                'timeRange': [trange[0], trange[1]], 
                'smooth': smooth,
                'vaknin': vaknin, 
                'saveFig': filename[:-4]+'_CSD_LFP_smooth-%d_%d-%d_vaknin-%d' % (smooth,trange[0], trange[1], vaknin), 
                'figSize': (4.1,8.2), 
                'dpi': 300, 
                'showFig': 0})
    


def fig_optuna_fitness():
    import seaborn as sns

    import optunaAnalysis as oa

    dataFolder = '../data/'
    batchSim = 'v34_batch14'
    loadFromFile = 1
    
    allpops = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'PV4', 'SOM4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']  #, 'IC']
    rateTimeRanges = ['1500_1750', '1750_2000', '2000_2250', '2250_2500']

    # set font size
    plt.rcParams.update({'font.size': 18})

    # get param labelsc
    paramLabels = oa.getParamLabels(dataFolder, batchSim)

    # load evol data from files
    df = oa.loadData(dataFolder, batchSim, pops=allpops, rateTimeRanges=rateTimeRanges, loadStudyFromFile=loadFromFile, loadDataFromFile=loadFromFile)

    trial_fitness_single = 1
    params_fitness_all = 1
    rates_fitness_all = 1
    rates_fitness_single = 1

    # PLOTS
    skipCols = rateTimeRanges

    # plot trial vs fitness
    if trial_fitness_single:
        excludeAbove = 1000

        if excludeAbove:
            df = df[df.value < excludeAbove]

        df.reset_index()
        df.number=range(len(df))

        dfcorr=df.corr('pearson')

        min=df.iloc[0]['value']
        minlist=[]
        for x in df.value:
            if x < min:
                minlist.append(x)
                min = x
            else:
                minlist.append(min)

        param = 'value'
        if not any([skipCol in param for skipCol in skipCols]): 
            print('Plotting scatter of %s vs %s param (R=%.2f) ...' %('fitness', param, dfcorr['value'][param]))
            #df.plot.scatter(param, 'value', s=4, c='number', colormap='viridis', alpha=0.5, figsize=(8, 8), colorbar=False)
            plt.figure(figsize=(8, 8))
            df.plot.scatter('number', param, s=6, c='blue', colormap='jet_r', alpha=0.5, figsize=(8, 8), colorbar=False, vmin=50, vmax=400)

            plt.plot(list(df['value'].rolling(window=10).mean()), c='orange', label='mean')
            plt.plot(minlist, c='red', label='min')
            # f = plt.gcf()
            # cax = f.get_axes()[1]
            # cax.set_ylabel('Fitness Error')
            plt.ylabel('Fitness Error')
            plt.xlabel('Trial')
            plt.ylim([0,1000])
            plt.legend()
            #plt.title('%s vs %s R=%.2f' % ('trial', param.replace('tune', ''), dfcorr['value'][param]))
            plt.savefig('%s/%s/%s_scatter_%s_%s.png' % (dataFolder, batchSim, batchSim, 'trial', param.replace('tune', '')), dpi=300)


    # plot params vs fitness
    if params_fitness_all:
        excludeAbove = 400
        ylim = None
        if excludeAbove:
            df = df[df.value < excludeAbove]

        df2 = df.drop(['value', 'number'], axis=1)
        fits = list(df['value'])
        plt.figure(figsize=(16, 8))

        paramsData = list(df2[paramLabels].items())

        for i, (k,v) in enumerate(paramsData):
            y = v 
            x = np.random.normal(i, 0.04, size=len(y))         # Add some random "jitter" to the x-axis
            s = plt.scatter(x, y, alpha=0.3, c=[int(f-1) for f in fits], cmap='jet_r') #)
        plt.colorbar(label = 'Fitness Error')
        plt.ylabel('Parameter value')
        plt.xlabel('Parameter')
        if ylim: plt.ylim(0, ylim)
        paramLabels = [x.replace('Gain','') for x in paramLabels]
        plt.xticks(range(len(paramLabels)), paramLabels, rotation=75)
        plt.subplots_adjust(top=0.99, bottom=0.30, right=1.05, left=0.05)
        plt.savefig('%s/%s/%s_scatter_params_%s.png' % (dataFolder, batchSim, batchSim, 'excludeAbove-'+str(excludeAbove) if excludeAbove else ''))
        #plt.show()


    # plot rates vs fitness
    if rates_fitness_all:
        excludeAbove = 400
        ylim = None
        if excludeAbove:
            df = df[df.value < excludeAbove]

        df2 = df[allpops]
        fits = list(df['value'])
        plt.figure(figsize=(16, 8))

        paramsData = list(df2[allpops].items())

        for i, (k,v) in enumerate(paramsData):
            y = v
            x = np.random.normal(i, 0.04, size=len(y))         # Add some random "jitter" to the x-axis
            s = plt.scatter(x, y, alpha=0.3, c=[int(f-1) for f in fits], cmap='jet_r') #)
        plt.colorbar(label = 'Fitness Error')
        plt.ylabel('Firing Rate (Hz)')
        plt.xlabel('Population')
        if ylim: plt.ylim(0, ylim)
        plt.xticks(range(len(allpops)), allpops, rotation=75)
        plt.subplots_adjust(top=0.99, bottom=0.30, right=1.05, left=0.05)
        plt.savefig('%s/%s/%s_scatter_rates_all%s.png' % (dataFolder, batchSim, batchSim, 'excludeAbove-'+str(excludeAbove) if excludeAbove else ''))
        #plt.show()


    # plot pop rate vs fitness
    if rates_fitness_single:
        excludeAbove = 1000
        if excludeAbove:
            df = df[df.value < excludeAbove]

        dfcorr=df.corr('pearson')

        for param in ['PV4']:
            if not any([skipCol in param for skipCol in skipCols]): 
                print('Plotting scatter of %s vs %s param (R=%.2f) ...' %('fitness', param, dfcorr['value'][param]))
                df.plot.scatter(param, 'value', s=6, c='blue', colormap='jet_r', alpha=0.5, figsize=(8, 8), colorbar=False, vmin=0, vmax=1000) #, colorbar=False)
                plt.ylabel('Fitness Error')
                plt.xlabel('PV4 Rate (Hz)')
                plt.xlim([0,30])
                plt.savefig('%s/%s/%s_scatter_%s_%s.png' % (dataFolder, batchSim, batchSim, 'fitness', param.replace('tune', '')), dpi=300)




# --------------------------
# Main
# --------------------------
if __name__ == '__main__':
    # Fig 2
    plot_conn()
    
    # Fig 3
    fig_raster('v34_batch27', 'v34_batch27_0_0', timeRange=[500,1500])
    fig_traces('v34_batch27', 'v34_batch27_0_0', timeRange=[500,1500])  
    fig_stats('v34_batch27', 'v34_batch27_0_0', timeRange=[500,1500])
    fig_CSD()
    fig_EEG()
    
    # Fig 4
    fig_CSD_comparison()    

    # Fig 5
    fig_LFP_PSD_matrix()

    # Supp Fig 2
    fig_optuna_fitness()
