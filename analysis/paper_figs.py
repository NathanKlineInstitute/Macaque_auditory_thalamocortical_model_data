"""
paper_figs.py 

Paper figures

Contributors: salvadordura@gmail.com, ericaygriffith@gmail.com
"""

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import pickle
from netpyne.analysis.utils import colorList
import collections
from netpyne import sim, specs
import utils

try:
    basestring
except NameError:
    basestring = str

import IPython as ipy


# ---------------------------------------------------------------------------------------------------------------
# Global population params
allpops = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM', 'IC']
excpops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6', 'TC', 'TCM', 'HTC', 'IC']
inhpops = ['NGF1', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'SOM6', 'PV6', 'VIP6', 'NGF6', 'IRE', 'IREM', 'TI', 'TIM']
supra =  ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3']
gran = ['ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4']
infra = ['IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6']
supraE =  ['IT2', 'IT3']
granE = ['ITP4', 'ITS4']
infraE = ['IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']
supraI =  ['NGF1', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'SOM3', 'PV3', 'VIP3', 'NGF3']
granI = ['SOM4', 'PV4', 'VIP4', 'NGF4']
infraI = ['SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'SOM6', 'PV6', 'VIP6', 'NGF6']
thal = ['TC', 'TCM', 'HTC']
cortexE = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']
cortexAll = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6']

popColor = {}
for i,pop in enumerate(allpops):
    popColor[pop] = colorList[i]


# ---------------------------------------------------------------------------------------------------
def loadSimData(dataFolder, batchLabel, simLabel):
    ''' helper function: load sim data file'''
    root = dataFolder+batchLabel+'/'
    sim,data,out = None, None, None
    if isinstance(simLabel, str): 
        filename = root+simLabel+'.pkl'
        print(filename)
        sim,data,out = utils.plotsFromFile(filename, raster=0, stats=0, rates=0, syncs=0, hist=0, psd=0, traces=0, grang=0, plotAll=0)
    
    return sim, data, out, root

# ---------------------------------------------------------------------------------------------------
def axisFontSize(ax, fontsize):
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize) 


# ---------------------------------------------------------------------------------------------------
def plot_empirical_conn():

    with open('../model/conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)
    pmat = connData['pmat']
    lmat = connData['lmat']
        
    popsPre = allpops
    popsPost = allpops  

    connMatrix = np.zeros((len(popsPre), len(popsPost)))

    d = 50  # for distance-dependent prob conns, assume fixed distance of 50 um
    for ipre, pre in enumerate(popsPre):
        for ipost, post in enumerate(popsPost):
            try:
                if pre in lmat and post in lmat[pre]:
                    connMatrix[ipre, ipost] = pmat[pre][post] * np.exp(-d / lmat[pre][post])** 2
                else:
                    connMatrix[ipre, ipost] = pmat[pre][post]
            except:
                connMatrix[ipre, ipost] = 0.0

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
    plt.colorbar(label='probability', shrink=0.8) #.set_label(label='Fitness',size=20,weight='bold')
    plt.xlabel('Post')
    plt.ylabel('Pre')
    title = 'Connection probability matrix' # empirical
    h.xaxis.set_label_coords(0.5, 1.11)
    plt.title(title, y=1.12, fontWeight='bold')
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.00)

    filename = 'fig2_prob_conn_matrix_empirical.png'
    plt.savefig('figs/'+filename, dpi=600)


# ---------------------------------------------------------------------------------------------------
def fig_raster(batchLabel, simLabel, timeRange=[500,1500]):
    dataFolder = '../data/'
    sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

    #raster
    include = allpops
    orderBy = ['pop'] #, 'y']
    fontsize = 16

    filename='%s%s_raster_%d_%d_%s.png'%(root, simLabel, timeRange[0], timeRange[1], orderBy)

    fig1 = sim.plotting.plotRaster(include=include, timeRange=timeRange, 
        popRates=False, orderInverse=True, linewidth=0, s=16, marker='.',  
        showFig=0, saveFig=filename, figSize=(9, 12),  orderBy=orderBy)# 

    ax = plt.gca()

    [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

    plt.xticks([1000, 1500, 2000], ['1.0', '1.5', '2.0'], fontsize=fontsize, fontname='Arial')
    plt.yticks([0, 2000, 4000, 6000, 8000, 10000, 12000], [0, 2000, 4000, 6000, 8000, 10000, 12000], fontsize=fontsize, fontname='Arial')
    plt.xlim(timeRange[0], timeRange[1])
    plt.ylim(12908, 0)
    plt.ylabel('Neuron ID', fontsize=fontsize, fontname='Arial') #Neurons (ordered by NCD within each pop)')
    plt.xlabel('Time (s)', fontsize=fontsize, fontname='Arial')
    
    leg = ax.get_legend()
    leg.set_title('')
    
    leg.set_bbox_to_anchor((1.05, 1.0))

    for text in leg.get_texts():
        text.set_fontsize(fontsize)
        text.set_fontname('Arial')

    plt.title('')
    plt.margins(x=0,y=0)

    plt.subplots_adjust(bottom=0.10, top=0.99, right=0.85, left=0.13)

    #filename='%s%s_raster_%d_%d_%s.png'%(root, simLabel, timeRange[0], timeRange[1], orderBy)
    filename='figs/fig3_spike_raster.png'
    plt.savefig(filename, dpi=300)


# ---------------------------------------------------------------------------------------------------
def fig_stats(batchLabel, simLabel, timeRange):

    dataFolder = '../data/'
    
    popColor['SOM'] = popColor['SOM2'] 
    popColor['PV'] = popColor['PV2'] 
    popColor['VIP'] = popColor['VIP2'] 
    popColor['NGF'] = popColor['NGF2'] 
    

    sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)

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

    timeRangeStats = [1000, 11000]
    fig1,modelData = sim.analysis.plotSpikeStats(include=statPops, stats=['rate'], timeRange=[1000, 11000], includeRate0=False,
        showFig=0, saveFig=0, figSize=(8.5, 13))
    modelData = modelData['statData']

    plt.figure(figsize=(9*0.8, 7*0.8))
    meanpointprops = dict(marker = (5, 1, 0), markeredgecolor = 'black', markerfacecolor = 'white')
    fontsize = 15

    bp=plt.boxplot(modelData[::-1], labels=labels[::-1], notch=False, sym='k+', meanprops=meanpointprops, whis=1.5, widths=0.6, vert=False, showmeans=True, showfliers=False, patch_artist=True)  #labels[::-1] #positions=np.array(range(len(statData)))+0.4,
    plt.xlabel('Rate (Hz)', fontsize=fontsize)
    plt.ylabel('', fontsize = fontsize)
    plt.subplots_adjust(left=0.1,right=0.95, top=0.95, bottom=0.1)

    icolor=0
    borderColor = 'k'
    for i in range(0, len(bp['boxes'])):
        icolor = i
        bp['boxes'][i].set_facecolor(colors[::-1][icolor])
        bp['boxes'][i].set_linewidth(2-0.5)
        # we have two whiskers!
        bp['whiskers'][i*2].set_color(borderColor)
        bp['whiskers'][i*2 + 1].set_color(borderColor)
        bp['whiskers'][i*2].set_linewidth(2-0.5)
        bp['whiskers'][i*2 + 1].set_linewidth(2-0.5)
        bp['medians'][i].set_color(borderColor)
        bp['medians'][i].set_linewidth(3-0.5)
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
    axisFontSize(ax, fontsize)

    plt.title('')
    #filename='%s%s_stats_%d_%d.png'%(root, simLabel, timeRange[0], timeRange[1])
    filename='figs/fig3_spike_stats.png'
    plt.savefig(filename, dpi=300)


# ---------------------------------------------------------------------------------------------------
def fig_traces(batchLabel, simLabel, timeRange):
    dataFolder = '../data/'

    sim, data, out, root = loadSimData(dataFolder, batchLabel, simLabel)
    allpops = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3',  'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B',  'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM', 'IC']
    firingpops = ['PV2', 'NGF2', 'IT3',  'SOM3',  'VIP3', 'ITP4', 'ITS4', 'PV4', 'IT5A', 'CT5A', 'SOM5A', 'VIP5A', 'IT5B', 'CT5B',  'PV5B', 'NGF5B', 'TC', 'HTC', 'IRE', 'TI']
    timeRanges = [timeRange]
    fontsize = 15

    cellNames = data['simData']['V_soma'].keys()
    popCells = {}
    for popName,cellName in zip(allpops,cellNames):
        popCells[popName] = cellName

    for timeRange in timeRanges:

        plt.figure(figsize=(8, 7)) 
        time = np.linspace(timeRange[0], timeRange[1], 10001)
        plt.ylabel('V (mV)', fontsize=fontsize)
        plt.xlabel('Time (s)', fontsize=fontsize)
        plt.xlim(1000, 2000)
        # plt.ylim(-80, -30)
        plt.ylim(-120*len(firingpops),20)
        plt.xticks([1000, 1500, 2000], ['1.0', '1.5', '2.0'], fontsize=fontsize)
        plt.yticks(np.arange(-120*len(firingpops)+60,60,120), firingpops[::-1], fontsize=fontsize)

        plt.margins(x=0,y=0)

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

        #filename='%s%s_%d_%d_firingTraces.png'%(root, simLabel, timeRange[0], timeRange[1])
        filename='figs/fig3_voltage_traces.png'
        plt.savefig(filename, facecolor = 'white', bbox_inches='tight' , dpi=300)



# ---------------------------------------------------------------------------------------------------
def fig_CSD(batchLabel, simLabel, timeRange):
    from netpyne import sim
    
    dataFolder = '../data/'
    filename = '%s/%s/%s.pkl' % (dataFolder, batchLabel, simLabel) 

    layer_bounds= {'L1': 100, 'L2': 160, 'L3': 950, 'L4': 1250, 'L5A': 1334, 'L5B': 1550, 'L6': 2000}

    fontsize = 16

    import matplotlib; matplotlib.rcParams.update({'font.size': fontsize})


    sim.load(filename, instantiate=False)
    
    figs, axs = sim.plotting.plotCSD(**{
        'spacing_um': 100, 
        'layer_lines': 1, 
        'layer_bounds': layer_bounds, 
        'overlay': 'LFP',
        'timeRange': timeRange, 
        'smooth': 60,
        'fontSize': fontsize,
        'saveFig': '', #filename[:-4]+'_CSD_LFP_%d-%d' % (timeRange[0], timeRange[1]), 
        'figSize': (9, 13), 
        'dpi': 300, 
        'showFig': 0})

    axs[0].set_xlabel('Time (s)', fontsize=fontsize, fontname='Arial')
    axs[0].set_xticks([1000, 1500, 2000])#, ['1.0', '1.5', '2.0'])#, fontsize=fontsize, fontname='Arial')
    axs[0].set_xticklabels(['1.0', '1.5', '2.0'])#, fontsize=fontsize, fontname='Arial')

    plt.margins(x=0,y=0)

    plt.subplots_adjust(bottom=0.10, top=0.95, right=0.85, left=0.13)

    #save_filename=filenames[0][:-4]+'_csd.png'
    save_filename='figs/fig3_LFP_CSD.png'

    plt.savefig(save_filename, dpi=300)


# ---------------------------------------------------------------------------------------------------
def fig_dipole(batchLabel, simLabel, timeRange):

    fontsize = 16

    from netpyne import sim
    dataFolder = '../data/'
    filename = '%s/%s/%s.pkl' % (dataFolder, batchLabel, simLabel) 
    sim.load(filename, instantiate=False)

    # remove artifact due to dipole calculation floating point error
    sim.allSimData['dipoleSum'][1203,:] = sim.allSimData['dipoleSum'][1204,:]  
    
    import matplotlib; matplotlib.rcParams.update({'font.size': fontsize})

    sim.analysis.plotDipole(saveFig=False, showFig=False, timeRange=timeRange, figSize=(9,3))

    plt.xticks([1000, 1500, 2000], ['1.0', '1.5', '2.0'])
    plt.xlim(1000,2000)
    plt.xlabel('Time (s)', fontsize=fontsize)

    plt.title('')
    plt.margins(x=0,y=0)

    plt.subplots_adjust(bottom=0.15, top=0.95, right=0.85, left=0.13)

    # filename='../data/v35_batch9/v35_batch9_0_0_dipole.png'
    filename='figs/fig3_dipole.png'
    plt.savefig(filename, dpi=300)


# ---------------------------------------------------------------------------------------------------
def plotEEG(sim, showCell=None, showPop=None,  timeRange=None, dipole_location='parietal_lobe', dpi=300, fontsize=16, figSize=(19,10), showFig=True, saveFig=True):

    from lfpykit.eegmegcalc import NYHeadModel

    nyhead = NYHeadModel()

    #dipole_location = 'parietal_lobe'  # predefined location from NYHead class
    nyhead.set_dipole_pos(dipole_location)
    M = nyhead.get_transformation_matrix()

    if showCell:
        p = sim.allSimData['dipoleCells'][showCell]    
    elif showPop:
        p = sim.allSimData['dipolePops'][showPop]    
    else:
        p = sim.allSimData['dipoleSum']

    if timeRange is None:
        timeRange = [0, sim.cfg.duration]

    timeSteps = [int(timeRange[0]/sim.cfg.recordStep), int(timeRange[1]/sim.cfg.recordStep)]

    p = np.array(p).T
    
    p=p[:, timeSteps[0]:timeSteps[1]]

    # We rotate current dipole moment to be oriented along the normal vector of cortex
    p = nyhead.rotate_dipole_to_surface_normal(p)
    eeg = M @ p * 1e6  # [mV] -> [nV] unit conversion

    # plot EEG daa
    x_lim = [-100, 100]
    y_lim = [-130, 100]
    z_lim = [-160, 120]

    t = np.arange(timeRange[0], timeRange[1], sim.cfg.recordStep)

    plt.close("all")
    fig = plt.figure(figsize=[9, 14])
    fig.subplots_adjust(top=0.96, bottom=0.05, hspace=0.17, wspace=0.3, left=0.1, right=0.99)
    ax_eeg = fig.add_subplot(313, xlabel="Time (s)", ylabel='nV', title='EEG at all electrodes')
    ax3 = fig.add_subplot(312)
    ax7 = fig.add_subplot(311, aspect=1, xlabel="x (mm)", ylabel='y (mm)', xlim=y_lim, ylim=x_lim)

    dist, closest_elec_idx = nyhead.find_closest_electrode()
    print("Closest electrode to dipole: {:1.2f} mm".format(dist))

    max_elec_idx = np.argmax(np.std(eeg, axis=1))
    time_idx = np.argmax(np.abs(eeg[max_elec_idx]))
    max_eeg = np.max(np.abs(eeg[:, time_idx]))
    max_eeg_idx = np.argmax(np.abs(eeg[:, time_idx]))
    max_eeg_pos = nyhead.elecs[:3, max_eeg_idx]

    for idx in range(eeg.shape[0]):
        ax_eeg.plot(t, eeg[idx, :], c='gray') 
    ax_eeg.plot(t, eeg[closest_elec_idx, :], c='green', lw=2)
    plt.setp(ax_eeg, xlim=[1000,2000], xticks=[1000, 1500, 2000], xticklabels=['1.0', '1.5', '2.0'])

    vmax = np.max(np.abs(eeg[:, time_idx]))
    v_range = vmax
    cmap = lambda v: plt.cm.bwr((v + vmax) / (2*vmax))
    threshold = 2

    xy_plane_idxs = np.where(np.abs(nyhead.cortex[2, :] - nyhead.dipole_pos[2]) < threshold)[0]

    for idx in range(eeg.shape[0]):
        c = cmap(eeg[idx, time_idx])
        ax7.plot(nyhead.elecs[1, idx], nyhead.elecs[0, idx], 'o', ms=10, c=c, 
                zorder=nyhead.elecs[2, idx])

    img = ax3.imshow([[], []], origin="lower", vmin=-vmax,
                    vmax=vmax, cmap=plt.cm.bwr)
    cbar=plt.colorbar(img, ax=ax7, shrink=0.5)
    cbar.ax.set_ylabel('nV', rotation=90)

    ax7.plot(nyhead.dipole_pos[1], nyhead.dipole_pos[0], '*', ms=15, color='orange', zorder=1000)

    plt.subplots_adjust(right=0.9)

    # save figure
    if saveFig:
        if isinstance(saveFig, basestring):
            filename = saveFig
        else:
            filename = sim.cfg.filename + '_EEG.png'
        try:
            plt.savefig(filename, dpi=dpi)
        except:
            plt.savefig('EEG_fig.png', dpi=dpi)


    # display figure
    if showFig is True:
        plt.show()


# ---------------------------------------------------------------------------------------------------
def fig_eeg(batchLabel, simLabel, timeRange):


    from netpyne import sim
    dataFolder = '../data/'
    filename = '%s/%s/%s.pkl' % (dataFolder, batchLabel, simLabel)
    sim.load(filename, instantiate=False)

    # remove artifact due to dipole calculation floating point error
    sim.allSimData['dipoleSum'][1203,:] = sim.allSimData['dipoleSum'][1204,:]  

    import matplotlib; matplotlib.rcParams.update({'font.size': 18})

    # plot EEG fig
    plotEEG(sim, showFig=False, saveFig='figs/fig3_EEG.png', timeRange=timeRange, figSize=(9,4.5))


# ---------------------------------------------------------------------------------------------------
def fig_CSD_comparison():

    '''
    Comparison of spontaneous, BBN and speech CSDs in experiment and model
    '''


    from netpyne import sim

    #-------------------------
    # Options to select experiment vs model and spont vs BBN vs speech input
    exp = 1
    model = 1

    spont = 1
    bbn = 0
    speech = 0


    # whether to plot the spike raster corresponding to the CSD 
    plotRasters = 0 

    #-------------------------
    # Experiment
    if exp:

        from expDataAnalysis import getflayers, loadfile, getbandpass

        # spont
        if spont:
            # spont
            datafile = '../data/NHPdata/spont/2-rb031032016@os_eye06_20.mat'   # LOCAL DIR 
            datafile_pkl = None

            dbpath = '../data/NHPdata/spont/21feb02_A1_spont_layers.csv'  # GCP # CHANGE ACCORDING TO MACHINE USED TO RUN

        # bbn
        elif bbn:
            datafile = '../data/NHPdata/BBN/1-rb067068029@os.mat'  # bbn 200 ms
            
            datafile_pkl = None
            dbpath = '../data/NHPdata/BBN/19jun21_A1_ERP_Layers.csv'   
            
        # speech
        elif speech:
            datafile = '../data/NHPdata/speech/2-bu037038046@os.mat'  #2-ma033034023@os.mat'# 2-bu037038046@os.mat'  
            datafile_pkl = '../data/NHPdata/speech/2-bu037038046@os_small.pkl'
            dbpath = '../data/NHPdata/speech/A1_speech_layers.csv'   

    
        vaknin = 0    

        # setup netpyne
        samprate = 11*1e3  # in Hz
        spacing = 100
        sim.initialize()
        sim.cfg.recordStep = 1000./samprate # in ms

        # load reduced file size
        if datafile_pkl:
            with open(datafile_pkl, 'rb') as f:
                dataLoad = pickle.load(f)

            [sampr, LFP_data, dt, tt, CSD_data, trigtimes] = dataLoad['data_small']

            if datafile_pkl == '../data/NHPdata/speech/2-bu037038046@os_small.pkl':
                LFP_data = getbandpass(LFP_data.T, sampr=samprate, minf=0.5, maxf=300)
        
        else:
            # load and reduce size by 4        
            [sampr,LFP_data,dt,tt,CSD_data,trigtimes] = loadfile(fn=datafile, samprds=samprate, spacing_um=100, vaknin=vaknin)
            
            save_reduced = False
            if save_reduced:
                reduced_size = int(len(tt) / 4.)
                LFP_data_small=LFP_data[:,1:reduced_size]
                CSD_data_small=CSD_data[:,1:reduced_size]
                tt_small=tt[1:reduced_size]
                trigtimes=[x for x in trigtimes if x<4977959]
                dataSave = {'data_small': [sampr,LFP_data_small,dt,tt_small,CSD_data_small,trigtimes]}
                
                with open ('../data/NHPdata/speech/2-ma033034023@os_small.pkl', 'wb') as f:
                    pickle.dump(dataSave,f)
        

        ##### GET LAYERS FOR OVERLAY #####
        s1low,s1high,s2low,s2high,glow,ghigh,i1low,i1high,i2low,i2high = getflayers(datafile,dbpath=dbpath,getmid=False,abbrev=False) # fullPath is to data file, dbpath is to .csv layers file 
        lchan = {}
        lchan['S'] = s2high
        lchan['G'] = ghigh
        lchan['I'] = CSD_data.shape[0]-1 #i2high
        print('s2high: ' + str(s2high))

        if bbn:
            depth = 2000
        else:
            depth = LFP_data.shape[0]*100 

        layer_bounds= {'S': (glow*spacing)+spacing, 'G': (i1low*spacing)+spacing, 'I': depth} 

        # plot CSD
        smooth = 30

        if spont:
            tranges = [[57400, 57600], [64100, 64300], [66400, 66600]]

        elif bbn:
            tranges = [[46705, 46705+200]]

        elif speech:
            tranges = [[41454, 41654], [93734, 93934]]
            
        num_elec = 21 # for speech

        for i,trange in enumerate(tranges): #speech 
            print(i)

            figname = 'figs/fig4_'+datafile.split('/')[-1][:-4]+'_CSD_LFP_smooth-%d_%d-%d_vaknin-%d_colorbar' % (smooth, trange[0], trange[1], vaknin)

            sim.plotting.plotCSD(**{
                'CSDData': CSD_data[:num_elec], 
                'LFPData': LFP_data[:num_elec].T,
                'vmin': -1.1, # spont: -1.1, BBN: -1.7; speech: -1.8
                'vmax': 1.1, #  spont: 1.1, BBN: 1.5; speech: 1.8
                'spacing_um': 100, 
                'dt': sim.cfg.recordStep,
                'ymax': depth,
                'layerBounds': layer_bounds, 
                'legendLabel': False,
                'overlay': 'LFP',
                'timeRange': [trange[0], trange[1]], 
                'smooth': smooth,
                'vaknin': vaknin,
                'saveFig': figname, 
                'figSize': (4.1,8.2), 
                'dpi': 300, 
                'colorbar': True,
                'showFig': 0})
            

    # ------------------------------------
    # Model
    if model:
        
        #layer_bounds= {'L1': 100, 'L2': 160, 'L3': 950, 'L4': 1250, 'L5A': 1334, 'L5B': 1550, 'L6': 2000}
        layer_bounds= {'S': 950, 'G': 1250, 'I': 2000}#, 'L6': 2000}
        
        if spont:
            filename = '../data/v34_batch56/v34_batch56_0_0_data.pkl' # spont 5*5 seeds
            tranges = [[3400, 3600], [2800, 3000], [5700, 5900]]

        elif bbn:
            filename = '../data/BBN_CINECA_v36_5656BF_SOA850/BBN_CINECA_v36_5656BF_850ms_data.pkl' # bbn
            tranges = [[4200, 4400]] 

        elif speech:
            filenames = ['../data/v34_batch58/v34_batch58_0_1_0_1_3_data.pkl',
                        '../data/v34_batch58/v34_batch58_0_0_1_1_3_data.pkl'] # speech
            tranges = [[3100, 3300], [3400, 3600]]


        vaknin = 0
        smooth = 30
        if not plotRasters:
            for filename,trange in zip(filenames,tranges):
                sim.load(filename, instantiate=False)
        
                figname = 'figs/fig4_'+filename.split('/')[-1][:-4]+'_CSD_LFP_smooth-%d_%d-%d_vaknin-%d_colorbar' % (smooth, trange[0], trange[1], vaknin)

                sim.plotting.plotCSD(**{
                    'vmin':-12, # spont: -12; BBN: -33; speech: -45
                    'vmax': 12, # spont: 12; BBN: 22; speech: 45
                    'spacing_um': 100, 
                    'dt': sim.cfg.recordStep,
                    'layerBounds': layer_bounds, 
                    'legendLabel': False,
                    'overlay': 'LFP',
                    'timeRange': [trange[0], trange[1]], 
                    'smooth': smooth,
                    'vaknin': vaknin,
                    'saveFig': figname, 
                    'figSize': (4.1,8.2), 
                    'dpi': 300, 
                    'colorbar': True,
                    'showFig': 0})
            
        if plotRasters:
            fontsize = 20
            for timeRange in tranges:           

                figname = 'figs/fig4_'+filename.split('/')[-1][:-4]+'_raster_%d_%d.png'%(timeRange[0], timeRange[1])

                sim.plotting.plotRaster(
                    **{ 'include': ['allCells'], 
                        'saveFig': None, 
                        'showFig': False, 
                        'popRates': 'minimal', 
                        'orderInverse': True, 
                        'timeRange': timeRange, 
                        'figSize': (3,8), 
                        'lw': 0.3, 
                        'markerSize': 3, 
                        'marker': '.', 
                        'dpi': 300})

                ax = plt.gca()

                [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

                ax.get_legend().remove()

                plt.title('')
                plt.margins(x=0,y=0)

                plt.subplots_adjust(bottom=0.10, top=0.98, right=0.85, left=0.15)

                plt.savefig(figname, dpi=300)            
   

# ---------------------------------------------------------------------------------------------------
def fig_BBN():

    filenameList = ['../data/BBN_CINECA_v36_5656BF_SOA850/BBN_CINECA_v36_5656BF_850ms_data.pkl'] # bbn

    plotRaster = 0
    plotSpikeHist = 1
    calculateAdaptation = 1

    stimDur = 100  # BBN stim duration (100 ms)

    meanRatesAll = {}
    peakRatesAll = {}

    for ifile, filename in enumerate(filenameList):

        sim.load(filename, instantiate=False)

        stimTimes = sim.cfg.ICThalInput['startTime']  # get BBN stim start times

        # --------------------------------------------
        # Raster
        if plotRaster:
            timeRange = [3000, 5000] 
            # plot full raster
            sim.plotting.plotRaster(
                **{ 'include': ['allCells'], 
                    'saveFig': filename[:-4]+'_raster_%d_%d.png'%(timeRange[0], timeRange[1]), 
                    'showFig': False, 
                    'popRates': 'minimal', 
                    'orderInverse': True, 
                    'timeRange': timeRange, 
                    'figSize': (18*1.5,10*1.5), 
                    'lw': 0.1, 
                    'markerSize': 0.1, 
                    'marker': '.', 
                    'dpi': 300})

            ax = plt.gca()

            [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

            ax.get_legend().remove()

            plt.title('')
            plt.margins(x=0,y=0)

            plt.subplots_adjust(bottom=0.10, top=0.98, right=0.95, left=0.10)

            save_filename='figs/fig4_%s_raster_%d_%d.png'%(filename.split('/')[-1][:-4], timeRange[0], timeRange[1])
            plt.savefig(save_filename, dpi=300)  


        # --------------------------------------------
        # Spike Hist
        if plotSpikeHist:
            from netpyne.analysis.spikes_legacy import plotSpikeHist

            includeList = [['IC', thal, cortexE]] 
            includeLabels = ['IC_thal_cxE'] 

            timeRange = [3000, 5000]   # for plot
            if calculateAdaptation: 
                timeRange = [2000, 11500]  # for adaptation measure 

            for include, includeLabel in zip(includeList, includeLabels):

                fontsize = 14
                binSize = 5
                measure = 'count'
                graphType = 'bar' #'line'                

                plt.set_cmap('tab10')
                ''''
                tab10:
                #1f77b4 - blue
                #ff7f0e - orange
                #2ca02c - green
                #d62728 - red
                #9467bd - purple
                #8c564b
                #e377c2
                #7f7f7f
                #bcbd22
                #17becf

                '''

                fig_hist, data_hist = plotSpikeHist(include=include, 
                                        popColors=['#02216b', '#ff7f0e', '#2ca02c'], #{'IC': 'blue', tuple(thal): 'green', tuple(cortexE): 'red'}, 
                                        lineWidth=2, 
                                        lw=2,
                                        fontSize=fontsize,
                                        timeRange=timeRange, 
                                        dpi=300, 
                                        figSize=(16,5), 
                                        legend=False, 
                                        showFig=0, 
                                        saveFig=0, 
                                        binSize=binSize, 
                                        graphType=graphType, 
                                        axis=True, 
                                        measure=measure, 
                                        scalebarLoc='upper center')

                ax = plt.gca()

                for line in ax.lines:
                    line.set_linewidth(1.5)

                [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

                ylim = ax.get_ylim()

                plt.ylim(0,115)
                plt.ylabel('Spike count', fontsize=fontsize, fontname='Arial') #Neurons (ordered by NCD within each pop)')
                plt.xlabel('Time (ms)', fontsize=fontsize, fontname='Arial')
                
                leg = ax.get_legend()
                leg.set_title('')
                
                leg.set_bbox_to_anchor((1.02, 1.0))

                for text in leg.get_texts():
                    text.set_fontsize(fontsize)
                    text.set_fontname('Arial')
                
                plt.title('')
                plt.margins(x=0,y=0)

                plt.subplots_adjust(bottom=0.1, top=0.98, right=0.90 , left=0.05)

                if not calculateAdaptation:
                    save_filename = 'figs/fig4_%s_spikehist_%d_%d_%s.png'%(filename.split('/')[-1][:-4], timeRange[0], timeRange[1], includeLabel)
                    print(save_filename)
                    plt.savefig(save_filename, dpi=300)

        
                # --------------------------------------------
                # Calculate adaptation
                if calculateAdaptation:
                    import scipy

                    popList = data_hist['include']
                    histoTime = data_hist['histoT']

                    # adaptation inter-stimuli
                    plt.figure(figsize=(8,6))
                    
                    for ipop, pop in enumerate(popList):
                        print(pop)

                        histoCount = data_hist['histoData'][ipop]

                        meanRatesInter = []
                        peakRatesInter = []

                        for t in stimTimes:
                            print('  ', t)
                            startIndex = np.where(histoTime == round(t,-1) + binSize/2.0)[0][0] # times based on bin center (eg 5ms for 0-10ms bin)
                            endIndex = int(startIndex+(150/binSize)) #look at 150 ms window even though stim is 100ms
                            meanRatesInter.append(np.mean(histoCount[startIndex:endIndex]))
                            peakRatesInter.append(max(histoCount[startIndex:endIndex]))

                        # plot
                        plt.plot(meanRatesInter, label='%s mean'%(pop))
                        plt.plot(peakRatesInter, label='%s peak'%(pop))
                    
                    plt.legend()
                    save_filename = '%s_adapatation_inter_%d_%d_%s.png'%(filename[:-4], timeRange[0], timeRange[1], includeLabel)
                    print(save_filename)
                    #plt.savefig(save_filename, dpi=300)

                    # adaptation inter-stimuli
                    plt.figure(figsize=(8,6))
                    
                    for ipop, pop in enumerate(popList):
                        print(pop)

                        histoCount = data_hist['histoData'][ipop]
                        
                        meanRatesIntra = []
                        peakRatesIntra = []

                        peakPres = []
                        peakPosts = []

                        avgPres = []
                        avgPosts = []

                        for t in stimTimes:
                            print('  ', t)
                            startIndex = np.where(histoTime == round(t, -1) + binSize/2.0)[0][0] # times based on bin center (eg 2.5ms for 0-5ms bin)
                            endIndex = int(startIndex+(100/binSize))  # 5 ms bins (look at 100 ms window)
                            meanRatesIntra.append(np.mean(histoCount[startIndex:startIndex+2]) - np.mean(histoCount[endIndex-2:endIndex]))
                            peakStart = max(histoCount[startIndex:startIndex+int(20/binSize)])
                            peakEnd = max(histoCount[endIndex-int(20/binSize):endIndex])

                            #print('%s %s: peakStart: %.2f, peakEnd: %.2f' % (pop, measure, peakStart, peakEnd))
                            peakRatesIntra.append(peakStart - peakEnd)

                            peakPost = max(histoCount[startIndex:startIndex+int(200/binSize)])
                            peakPre = max(histoCount[startIndex-int(200/binSize):startIndex-1])
                            peakPosts.append(peakPost)
                            peakPres.append(peakPre)

                            avgPost = np.mean(histoCount[startIndex:startIndex+int(200/binSize)])
                            avgPre = np.mean(histoCount[startIndex-int(200/binSize):startIndex-1])
                            avgPosts.append(avgPost)
                            avgPres.append(avgPre)

                        # plot
                        plt.plot(meanRatesIntra, label='%s mean'%(pop))
                        plt.plot(peakRatesIntra, label='%s peak'%(pop))

                        # For Thal an cx use 200ms peak response post - pre stim   
                        peakPrePostDiff = np.array(peakPosts)-np.array(peakPres)                      
                        
                        print('%s %s: PEAK pre: mean=%.2f, std=%.2f; PEAK post: mean=%.2f, std=%.2f; PEAK post-pre: mean=%.2f, std=%.2f;' 
                              % (pop, measure, 
                                 np.mean(peakPres), np.std(peakPres),
                                 np.mean(peakPosts), np.std(peakPosts),
                                 np.mean(peakPrePostDiff), np.std(peakPrePostDiff)))

                        print(scipy.stats.mannwhitneyu(peakPres, peakPosts))


                        avgPrePostDiff = np.array(peakPosts)-np.array(peakPres)   
                        print('%s %s: AVG pre: mean=%.2f, std=%.2f; AVG post: mean=%.2f, std=%.2f; AVG post-pre: mean=%.2f, std=%.2f;' 
                              % (pop, measure, 
                                 np.mean(avgPres), np.std(avgPres),
                                 np.mean(avgPosts), np.std(avgPosts),
                                 np.mean(avgPrePostDiff), np.std(avgPrePostDiff)))

                        print(scipy.stats.mannwhitneyu(avgPres, avgPosts))
                        
                    plt.legend()
                    save_filename = '%s_adapatation_intra_%d_%d_%s_%s.png'%(filename[:-4], timeRange[0], timeRange[1], includeLabel, measure)
                    print(save_filename)
                    #plt.savefig(save_filename, dpi=300)

# ---------------------------------------------------------------------------------------------------
def fig_BBN_ERP(usemodel=0):
    from avgERP import drawDataAvgERP

    # usemodel = 0 # 0=NHP, 1=model

    if usemodel:
        ### SIM ###	
        regions = [6,11,12] 

        basedirSim = '../data/'
        dataFiles = ['BBN_CINECA_v36_5656BF_SOA850/BBN_CINECA_v36_5656BF_850ms_data.pkl'] 
        
        for region in regions:
            for filePath in dataFiles:
                dataFileFullPath = basedirSim + filePath
                drawDataAvgERP(usemodel=usemodel, region=region, dataFile=dataFileFullPath, saveDir='figs/fig4_', lfpSmooth=1, swindowms=0, ewindowms=300, ERPdata='CSD') 


    else:
        ## NHP ###
        regions = [4,10,13] 

        NHPdataFiles = {'1-rb067068029@os.mat': {'stimType': 'BBN'}}

        basedirNHP = '../data/NHPdata/BBN/'
        dataFiles = ['1-rb067068029@os.mat']


        for region in regions:
            for dataFile in dataFiles:
                dataFileFullPath = basedirNHP + dataFile
                drawDataAvgERP(usemodel=usemodel, region=region, dataFile=dataFileFullPath, saveDir='figs/fig4_', lfpSmooth=1, swindowms=0, ewindowms=300, chan=None, dataFilesInfoDict=NHPdataFiles, ERPdata='CSD')



# ---------------------------------------------------------------------------------------------------
def plot_LFP_PSD_combined(dataFile, norm=False):
    NHP = dataFile[:-4]
    
    with open(dataFile, 'rb') as f:
        loadedData = pickle.load(f)
        allData = loadedData['allData']
        freqs = loadedData['allFreqs'][0]

    plt.figure(figsize=(8,12))
    fontSize = 12
    lw = 1

    elecLabels = ['All layers', 'Supragranular', 'Granular', 'Infragranular']

    for i in range(len(elecLabels)):
        plt.subplot(4, 1, i+1)
        for itime in range(len(allData)):
            signal = allData[itime][i]
            if norm:
                signal = signal / max(signal)
            plt.plot(freqs, signal, linewidth=lw) #, color=color)
        
        plt.title(elecLabels[i])
        plt.xlim([0, 50])
        if norm:
            plt.ylabel('Normalized Power', fontsize=fontSize)
        else:
            plt.ylabel('Power (mV^2/Hz)', fontsize=fontSize)
        plt.xlabel('Frequency (Hz)', fontsize=fontSize)
    
    # format plot 
    plt.tight_layout()
    plt.suptitle('LFP PSD - %s' % (NHP), fontsize=fontSize, fontweight='bold') # add yaxis in opposite side
    plt.subplots_adjust(bottom=0.08, top=0.92)

    if norm:
        plt.savefig('figs/fig5_%s_PSD_combined_norm.png' % (NHP.split('/')[-1][:-4]))
    else:
        plt.savefig('figs/fig5_%s_PSD_combined.png' % (NHP.split('/')[-1][:-4]))


# ---------------------------------------------------------------------------------------------------
def fig_LFP_PSD_comparison():
    import seaborn as sns
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler as Sc
    from scipy.spatial.distance import cdist, pdist
    import random

    dataFiles = ['../data/NHPdata/spont/spont_LFP_PSD/1-bu031032017@os_eye06_20_10sec_allData.pkl', 
                '../data/NHPdata/spont/spont_LFP_PSD/2-ma031032023@os_eye06_20_10sec_allData.pkl', 
                '../data/NHPdata/spont/spont_LFP_PSD/2-rb031032016@os_eye06_20_10sec_allData.pkl', 
                '../data/NHPdata/spont/spont_LFP_PSD/2-rb045046026@os_eye06_20_10sec_allData.pkl',
                '../data/v34_batch57/v34_batch57_10sec_allData.pkl']

    plot_combined_PSD = 1
    plotPCA = 1
    plotMatrix = 1

    psd = []
    clusters = []
    start = 0
    elec = 0
    normalize = 1

    for dataFile in dataFiles:            
        NHP = dataFile[:-4]
    
        file_psds = []

        with open(dataFile, 'rb') as f:
            loadedData = pickle.load(f)
            allData = loadedData['allData']
            freqs = loadedData['allFreqs'][0]

        for itime in range(len(allData)):
            signal = allData[itime][elec]
            if normalize:
                signal = signal / max(signal)
            file_psds.append(signal)

        psd.extend(file_psds)

        clusters.append(range(start, start+len(allData)))
        start = start + len(allData)

        if plot_combined_PSD:
            plot_LFP_PSD_combined(dataFile, norm=normalize)


    labels = ['Macaque 1', 'Macaque 2', 'Macaque 3 (Day 1)', 'Macaque 3 (Day 2)', 'Model', 'Model shuffled']
    df = pd.DataFrame(psd)
    
    # add model shuffled
    model_shuffled = []
    random.seed(6)
    for row in df.iloc[-25:].iterrows():
        rowshuff = list(row[1]) 
        random.shuffle(rowshuff)
        model_shuffled.append(rowshuff)
    clusters.append(range(start, start+len(model_shuffled)))
    
    model_shuffled = pd.DataFrame(model_shuffled)
    df=pd.concat([df, model_shuffled])#.reset_index()
    

    # PCA
    if plotPCA:

        dfnorm = df.div(df.max(axis=1), axis=0)   
        
        pca = PCA(n_components=2)
        X_pca = pca.fit_transform(dfnorm)

        print('Explained variance ratio: ', pca.explained_variance_ratio_)

        fig = plt.figure(figsize=(7, 5))
        colors = ['blue', 'r', 'm', 'y',  'limegreen', 'k']

        for cluster, label, c in zip(clusters, labels, colors): # 'rgbcmykw'):
            plt.scatter(X_pca[cluster, 0], X_pca[cluster, 1], c=c, label=label)

        cluster_points = {}
        for cluster, label, c in zip(clusters, labels, colors): # 'rgbcmykw'):
            cluster_points[label] = np.array(list(zip(X_pca[cluster, 0], X_pca[cluster, 1])))

        # distance between Macaque 2 and Macaque 3 vs Macaque2 vs model
        print('Means:')
        print('Macaque 2 vs Macaque 3 (Day 1)', np.mean(cdist(cluster_points['Macaque 2'], cluster_points['Macaque 3 (Day 1)'])))
        print('Macaque 3 Day 1 vs Day 2', np.mean(cdist(cluster_points['Macaque 3 (Day 1)'], cluster_points['Macaque 3 (Day 2)'])))
        print('Macaque 2 vs Model', np.mean(cdist(cluster_points['Macaque 2'], cluster_points['Model'])))

        pca_mean_distances = np.zeros((len(labels), len(labels)))
        pca_std_distances = np.zeros((len(labels), len(labels)))

        for i,ilabel in enumerate(labels):
            for j,jlabel in enumerate(labels):
                pca_mean_distances[i,j] = np.mean(cdist(cluster_points[ilabel], cluster_points[jlabel]))
                pca_std_distances[i,j] = np.std(cdist(cluster_points[ilabel], cluster_points[jlabel]))

        within_exp = [pca_mean_distances[0,0], pca_mean_distances[1,1], pca_mean_distances[2,2], pca_mean_distances[3,3]] 
        print(np.mean(within_exp), np.std(within_exp))

        across_exp = [pca_mean_distances[x[0], x[1]] for x in
                    [(0,1), (0,2), (0,3), (0,4),(1,2), (1,3), (1,4)]]
        print(np.mean(across_exp), np.std(across_exp))

        nhps = ['Macaque 1', 'Macaque 2', 'Macaque 3 (Day 1)', 'Macaque 3 (Day 2)']

        # distance between macaques and model vs macaques and shuffled
        distToModel = 0
        distToModelShuffled = 0
        
        for nhp in nhps:
            distToModel += np.mean(cdist(cluster_points[nhp], cluster_points['Model']))
            distToModelShuffled += np.mean(cdist(cluster_points[nhp], cluster_points['Model shuffled']))

        distToModel /= 4
        distToModelShuffled /= 4

        print('distToModel: ', distToModel)
        print('distToModelShuffled: ', distToModelShuffled)

        plt.xlabel('PC 1')
        plt.ylabel('PC 2')
        plt.subplots_adjust(left=0.1,right=0.7, top=0.9, bottom=0.1)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig('figs/fig5_PSD_norm_PCA_%d.png'%(normalize))



    if plotMatrix:
        # corr matrix
        dfT = df.T
        corrMatrix = dfT.corr()
        plt.figure(figsize=(10,10))
        im = sns.heatmap(corrMatrix, cmap='viridis', square=True)    
        ax = plt.gca()
        ax.hlines([34, 63, 101, 123, 148], *ax.get_xlim(), colors=['white']*5, linestyle='dotted', linewidth=1.5)
        ax.vlines([34, 63, 101, 123, 148], *ax.get_ylim(), colors=['white']*5, linestyle='dotted', linewidth=1.5)
        plt.savefig('figs/fig5_PSD_corr_matrix_%d.png'%(normalize))

        # mean correlation
        meanCorr = {}
        for cluster, label in zip(clusters, labels):
            meanCorr[label] = corrMatrix.iloc[cluster, range(clusters[0][0], clusters[3][-1])].mean().mean()

        # stats
        import scipy
        modelCorr=corrMatrix.iloc[clusters[-2], range(clusters[0][0], clusters[3][-1])]
        modelShuffCorr=corrMatrix.iloc[clusters[-1], range(clusters[0][0], clusters[3][-1])]    
        ranksum1 = scipy.stats.mannwhitneyu(modelCorr, modelShuffCorr)    

        print(meanCorr)
        print('ranksum model-exp vs shuff-exp: ', ranksum1)

        print(modelCorr.values.mean(), modelCorr.values.std())
        print(modelShuffCorr.values.mean(), modelShuffCorr.values.std())

        modelvsShuffCorr = corrMatrix.iloc[clusters[-2], clusters[-1]]
        # ranksum2 = scipy.stats.mannwhitneyu(modelvsShuffCorr)    
        # print('ranksum model vs shuff: ', ranksum2)

        print(modelvsShuffCorr.values.mean(), modelvsShuffCorr.values.std())



# ---------------------------------------------------------------------------------------------------
def fig_osc_events():
    from osc_events import plotOscEvents

    ### Dict with paths to data directories ### 
    dataPaths = {'sim': '../data/v34_batch57/osc_events_data/', 'nhp': '../data/NHPdata/spont/osc_events_data/'} 

    ### Dict with oscillation events info ### 
    oscEventsInfo = {'gamma': 
                        {'sim': {'subjName': 'v34_batch57_4_4_data_timeRange_0_6', 'chan': 12, 'eventIdx': 1423, 'specrange': (0,30), 'ylspec': (30,95)},  
                        'nhp':{'subjName': '2-bu027028013_timeRange_0_40', 'chan': 14, 'eventIdx': 2658, 'specrange': (0,30), 'ylspec': (30,95)}}, 
                    'beta': 
                        {'sim': {'subjName': 'v34_batch57_3_2_data_timeRange_0_6', 'chan': 14, 'eventIdx': 1710, 'specrange': (0,30), 'ylspec': (10,50)}, 
                        'nhp': {'subjName': '2-rb031032016_timeRange_40_80', 'chan': 14, 'eventIdx': 2342, 'specrange': (0,25), 'ylspec': (10,50)}}, 
                    'alpha': 
                        {'sim': {'subjName': 'v34_batch57_3_2_data_timeRange_6_11', 'chan': 9, 'eventIdx': 863, 'specrange': (0,20), 'ylspec': (1,30)}, 
                        'nhp':{'subjName': '2-bu027028013_timeRange_80_120', 'chan': 7, 'eventIdx': 784, 'specrange': (0,20), 'ylspec': (1,30)}}, 
                    'theta': 
                        {'sim': {'subjName': 'v34_batch57_3_3_data_timeRange_0_6', 'chan': 8, 'eventIdx': 973, 'specrange': (0,20), 'ylspec': (1,12)}, 
                        'nhp':{'subjName': '2-rb031032016_timeRange_40_80', 'chan': 6, 'eventIdx': 747, 'specrange': (0,15), 'ylspec': (1,12)}}, 
                    'delta': 
                        {'sim': {'subjName': 'v34_batch57_3_4_data_timeRange_0_6', 'chan': 14, 'eventIdx': 1666, 'specrange': (0,30), 'ylspec': (1,10)}, 
                        'nhp':{'subjName': '2-rb031032016_timeRange_160_200', 'chan': 18, 'eventIdx': 3020, 'specrange': (0,30), 'ylspec': (1,10)}}}
    
    plotOscEvents(oscEventsInfo, dataPaths, ['gamma', 'beta', 'alpha', 'theta', 'delta'], eventTypes=['nhp'], saveFig=1)  


# ---------------------------------------------------------------------------------------------------
def fig_beta_osc_event():
    from osc_events import plotOscEvents

    dataPaths_beta = {'sim': '../data/v34_batch67/', 'nhp': '../data/NHP_data/spont/osc_events_data/'}

    ### Dict with beta  ###
    betaOscEventInfo = {'beta': 
                            {'sim': {'subjName': 'v34_batch67_CINECA_0_0_data', 'chan': 19, 'eventIdx': 4336, 'specrange': (0,20), 'ylspec': (2,50)},
                            'nhp': {'subjName': '2-bu027028013_timeRange_40_80', 'chan': 14, 'eventIdx': 2241, 'specrange': (0,20), 'ylspec': (2,50)}}}


    plotOscEvents(betaOscEventInfo, dataPaths_beta, ['beta'], eventTypes=['sim'], saveFig=1)


# ---------------------------------------------------------------------------------------------------
def fig_osc_stats():
    from osc_stats import plotStats
    plotStats()

# ---------------------------------------------------------------------------------------------------
def fig_CSD_osc_spiking():

    from netpyne import sim

    timeRange = [2149-50, 2332+50]
    fontsize = 20

    plotRaster = 1
    plotSpikeHist = 1 


    filename = '../data/v34_batch67/v34_batch67_CINECA_0_0_data.pkl'
    sim.load(filename, instantiate=False, setupRecording=1)

    if plotRaster:
        sim.plotting.plotRaster(**{'include': ['allCells'], 
                                'saveFig': filename[:-4]+'_raster_%d_%d.png'%(timeRange[0], timeRange[1]), 
                                'showFig': False, 
                                'popRates': 'minimal', 
                                'orderInverse': True, 
                                'timeRange': timeRange, 
                                'figSize': (10,8), 
                                'lw': 0.3, 
                                'markerSize': 3, 
                                'marker': '.', 
                                'dpi': 300})

        ax = plt.gca()

        [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

        # 170 = 0 ms  (total 283) ticks at 0, -50, -100, -150, 50, 10
        center = 2149-50 + 170
        plt.xticks([center-100, center, center+100], ['-100', '0', '100'], fontsize=fontsize, fontname='Arial')
        plt.yticks([0, 2000, 4000, 6000, 8000, 10000, 12000], [0, 2000, 4000, 6000, 8000, 10000, 12000], fontsize=fontsize, fontname='Arial')
        plt.xlim(timeRange[0], timeRange[1])
        plt.ylim(12908, 0)
        plt.ylabel('Neuron ID', fontsize=fontsize, fontname='Arial') #Neurons (ordered by NCD within each pop)')
        plt.xlabel('Time (ms)', fontsize=fontsize, fontname='Arial')
        
        leg = ax.get_legend()
        leg.set_title('')
        
        leg.set_bbox_to_anchor((1.05, 1.0))

        for text in leg.get_texts():
            text.set_fontsize(fontsize)
            text.set_fontname('Arial')
        
        plt.title('')
        plt.margins(x=0,y=0)

        plt.subplots_adjust(bottom=0.10, top=0.98, right=0.85, left=0.15)

        save_filename='figs/fig7_%s_raster_%d_%d.png'%(filename.split('/')[-1][:-4], timeRange[0], timeRange[1])
        plt.savefig(save_filename, dpi=300)


    if plotSpikeHist:
        # --------------------------------------------
        # Spike Hist
        from netpyne.analysis.spikes_legacy import plotSpikeHist,plotRatePSD

        fontsize = 30
        binSize = 10
        measure = 'count'
        graphType = 'line'

        popColors = {'CT6': colorList[1], 'IT6': colorList[0], 'PT5B': colorList[10]}

        include = ['CT6', 'IT6']

        fig_hist = plotSpikeHist(include=include, 
                                popColors=popColors, 
                                lineWidth=3, 
                                lw=3,
                                fontSize=fontsize,
                                timeRange=timeRange, 
                                dpi=300, 
                                figSize=(12,8), 
                                legend=False, 
                                showFig=0, 
                                saveFig=filename[:-4]+'_spikehist_%d_%d.png'%(timeRange[0], timeRange[1]),
                                binSize=binSize, 
                                graphType=graphType, 
                                axis=True, 
                                measure=measure, 
                                scalebarLoc='upper center')

        ax = plt.gca()

        for line in ax.lines:
            line.set_linewidth(4.)

        [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

        # 170 = 0 ms  (total 283) ticks at 0, -50, -100, -150, 50, 10
        center = 2149-50 + 170
        plt.xticks([center-100, center, center+100], ['-100', '0', '100'], fontsize=fontsize, fontname='Arial')
        plt.xlim(timeRange[0], timeRange[1])
        plt.ylabel('Spike count', fontsize=fontsize, fontname='Arial') #Neurons (ordered by NCD within each pop)')
        plt.xlabel('Time (ms)', fontsize=fontsize, fontname='Arial')
        
        leg = ax.get_legend()
        leg.set_title('')
        
        leg.set_bbox_to_anchor((0.83, 0.97))

        for text in leg.get_texts():
            text.set_fontsize(fontsize)
            text.set_fontname('Arial')
        
        plt.title('')
        plt.margins(x=0,y=0)

        plt.subplots_adjust(bottom=0.12, top=0.98, right=0.85, left=0.1)

        save_filename='figs/fig7_%s_spikehist_%d_%d.png'%(filename.split('/')[-1][:-4], timeRange[0], timeRange[1])
        plt.savefig(save_filename, dpi=300)


    # --------------------------------------------
    # Rate PSD
    fig_psd = plotRatePSD(include=['IT6', 'CT6'],  
                          popColors=popColors, 
                          lineWidth=3,
                          fontSize=fontsize, 
                          timeRange=timeRange, 
                          figSize=(12,8), 
                          showFig=0, 
                          saveFig=filename[:-4]+'_psd_%d_%d.png'%(timeRange[0], timeRange[1]), 
                          minFreq=1, 
                          maxFreq=50, 
                          binSize=5)
    
    ax = plt.gca()

    for line in ax.lines:
        line.set_linewidth(4.)

    [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

    plt.ylabel('Power', fontsize=fontsize, fontname='Arial') #Neurons (ordered by NCD within each pop)')
    plt.xlabel('Frequency (Hz)', fontsize=fontsize, fontname='Arial')
    plt.xticks(fontsize=fontsize, fontname='Arial')
    plt.yticks(fontsize=fontsize, fontname='Arial')

    
    leg = ax.get_legend()
    leg.set_title('')
    
    leg.set_bbox_to_anchor((0.7, 0.95))

    for text in leg.get_texts():
        text.set_fontsize(fontsize)
        text.set_fontname('Arial')
    
    plt.title('')
    plt.margins(x=0,y=0)

    plt.subplots_adjust(bottom=0.12, top=0.98, right=0.85, left=0.12)

    save_filename='%s_spikePSD_%d_%d.png'%(filename[:-4], timeRange[0], timeRange[1])
    plt.savefig(save_filename, dpi=300)


# ---------------------------------------------------------------------------------------------------
def fig_CSD_contrib_PSD():
    from csd_analysis import getCSDDataFrames
    fontsiz = 16
    
    dataFile = '../data/v34_batch67/v34_batch67_CINECA_0_0_data.pkl'
    timeRange = [2149.66, 2332.71]
    betaOscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': -2270.0, 'left': 42993, 'right': 46654, 'w2':1098}
    
    ## Get peak and avg dataframes  
    dfCSDPeak, dfCSDAvg = getCSDDataFrames(dataFile=dataFile, timeRange=timeRange, oscEventInfo=betaOscEventInfo) ## going to need oscEventInfo here 


    dfCSDAvg[19].abs().sort_values(ascending=False)[0:13].plot.bar(fontsize = fontsiz)
    plt.subplots_adjust(bottom=0.18, top=0.98, right=0.85, left=0.1)
    plt.savefig('figs/fig7_pop_contrib.png', dpi=300)


# ---------------------------------------------------------------------------------------------------
def fig_fI_curve():

    from batchAnalysisPlotCombined import fIAnalysis

    dataFolder = '../data/'
    batchLabel = 'v22_batch14'  # 'v50_batch1' #
    loadAll = 0

    fIAnalysis(dataFolder, batchLabel, loadAll)


# ---------------------------------------------------------------------------------------------------
def fig_optuna_fitness():
    import optunaAnalysis as oa
    import seaborn as sns

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

            #import IPython; IPython.embed()

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
            y = v #(np.array(v)-min(v))/(max(v)-min(v)) # normalize
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
            y = v #(np.array(v)-min(v))/(max(v)-min(v)) # normalize
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
                #plt.title('%s vs %s R=%.2f' % ('fitness', param.replace('tune', ''), dfcorr['value'][param]))
                plt.savefig('%s/%s/%s_scatter_%s_%s.png' % (dataFolder, batchSim, batchSim, 'fitness', param.replace('tune', '')), dpi=300)



# ---------------------------------------------------------------------------------------------------
def fig_speech():

    plotWav = 1
    plotICrates = 1
    plotRaster = 0
    plotSpikeHist = 1

    if plotWav:
            
        # plot .wav
        from scipy.io.wavfile import read
        from filter import bandpass

        def getbandpass (lfps,sampr,minf=0.05,maxf=300):
            datband = []
            for i in range(len(lfps)): datband.append(bandpass(lfps[:,i],minf,maxf,df=sampr,zerophase=True))
            datband = np.array(datband)
            return datband

        # read audio samples
        plt.figure(figsize=(14,6))
        fs = 100000
        input_data = read("../data/ICoutput/01_ba_peter.wav")
        audio = input_data[1][:,0]
        audio_filt = bandpass(audio, 9600,10400, fs)
        plt.plot(audio, lw=1)
        plt.ylabel("Amplitude")
        plt.xlabel("Time")
        plt.savefig('figs/figS5_01_ba_peter_wav.png', dpi=300)

        plt.figure(figsize=(14,6))
        plt.plot(audio_filt, lw=1)
        plt.ylabel("Amplitude")
        plt.xlabel("Time")
        plt.savefig('figs/figS5_01_ba_peter_wav_filtered.png', dpi=300)

    if plotICrates:
        #plot IC rates
        import scipy.io
        plt.figure(figsize=(14,6))
        mat = scipy.io.loadmat("../data/ICoutput/ICoutput_CF_9600_10400_wav_01_ba_peter.mat")
        audio = mat['BE_sout_population']
        plt.plot(audio.T, lw=1)
        plt.ylabel("Amplitude")
        plt.xlabel("Time")
        plt.savefig('figs/figS5_01_ba_peter_ICoutput.png', dpi=300)



    # --------------------------------------------
    # Raster
    if plotRaster:
        
        filename = '../data/v34_batch58/v34_batch58_0_0_0_0_0_data.pkl' # speech
        sim.load(filename, instantiate=False)

        timeRange = [2500, 2500+1225] #[2000, 4000] #[2000, 11500]
        # plot full raster
        sim.plotting.plotRaster(
            **{ 'include': ['IC'], #['allCells'], 
                'saveFig': filename[:-4]+'_raster_%d_%d.png'%(timeRange[0], timeRange[1]), 
                'popColors': {'IC': '#02216b'},
                'showFig': False, 
                'popRates': 'minimal', 
                'orderInverse': True, 
                'timeRange': timeRange, 
                'figSize': (16,5), 
                'lw': 1.0, 
                'markerSize': 1.0, 
                'marker': '.', 
                'dpi': 300})

        ax = plt.gca()

        [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

        #ax.get_legend().remove()

        plt.title('')
        plt.margins(x=0,y=0)

        plt.subplots_adjust(bottom=0.10, top=0.98, right=0.95, left=0.10)

        save_filename='%s_raster_%d_%d.png'%(filename[:-4], timeRange[0], timeRange[1])
        plt.savefig(save_filename, dpi=300)  

    # --------------------------------------------
    # Spike Hist
    if plotSpikeHist:
        from netpyne.analysis.spikes_legacy import plotSpikeHist,plotRatePSD

        filename = '../data/v34_batch58/v34_batch58_0_0_0_0_0_data.pkl' # speech
        sim.load(filename, instantiate=False)

        includeList = [['IC', thal, cortexE]] #[allpops, excpops, ['allCells'], [supra, gran, infra], [supraE, granE, infraE]]
        includeLabels = ['IC_thal_cxE'] #['eachpop', 'excpops', 'allCells', 'supra_gran_infra', 'E_supra_gran_infra']

        timeRange = [2500, 2500+1225]  # for adaptation measure 

        for include, includeLabel in zip(includeList, includeLabels):

            fontsize = 14
            binSize = 5
            measure = 'count'
            graphType = 'bar' #'line'                

            plt.set_cmap('tab10')
            ''''
            tab10:
            #1f77b4 - blue
            #ff7f0e - orange
            #2ca02c - green
            #d62728 - red
            #9467bd - purple
            #8c564b
            #e377c2
            #7f7f7f
            #bcbd22
            #17becf

            '''

            fig_hist, data_hist = plotSpikeHist(include=include, 
                                    popColors=['#02216b', '#ff7f0e', '#2ca02c'], #{'IC': 'blue', tuple(thal): 'green', tuple(cortexE): 'red'}, 
                                    lineWidth=2, 
                                    lw=2,
                                    fontSize=fontsize,
                                    timeRange=timeRange, 
                                    dpi=300, 
                                    figSize=(16,5), 
                                    legend=False, 
                                    showFig=0, 
                                    saveFig=0, #filename[:-4]+'_spikehist_%d_%d.png'%(timeRange[0], timeRange[1]),
                                    binSize=binSize, 
                                    graphType=graphType, 
                                    axis=True, 
                                    measure=measure, 
                                    scalebarLoc='upper center')

            ax = plt.gca()

            for line in ax.lines:
                line.set_linewidth(1.5)

            [i.set_linewidth(0.5) for i in ax.spines.values()] # make border thinner

            ylim = ax.get_ylim()

            plt.ylim(0,150)

            #plt.xticks([center-100, center, center+100], ['-100', '0', '100'], fontsize=fontsize, fontname='Arial')
            #plt.xlim(timeRange[0], timeRange[1])
            plt.ylabel('Spike count', fontsize=fontsize, fontname='Arial') #Neurons (ordered by NCD within each pop)')
            plt.xlabel('Time (ms)', fontsize=fontsize, fontname='Arial')
            
            leg = ax.get_legend()
            leg.set_title('')
            
            leg.set_bbox_to_anchor((1.02, 1.0))

            for text in leg.get_texts():
                text.set_fontsize(fontsize)
                text.set_fontname('Arial')
            
            plt.title('')
            plt.margins(x=0,y=0)

            plt.subplots_adjust(bottom=0.1, top=0.98, right=0.90 , left=0.05)

            save_filename = 'figs/figS5_%s_spikehist_%d_%d_%s.png'%(filename.split('/')[-1][:-4], timeRange[0], timeRange[1], includeLabel)
            print(save_filename)
            plt.savefig(save_filename, dpi=300)






# --------------------------
# Main
#
# Select what figure to plot by uncommenting the corresponding lines
#
# --------------------------
if __name__ == '__main__':

    # -------------------------------------------
    # Figure 2 (conn matrix)
    # plot_empirical_conn()
  
    # -------------------------------------------    
    # Figure 3 (multiscale measures)
    fig_raster('v35_batch9', 'v35_batch9_0_0_data', timeRange=[1000,2000])
    # fig_traces('v34_batch27', 'v34_batch27_0_0', timeRange=[1000,2000])
    # fig_stats('v35_batch9', 'v35_batch9_0_0_data', timeRange=[1000,2000])
    # fig_CSD('v34_batch27', 'v34_batch27_0_0', timeRange=[1000,2000])
    # fig_dipole('v35_batch9', 'v35_batch9_0_0_data', timeRange=[1000,2000])
    # fig_eeg('v35_batch9', 'v35_batch9_0_0_data', timeRange=[1000,2000])


    # -------------------------------------------
    # Figure 4 (laminar LFP/CSD spont and BBN)
    # fig_CSD_comparison()          # Fig 4A,D 
    # fig_BBN()                     # Fig 4B
    # fig_BBN_ERP(usemodel = 1)     # Fig 4C
    # fig_BBN_ERP(usemodel = 0)     # Fig 4C

    # -------------------------------------------
    # Figure 5 (LFP PSD comparison)
    # fig_LFP_PSD_comparison()    

    # -------------------------------------------
    # Figure 6 (CSD osc event examples and stats)
    # fig_osc_events()  # Fig 6A
    # fig_osc_stats()   # Fig 6B

    # -------------------------------------------
    # Figure 7 (analysis of CSD osc event)
    # fig_beta_osc_event()                   # Fig 7A,B
    # fig_CSD_contrib_PSD()                 # Fig 7C
    # fig_CSD_osc_spiking()                 # Fig 7D
    # /input_analysis/scripts/plotSPTH.py   # Fig 7F


    # -------------------------------------------
    # Figures S1 (fI curve)
    # run batch using b = fIcurve(pops=[...]) in batch.py
    # fig_fI_curve()

    # -------------------------------------------
    # Figures S2
    # fig_optuna_fitness()

    # -------------------------------------------
    # Figures S3 (UR_EAR screenshot)

    # -------------------------------------------
    # Figure S4 (supplemental; speech)
    # fig_speech()
    # fig_CSD_comparison()  # with speech = 1

    # -------------------------------------------
    # Figure S5 (supplemental; osc event analysis)
    # import csd_analysis