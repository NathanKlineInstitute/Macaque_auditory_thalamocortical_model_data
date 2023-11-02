"""
dataAnalysis.py 

Code to anlayze non-human primate data 

Contributors: ericaygriffith@gmail.com, samnemo@gmail.com
"""


try:
    basestring
except NameError:
    basestring = str


## IMPORTS ## 
import sys
import os
import shutil 
import h5py									              # for rdmat() and getTriggerTimes()
import numpy as np
import downsample
from collections import OrderedDict
from filter import lowpass,bandpass 		  # for getbandpass()
import scipy                              # for plotCSD()
import matplotlib                         # for plotCSD()
from matplotlib import pyplot as plt      # for plotCSD() 
# trying to get sort()
from pylab import *
import pickle as pkl
from netpyne import sim

## PRE-PROCESSING FUNCTIONS ## 
### Originally in rdmat.py ### 
def rdmat (fn,samprds=0):  
  fp = h5py.File(fn,'r') # open the .mat / HDF5 formatted data
  sampr = fp['craw']['adrate'][0][0] # original sampling rate
  print('fn:',fn,'sampr:',sampr,'samprds:',samprds)
  dt = 1.0 / sampr # time-step in seconds
  dat = fp['craw']['cnt'] # cnt record stores the electrophys data
  npdat = np.zeros(dat.shape)
  tmax = ( len(npdat) - 1.0 ) * dt # use original sampling rate for tmax - otherwise shifts phase
  dat.read_direct(npdat) # read it into memory; note that this LFP data usually stored in microVolt
  npdat *= 0.001 # convert microVolt to milliVolt here
  fp.close()
  if samprds > 0.0: # resample the LFPs
    dsfctr = sampr/samprds
    dt = 1.0 / samprds
    siglen = max((npdat.shape[0],npdat.shape[1]))
    nchan = min((npdat.shape[0],npdat.shape[1]))
    npds = [] # zeros((int(siglen/float(dsfctr)),nchan))
    # print dsfctr, dt, siglen, nchan, samprds, ceil(int(siglen / float(dsfctr))), npds.shape
    for i in range(nchan): 
      print('resampling channel', i)
      npds.append(downsample.downsample(npdat[:,i], sampr, samprds))
    npdat = np.array(npds)
    npdat = npdat.T
    sampr = samprds
  tt = np.linspace(0,tmax,len(npdat)) # time in seconds
  return sampr,npdat,dt,tt # npdat is LFP in units of milliVolt


# bandpass filter the items in lfps. lfps is a list or numpy array of LFPs arranged spatially by column
### Originally in load.py ### 
def getbandpass (lfps,sampr,minf=0.05,maxf=300):
  datband = []
  for i in range(len(lfps[0])): datband.append(bandpass(lfps[:,i],minf,maxf,df=sampr,zerophase=True))
  datband = np.array(datband)
  return datband


# Vaknin correction for CSD analysis
# Allows CSD to be performed on all N contacts instead of N-2 contacts
# See Vaknin et al (1989) for more details
### Originally in load.py ### 
def Vaknin(x):
    # Preallocate array with 2 more rows than input array
    x_new = np.zeros((x.shape[0]+2, x.shape[1]))

    # Duplicate first and last row of x into first and last row of x_new
    x_new[0, :] = x[0, :]
    x_new[-1, :] = x[-1, :]

    # Duplicate all of x into middle rows of x_new
    x_new[1:-1, :] = x

    return x_new


def removemean(x, ax=1):
    """
    Function to subtract the mean from an array or list

    Parameters
    ----------
    x : array 
        Data to be processed.
        **Default:** *required*

    ax : int
        The axis to remove the mean across.
        **Default:** ``1``


    Returns
    -------
    data : array
        The processed data.

    """
  
    mean = np.mean(x, axis=ax, keepdims=True)
    x -= mean


# get CSD - first do a lowpass filter. lfps is a list or numpy array of LFPs arranged spatially by column
def getCSD (lfps,sampr,spacing_um,minf=0.05,maxf=300, vaknin=False, norm=True):

  # convert from uV to mV ## Does this already happen in rdmat() ?? 
  # lfps = lfps/1000

  # Bandpass filter the LFP data with getbandpass() fx defined above
  datband = getbandpass(lfps,sampr,minf,maxf)
  
  # Take CSD along smaller dimension
  if datband.shape[0] > datband.shape[1]: # take CSD along smaller dimension
    ax = 1
  else:
    ax = 0

  # Vaknin correction
  if vaknin:
    datband = Vaknin(datband)

  if norm:
    removemean(datband,ax=ax)


  # Convert spacing from microns to um 
  spacing_mm = spacing_um/1000
  # when drawing CSD make sure that negative values (depolarizing intracellular current) drawn in red,
  # and positive values (hyperpolarizing intracellular current) drawn in blue
  CSD = -np.diff(datband,n=2,axis=ax)/spacing_mm**2 # now each column (or row) is an electrode -- CSD along electrodes

  return CSD, datband


#
def getTriggerTimes (fn):
  fp = h5py.File(fn,'r')
  hdf5obj = fp['trig/anatrig']
  x = np.array(fp[hdf5obj.name])
  val = [y[0] for y in fp[x[0,0]].value]
  fp.close()
  return val  

def getTriggerTimesNEW (fn):
  fp = h5py.File(fn,'r')
  hdf5obj = fp['trig/anatrig']
  #x = np.array(fp[hdf5obj.name])   # fp[hdf5obj.name] == hdf5obj
  #x = np.array(hdf5obj)  # NEW LINE
  x = hdf5obj
  z = np.array(fp[x[0,0]])  # NEW LINE 
  #val = [y[0] for y in fp[x[0,0]].value]
  #z[0,0] and z[1,0] -- thru z[799,0]
  val = []
  for i in range(len(z)):
    val.append(z[i,0])
  #val = [y[0] for y in z[0]]
  fp.close()
  return val  

#
def loadfile (fn,samprds,spacing_um=100,vaknin=False):
  # load a .mat data file (fn) using sampling rate samprds (should be integer factor of original sampling rate (44000)),
  # returns: sampr is sampling rate after downsampling
  #          LFP is laminar local field potential data
  #          dt is time-step (redundant with sampr)
  #          tt is time array (in seconds)
  #          CSD is laminar current source density signal
  #          trigtimes is array of stimulus trigger indices (indices into arrays)
  sampr,LFP,dt,tt=rdmat(fn,samprds=samprds) # # samprds = 11000.0 # downsampling to this frequency
  sampr,dt,tt[0],tt[-1] # (2000.0, 0.001, 0.0, 1789.1610000000001)
  CSD, LFP_filtered = getCSD(LFP,sampr,spacing_um, vaknin=vaknin)
  divby = 44e3 / samprds
  trigtimes = None
  try: # not all files have stimuli
    trigtimes = [int(round(x)) for x in np.array(getTriggerTimesNEW(fn)) / divby] # divby since downsampled signals by factor of divby
  except:
    pass
  #trigIDs = getTriggerIDs(fn)
  LFP = LFP.T # make sure each row is a channel
  return sampr,LFP_filtered,dt,tt,CSD,trigtimes


### AVERAGING FUNCTIONS ###
#### NOTE: should also make these available for use for sim data as well in netpyne 
def ms2index (ms, sampr): return int(sampr*ms/1e3)


### REMOVE BAD EPOCHS ### 
def calPosThresh(dat, sigmathresh):
  #return dat.mean() + sigmathresh * dat.std() # this was originally commented out, reversing it to try (EYG)
  return 100

def calNegThresh(dat, sigmathresh):
  #return dat.mean() - sigmathresh * dat.std() # this was originally commented out, reversing it to try (EYG)
  return -100

# remove noise, where noise < negthres < dat < posthres < noise
def badEpoch (dat, sigmathresh):
  badValues = len(np.where(dat <= calNegThresh(dat, sigmathresh))[0]) + \
              len(np.where(dat >= calPosThresh(dat, sigmathresh))[0])
  if badValues > 0:
    return True
  else:
    return False

def removeBadEpochs (dat, sampr, trigtimes, swindowms, ewindowms, sigmathresh):
  nrow = dat.shape[0]
  swindowidx = ms2index(swindowms,sampr) # could be negative
  ewindowidx = ms2index(ewindowms,sampr) # should this be sampr or samprds -- question from EYG (12/11/2020)

  # trigByChannel could be returned for removing different epochs on each channel
  trigByChannel = [x for x in range(nrow)]
  badEpochs = []
  for chan in range(nrow): # go through channels
    trigByChannel[chan] = []
    for trigidx in trigtimes: # go through stimuli
      sidx = max(0,trigidx+swindowidx)
      eidx = min(dat.shape[1],trigidx+ewindowidx)
      if not badEpoch(dat[chan, sidx:eidx], sigmathresh):
        trigByChannel[chan].append(trigidx)
      else:
        badEpochs.append(trigidx)
    print('Found %d bad epochs in channel %d. Range: [%.2f, %.2f]'%
          (len(trigtimes) - len(trigByChannel[chan]), chan,
           calNegThresh(dat[chan, sidx:eidx], sigmathresh),
           calPosThresh(dat[chan, sidx:eidx], sigmathresh)))

  # combine bad epochs into a single sorted list (without duplicates)
  badEpochs = sort(list(set(badEpochs)))
  print('%d bad epochs:'%len(badEpochs),[x for x in badEpochs])

  # remove the associated trigger times before returning
  trigtimes = np.delete(trigtimes,[trigtimes.index(x) for x in badEpochs])

  return trigtimes

### 




# get the average ERP (dat should be either LFP or CSD)
# originally from load.py 
def getAvgERP (dat, sampr, trigtimes, swindowms, ewindowms):
  nrow = dat.shape[0] # number of channels 
  tt = np.linspace(swindowms, ewindowms,ms2index(ewindowms - swindowms,sampr))
  avgERP = np.zeros((nrow,len(tt))) # set up array for averaged values 
  swindowidx = ms2index(swindowms,sampr) # could be negative
  ewindowidx = ms2index(ewindowms,sampr)
  for chan in range(nrow): # go through channels
    for trigidx in trigtimes: # go through stimuli
      sidx = max(0,trigidx+swindowidx)
      eidx = min(dat.shape[1],trigidx+ewindowidx)
      avgERP[chan,:] += dat[chan, sidx:eidx] # add together data points from each time window
    avgERP[chan,:] /= float(len(trigtimes)) # average each data point 
  return tt,avgERP


def getIndividualERP(dat,sampr,trigtimes,swindowms,ewindowms,ERPindex):
  nrow = dat.shape[0] # number of channels 
  tt = np.linspace(swindowms, ewindowms,ms2index(ewindowms - swindowms,sampr))
  individualERP = np.zeros((nrow,len(tt))) # set up array for averaged values 
  swindowidx = ms2index(swindowms,sampr) # could be negative
  ewindowidx = ms2index(ewindowms,sampr)
  trigidx = trigtimes[ERPindex] # which stimulus to look at 
  for chan in range(nrow): # go through channels
    sidx = max(0,trigidx+swindowidx)
    eidx = min(dat.shape[1],trigidx+ewindowidx)
    individualERP[chan,:] = dat[chan, sidx:eidx] # (don't?) add together data points from each time window
  return tt,individualERP


### PLOTTING FUNCTIONS ### 
# PLOT CSD 
def plotCSD(dat,tt,fn=None,saveFolder=None,overlay=None,LFP_data=None,timeRange=None,trigtimes=None,saveFig=True,showFig=True, fontSize=12, figSize=(8,8), layerLines=False, layerBounds=None, smooth=None, dpi=300):
  ## dat --> CSD data as numpy array
  ## timeRange --> time range to be plotted (in ms)
  ## trigtimes --> trigtimes from loadfile() (indices -- must be converted)
  ## tt --> numpy array of time points (time array in seconds)
  ## fn --> filename -- string -- used for saving! 
  ## saveFolder --> string -- path to directory where figs should be saved
  ## overlay --> can be 'LFP' or 'CSD' to overlay time series of either dataset 
  ## LFP_data --> numpy array of LFP data 
  ## layerLines : bool 
    #     Whether to plot horizontal lines over CSD plot at layer boundaries. 
    #     **Default:** ``False`` 

  ## layerBounds : dict
    #     Dictionary containing layer labels as keys, and layer boundaries as values, e.g. {'L1':100, 'L2': 160, 'L3': 950, 'L4': 1250, 'L5A': 1334, 'L5B': 1550, 'L6': 2000}
    #     **Default:** ``None``
  
  if trigtimes is not None:
    trigtimesMS = []                # GET TRIGGER TIMES IN MS -- convert trigtimes to trigtimesMS (# NOTE: SHOULD MAKE THIS A FUNCTION)
    for idx in trigtimes:
      trigtimesMS.append(tt[idx]*1e3)
  else:
    trigtimesMS = None



  tt = tt*1e3 # Convert units from seconds to ms 
  dt = tt[1] - tt[0] # dt is the time step of the recording # UNITS: in ms after converstion
  print('tt --> last time point: ' + str(tt[-1]))

  if timeRange is None:
    timeRange = [0,tt[-1]] # if timeRange is not specified, it takes the entire time range of the recording (ms)
  else:
    dat = dat[:,int(timeRange[0]/dt):int(timeRange[1]/dt)] # SLICE CSD DATA APPROPRIATELY
    tt = tt[int(timeRange[0]/dt):int(timeRange[1]/dt)] # DO THE SAME FOR TIME POINT ARRAY 
    #tt = np.arange(timeRange[0], timeRange[1], dt)
    print('tt --> new, last time point: ' + str(tt[-1]))

  # PLOTTING
  X = tt 
  Y = np.arange(dat.shape[0]) 

  # interpolation 
  CSD_spline = scipy.interpolate.RectBivariateSpline(Y,X,dat)
  Y_plot = np.linspace(0,dat.shape[0]-1,num=1000)  # dat.shape[0] - 1 --> cuts off the noise at the bottom 
  Z = CSD_spline(Y_plot,X) 

  # (i) Set up axes
  plt.rcParams.update({'font.size': fontSize}) # this is from plotCSD in netpyne -- is fontSize valid here? 
  xmin = int(X[0])
  xmax = int(X[-1]) + 1 
  ymin = 0 #1   # 0 in csd.py in netpyne 
  ymax = dat.shape[0]-1 #21 #24 #20   #22 # 24 in csd_verify.py, but it is spacing in microns in csd.py in netpyne --> WHAT TO DO HERE? TRY 24 FIRST! 
  extent_xy = [xmin, xmax, ymax, ymin]

  # (ii) Set up figure
  fig = plt.figure(figsize=figSize) # plt.figure(figsize = figSize) <-- would have to add figSize as an argument above 
  ax = plt.gca()
  ax.grid(False)

  # (iii) Create plots w/ common axis labels and tick marks
  axs = []
  numplots = 1 # HOW TO DETERMINE THIS? WHAT IF MORE THAN ONE? 
  gs_outer = matplotlib.gridspec.GridSpec(1,1) #matplotlib.gridspec.GridSpec(2, 2, figure=fig, wspace=0.4, hspace=0.2, height_ratios = [20, 1])

  for i in range(numplots):
    axs.append(plt.Subplot(fig,gs_outer[i*2:i*2+2]))
    fig.add_subplot(axs[i])
    axs[i].set_xlabel('Time (ms)',fontsize=fontSize) #fontsize=12)
    #axs[i].tick_params(axis='y', which='major', labelsize=8)
    axs[i].set_yticks(np.arange(ymin, ymax, step=3))
    #axs[i].tick_params(axis='y', which='major', labelsize=fontSize)
    axs[i].tick_params(axis='x', which='major', labelsize=fontSize)


  # SMOOTHING
  if smooth:
    Z = scipy.ndimage.filters.gaussian_filter(Z, smooth, mode='nearest')

  # (iv) PLOT INTERPOLATED CSD COLOR MAP
  spline=axs[0].imshow(Z, extent=extent_xy, interpolation='none', aspect='auto', origin='upper', cmap='jet_r', alpha=0.9) # alpha controls transparency -- set to 0 for transparent, 1 for opaque
  axs[0].set_ylabel('Channel', fontsize = fontSize) # Contact depth (um) -- convert this eventually 


  # (v) OVERLAY -- SETTING ASIDE FOR NOW -- THAT IS NEXT GOAL 
  if overlay is None:
    print('No data being overlaid')
    axs[0].set_title('NHP Current Source Density (CSD)', fontsize=fontSize)

  elif overlay == 'CSD' or overlay == 'LFP':
    nrow = dat.shape[0] # number of channels
    gs_inner = matplotlib.gridspec.GridSpecFromSubplotSpec(nrow, 1, subplot_spec=gs_outer[0:2], wspace=0.0, hspace=0.0) 
    subaxs = []

    # go down grid and add data from each channel 
    if overlay == 'CSD':
      axs[0].set_title('NHP CSD with time series overlay', fontsize=fontSize)
      for chan in range(nrow):
          subaxs.append(plt.Subplot(fig,gs_inner[chan],frameon=False))
          fig.add_subplot(subaxs[chan])
          subaxs[chan].margins(0.0,0.01)
          subaxs[chan].get_xaxis().set_visible(False)
          subaxs[chan].get_yaxis().set_visible(False)
          subaxs[chan].plot(X,dat[chan,:],color='green',linewidth=0.3)

    elif overlay == 'LFP':
      if LFP_data is None:
        print('no LFP data provided!')
      else:
        axs[0].set_title('NHP CSD with LFP overlay', fontsize=fontSize)
        LFP_data = LFP_data[:,int(timeRange[0]/dt):int(timeRange[1]/dt)] # slice LFP data according to timeRange
        for chan in range(nrow):
          subaxs.append(plt.Subplot(fig,gs_inner[chan],frameon=False))
          fig.add_subplot(subaxs[chan])
          subaxs[chan].margins(0.0,0.01)
          subaxs[chan].get_xaxis().set_visible(False)
          subaxs[chan].get_yaxis().set_visible(False)
          subaxs[chan].plot(X,LFP_data[chan,:],color='gray',linewidth=0.3) #LFP_data[:,chan] in netpyne csd 


  ## ADD ARROW POINTING TO STIMULUS TIMES      ## NOTE: this was taken from (v) in plotAvgCSD() & main code block
  if trigtimesMS is not None:
    for time in trigtimesMS:
      if time >= xmin and time <= xmax:
        axs[0].annotate(' ', xy = (time,24), xytext = (time,24), arrowprops = dict(facecolor='red', shrink=0.1, headwidth=4,headlength=4),annotation_clip=False)
        axs[0].vlines(time,ymin=ymin,ymax=ymax,linestyle='dashed')


  if layerLines:
    if layerBounds is None:
      print('No layer boundaries given -- will not overlay layer boundaries on CSD plot')
    else:
      layerKeys = []
      for i in layerBounds.keys():
        axs[0].hlines(layerBounds[i], xmin, xmax, colors='black', linewidth=1, linestyles='dotted')
        layerKeys.append(i) # makes a list with names of each layer, as specified in layerBounds dict argument 

      for n in range(len(layerKeys)): # label the horizontal layer lines with the proper layer label 
        if n==0:
          axs[0].text(xmax+5, layerBounds[layerKeys[n]]/2, layerKeys[n], color='black', fontsize=fontSize)
        else:
          axs[0].text(xmax+5, (layerBounds[layerKeys[n]] + layerBounds[layerKeys[n-1]])/2, layerKeys[n], color='black', fontsize=fontSize, verticalalignment='center')



  # SAVE FIGURE
  ax = plt.gca()
  ax.grid(False)

  if saveFig:
    if isinstance(saveFig, basestring):
      filename = saveFig
    else:
      filename = 'CSD_fig.png'
    try:
      plt.savefig(filename, dpi=dpi)
    except:
      plt.savefig('CSD_fig.png', dpi=dpi)

  # if saveFig:
  #   if fn is None:
  #     if overlay == 'LFP':
  #       figname = 'NHP_CSD_withLFP.png'
  #     elif overlay == 'CSD':
  #       figname = 'NHP_CSD_csdOverlay.png'
  #     else:
  #       figname = 'NHP_CSD_fig.png'
  #   else:
  #     #filename = fn[31:-4] # takes out the .mat from the filename given as arg
  #     filename = fn.split("/")[6][:-4] # filename without any of the path info or the .mat suffix 
  #     if overlay == 'LFP':
  #       figname = 'NHP_CSD_withLFP_%s.png' % filename
  #     elif overlay == 'CSD':
  #       figname = 'NHP_CSD_csdOverlay_%s.png' % filename
  #     else:
  #       figname = 'NHP_CSD_fig_%s.png' % filename
  #   try:
  #     plt.savefig(figname, dpi=dpi) # dpi
  #   except:
  #     plt.savefig('NHP_CSD_fig.png', dpi=dpi)

    # move saved fig to appropriate directory 
    if saveFolder is not None:
      if not os.path.isdir(saveFolder):
        cwd = os.getcwd()
        origFigPath = cwd + '/' + figname
        newFigPath = saveFolder + figname
        os.mkdir(saveFolder)
        shutil.move(origFigPath, newFigPath)




  # DISPLAY FINAL FIGURE
  if showFig is True:
    plt.show()
    #plt.close()



### PLOT CSD OF AVERAGED ERP ### 
def plotAvgCSD(dat,tt,trigtimes=None,fn=None,saveFolder=None,overlay=True,saveFig=True,showFig=True):
  ## dat --> CSD data as numpy array (from getAvgERP)
  ## tt --> numpy array of time points (from getAvgERP)
  ## trigtimes --> should be relativeTrigTimesMS (NOTE: trigtimes changed from relativeTrigTimes=None)
  ## fn --> string --> input filename --> used for saving! 
  ## saveFolder --> string --> directory where figs should be saved
  ## overlay --> Default TRUE --> plots avgERP CSD time series on top of CSD color map 

  # INTERPOLATION
  X = tt 
  Y = np.arange(dat.shape[0]) # number of channels
  CSD_spline = scipy.interpolate.RectBivariateSpline(Y,X,dat)
  Y_plot = np.linspace(0,dat.shape[0],num=1000) # ,num=1000 is included in csd.py in netpyne --> hmm. necessary? 
  Z = CSD_spline(Y_plot,X)

  # (i) Set up axes
  xmin = int(X[0])
  xmax = int(X[-1]) + 1 
  ymin = 1  # 0 in csd.py in netpyne 
  ymax = 21 # 22  # dat.shape[0] # 24 in csd_verify.py, but it is spacing in microns in csd.py in netpyne --> WHAT TO DO HERE? TRY 24 FIRST! 
  extent_xy = [xmin, xmax, ymax, ymin]

  # (ii) Set up figure
  fig = plt.figure(figsize=(5,8))

  # (iii) Create plots w/ common axis labels and tick marks
  axs = []

  numplots = 1 # HOW TO DETERMINE THIS? WHAT IF MORE THAN ONE? 

  gs_outer = matplotlib.gridspec.GridSpec(2, 2, figure=fig, wspace=0.4, hspace=0.2, height_ratios = [20, 1])

  for i in range(numplots):
    axs.append(plt.Subplot(fig,gs_outer[i*2:i*2+2]))
    fig.add_subplot(axs[i])
    axs[i].set_xlabel('Time (ms)',fontsize=12)
    #axs[i].tick_params(axis='y', which='major', labelsize=8)
    axs[i].set_yticks(np.arange(ymin, ymax, step=1))

  # (iv) PLOT INTERPOLATED CSD COLOR MAP
  spline=axs[0].imshow(Z, extent=extent_xy, interpolation='none', aspect='auto', origin='upper', cmap='jet_r', alpha=0.9) # alpha controls transparency -- set to 0 for transparent, 1 for opaque
  axs[0].set_ylabel('Channel', fontsize = 12) # Contact depth (um) -- convert this eventually 

  #########
  # (v) ADD ARROW POINTING TO STIMULUS TIMES
    # OLD CODE #### DEPRECATED #### ## ADD ARG THEN IF STATEMENT
    # for time in relativeTrigTimes:
    #   if time <= xmax:
    #     axs[0].annotate(' ', xy = (time,24), xytext = (time,24), arrowprops = dict(facecolor='red', shrink=0.1, headwidth=4,headlength=4),annotation_clip=False)
    #     axs[0].vlines(time,ymin=ymin,ymax=ymax,linestyle='dashed')


  # (v) ADD ARROW POINTING TO STIMULUS TIMES        ## NOTE: this was taken from (v) in plotAvgCSD() & main code block
  # if trigtimes: ## FIGURE THIS OUT 
  
  # ## GET TRIGGER TIMES IN MS
  # trigTimesMS = []
  # for idx in trigtimes:
  #   trigTimesMS.append(tt[idx]*1e3)
  # #print(trigTimesMS) # USEFUL FOR KNOWING ABSOLUTE VALUE OF STIMULUS TIMES

  # ## GET RELATIVE TRIGGER TIME INDICES
  # relativeTrigTimes = []
  # for idx in trigtimes:
  #   relativeTrigTimes.append(idx - trigtimes[0])#80143)

  # ## GET RELATIVE TRIGGER TIMES IN SECONDS 
  # relativeTrigTimesMS = []
  # for time in trigTimesMS:
  #   relativeTrigTimesMS.append(time-trigTimesMS[0])

    # for time in trigtimesMS:
    #   if time >= xmin and time <= xmax:
    #     axs[0].annotate(' ', xy = (time,24), xytext = (time,24), arrowprops = dict(facecolor='red', shrink=0.1, headwidth=4,headlength=4),annotation_clip=False)
    #     axs[0].vlines(time,ymin=ymin,ymax=ymax,linestyle='dashed')

  #########

  # (vi) SET TITLE AND OVERLAY AVERAGE ERP TIME SERIES (OR NOT)
  ## NOTE: add option for overlaying average LFP...??
  if overlay:
    nrow = dat.shape[0] # number of channels 
    gs_inner = matplotlib.gridspec.GridSpecFromSubplotSpec(nrow, 1, subplot_spec=gs_outer[0:2], wspace=0.0, hspace=0.0)
    subaxs = []

    # set title
    axs[0].set_title('NHP CSD with CSD time series overlay',fontsize=12)
    # go down grid and add data from each channel
    for chan in range(nrow):
        subaxs.append(plt.Subplot(fig,gs_inner[chan],frameon=False))
        fig.add_subplot(subaxs[chan])
        subaxs[chan].margins(0.0,0.01)
        subaxs[chan].get_xaxis().set_visible(False)
        subaxs[chan].get_yaxis().set_visible(False)
        subaxs[chan].plot(X,dat[chan,:],color='red',linewidth=0.3)

  else:
    axs[0].set_title('NHP Current Source Density (CSD)', fontsize=14)


  #if overlay_individualERP:


  # SAVE FIGURE
  ## make this a little more explicable 
  if saveFig:
    if fn is None: 
      if overlay:
        figname = 'NHP_avgCSD_csdOverlay.png'
      else:
        figname = 'NHP_avgCSD.png'
    else: 
      filename = fn.split("/")[6][:-4] # filename without any of the path info or the .mat suffix 
      print(filename)
      if overlay:  
        figname = 'NHP_avgCSD_csdOverlay_%s.png' % filename
      else:
        figname = 'NHP_avgCSD_%s.png' % filename
    try:
      plt.savefig(figname) #dpi
    except:
      plt.savefig('NHP_avgCSD.png')


  # DISPLAY FINAL FIGURE
  if showFig is True:
    plt.show()
    #plt.close()


def plotIndividualERP(dat,tt,saveFig=False,showFig=True): # trigtimes,
  # trigtines unnecessary argument? 
  # Add fn as argument? 

  # dat should be individual ERP from getIndividual ERP
  # tt should be ttavg from getAvgERP
  # trigtimes should be relativeTrigTimesMS -- NOT EVEN USED RIGHT NOW; Pre-emptive code below 


  # Convert trigtimes (list of indices in tt that correspond to onset of trigger stimuli) to trigtimesMS and/or relativeTrigTimesMS
  # (i) ## GET TRIGGER TIMES IN MS
  # trigTimesMS = []
  # for idx in trigtimes:
  #   trigTimesMS.append(tt[idx]*1e3)

  # (ii) ## GET RELATIVE TRIGGER TIME INDICES
  # relativeTrigTimes = []
  # for idx in trigtimes:
  #   relativeTrigTimes.append(idx - trigtimes[0])#80143)

  # (iii) ## GET RELATIVE TRIGGER TIMES IN MILLISECONDS 
  # relativeTrigTimesMS = []
  # for time in trigTimesMS:
  #   relativeTrigTimesMS.append(time-trigTimesMS[0])

  ### CODE FROM AVG CSD #### 
  # INTERPOLATION
  X = tt 
  Y = np.arange(dat.shape[0]) # number of channels
  CSD_spline = scipy.interpolate.RectBivariateSpline(Y,X,dat)
  Y_plot = np.linspace(0,dat.shape[0],num=1000) # ,num=1000 is included in csd.py in netpyne --> hmm. necessary? 
  Z = CSD_spline(Y_plot,X)

  # (i) Set up axes
  xmin = int(X[0])
  xmax = int(X[-1]) + 1 
  ymin = 1  # 0 in csd.py in netpyne 
  ymax = 21 # 22  # dat.shape[0] # 24 in csd_verify.py, but it is spacing in microns in csd.py in netpyne --> WHAT TO DO HERE? TRY 24 FIRST! 
  extent_xy = [xmin, xmax, ymax, ymin]

  # (ii) Set up figure
  fig = plt.figure(figsize=(5,8))

  # (iii) Create plots w/ common axis labels and tick marks
  axs = []

  numplots = 1 # HOW TO DETERMINE THIS? WHAT IF MORE THAN ONE? 

  gs_outer = matplotlib.gridspec.GridSpec(2, 2, figure=fig, wspace=0.4, hspace=0.2, height_ratios = [20, 1])

  for i in range(numplots):
    axs.append(plt.Subplot(fig,gs_outer[i*2:i*2+2]))
    fig.add_subplot(axs[i])
    axs[i].set_xlabel('Time (ms)',fontsize=12)
    #axs[i].tick_params(axis='y', which='major', labelsize=8)
    axs[i].set_yticks(np.arange(ymin, ymax, step=1))

  # (iv) PLOT INTERPOLATED CSD COLOR MAP
  spline=axs[0].imshow(Z, extent=extent_xy, interpolation='none', aspect='auto', origin='upper', cmap='jet_r', alpha=0.9) # alpha controls transparency -- set to 0 for transparent, 1 for opaque
  axs[0].set_ylabel('Channel', fontsize = 12) # Contact depth (um) -- convert this eventually 

  ####### ### OLD CODE BEGINNING ###### 

  # ### MAKE FIGURE OF INDIVIDUAL ERP ### --- NOTE: are we sure this is what that does? VERIFY. 
 
  # fig = plt.figure()
  # nrow = dat.shape[0] # number of channels
  # axs = []


  # gs_outer = matplotlib.gridspec.GridSpec(2, 2, figure=fig, wspace=0.4, hspace=0.2, height_ratios = [20, 1])
  # gs_inner = matplotlib.gridspec.GridSpecFromSubplotSpec(nrow, 1, subplot_spec=gs_outer[0:2], wspace=0.0, hspace=0.0)

  # for chan in range(nrow):
  #   axs.append(plt.Subplot(fig,gs_inner[chan],frameon=False))
  #   fig.add_subplot(axs[chan])
  #   axs[chan].margins(0.0,0.01)
  #   axs[chan].get_xaxis().set_visible(False)
  #   axs[chan].get_yaxis().set_visible(False)
  #   axs[chan].plot(tt,dat[chan,:],color='red',linewidth=0.3)

  # plt.xlabel('Time (ms)')
  # plt.ylabel('Channel')  
  # #plt.show()

  if saveFig:
    figname = 'NHP_individualERP.png'
    plt.savefig(figname) #dpi
    # if fn is None: 
    #   if overlay:
    #     figname = 'NHP_avgCSD_csdOverlay.png'
    #   else:
    #     figname = 'NHP_avgCSD.png'
    # else: 
    #   filename = fn.split("/")[6][:-4] # filename without any of the path info or the .mat suffix 
    #   print(filename)
    #   if overlay:  
    #     figname = 'NHP_avgCSD_csdOverlay_%s.png' % filename
    #   else:
    #     figname = 'NHP_avgCSD_%s.png' % filename
    # try:
    #   plt.savefig(figname) #dpi
    # except:
    #   plt.savefig('NHP_avgCSD.png')


  # DISPLAY FINAL FIGURE
  if showFig is True:
    plt.show()
    #plt.close()



#### PLOT LFP FUNCTIONS ####

# return first line matching s if it exists in file fn
def grepstr (fn, s):
  try:
    fp = open(fn,'r', encoding='utf-8')
    lns = fp.readlines()
    for ln in lns:
      if ln.count(s) > 0:
        fp.close()
        return ln.strip()
    fp.close()
  except:
    pass
  return False

#
def monoinc (lx):
  if len(lx) < 2: return True
  for i,j in zip(lx,lx[1:]):
    if i > j:
      return False
  return True

# this function gets the CSD/LFP channel ranges for the .mat cortical recording:
# s1: supragranular source
# s2: supragranular sink
# g: granular sink
# i1: infragranular sink
# i2: infragranular source
# each of these values have a range, by default will pick the middle value as s1,s2,g,i1,i2
#
# note that indices in dbpath file are Matlab based so subtracts 1 first
# since not all files have layers determined, returns empty values (-1) when not found
# when abbrev==True, only get s2,g,i1
def getflayers (fn, dbpath,getmid=True,abbrev=False):
  if dbpath is None or len(dbpath)==0: dbpath = findcsvdbpath(fn)
  s = grepstr(dbpath,os.path.split(fn)[-1])
  if s == False:
    if abbrev:
      return [-1 for i in range(3)]
    else:
      return [-1 for i in range(5)]
  ls = s.split(',')
  print(ls)
  try:
    lint = [int(x)-1 for x in ls[2:]] #[int(x)-1 for x in ls[2:]] # TRY NOT SUBTRACTING 1 
    if not monoinc(lint):
      if abbrev:
        return [-1 for i in range(3)]
      else:
        return [-1 for i in range(5)]      
    s1low,s1high,s2low,s2high,glow,ghigh,i1low,i1high,i2low,i2high = lint
    if getmid:
      s1 = int((s1low+s1high)/2.0)
      s2 = int((s2low+s2high)/2.0)
      g = int((glow+ghigh)/2.0)
      i1 = int((i1low+i1high)/2.0)
      i2 = int((i2low+i2high)/2.0)
      print(s1low,s1high,s2low,s2high,glow,ghigh,i1low,i1high,i2low,i2high,s1,s2,g,i1,i2)
      if abbrev:
        return s2,g,i1
      else:
        return s1,s2,g,i1,i2
    else:
      return s1low,s1high,s2low,s2high,glow,ghigh,i1low,i1high,i2low,i2high
  except:
    if abbrev:
      return [-1 for i in range(3)]
    else:
      return [-1 for i in range(5)]


# Function to convert   
def listToString(s):  
    # initialize an empty string 
    str1 = ""  
    # traverse in the string   
    for elec in s:  
        str1 += str(elec)   
    # return string   
    return str1  
        
      

# -------------------------------------------------------------------------------------------------------------------
## Plot LFP (time-resolved, power spectral density, time-frequency and 3D locations)
# -------------------------------------------------------------------------------------------------------------------
## ADAPTED FROM NETPYNE ANALYSIS lfp.py 
## may change electrodes 

def plotLFP(dat,tt,timeRange=None,trigtimes=None, triglines=False, electrodes=['avg', 'all'], plots=['timeSeries','spectrogram'], dbpath=None, fn=None, NFFT=256, noverlap=128, nperseg=256, minFreq=1, maxFreq=100, stepFreq=1, smooth=0, separation=1.0, includeAxon=True, logx=False, logy=False, normSignal=False, normPSD=False, normSpec=False, filtFreq=False, filtOrder=3, detrend=False, transformMethod='morlet', maxPlots=8, overlay=False, colors=None, figSize=(8, 8), fontSize=14, lineWidth=1.5, dpi=200, saveData=None, saveFig=None, showFig=True):
  ## dat --> LFP data as numpy array
  ## tt --> numpy array of time points (time array in seconds)
  ## timeRange --> time range to be plotted (in ms)
  ## trigtimes --> trigtimes from loadfile() (indices -- must be converted)
  ## triglines --> boolean; indicates whether or not you want a line at the stimulus onset trigger
  ## electrodes --> list of electrodes to include -- can be numbers, or 'all', 'avg', or 'supra'/'gran'/'infra' if there is a .csv included with dbpath
  ## dbpath --> path to .csv layer file with relevant .mat filename, needed if electrodes includes 'supra', 'infra', or 'gran'
  ## fn --> path to .mat filename, needed if electrodes includes 'supra', 'infra', or 'gran', also works for saving
    #### ARGS NOT YET ADDED / USED: 
    ## saveFolder --> string -- path to directory where figs should be saved
    ## plotByLayer -- break the graphs up by supra/infra/gran or plot all on the same graph 

  try:
      from netpyne.plotting.plotter import add_scalebar, colorList
  except:
      from netpyne.analysis.utils import colorList
      from netpyne.support.scalebar import add_scalebar
  from numbers import Number


  print('Plotting LFP ...')

  if not colors: colors = colorList

  # set font size
  plt.rcParams.update({'font.size': fontSize})

  # Get trigger times (in ms) iff file is speech or click 
  if trigtimes is not None:
    trigtimesMS = []
    for idx in trigtimes:
      trigtimesMS.append(tt[idx]*1e3) # tt already converted to ms


  tt = tt*1e3        # Convert time array units from seconds to ms 
  dt = tt[1] - tt[0] # dt is the time step of the recording (equivalent to sim.cfg.recordStep) # UNITS: in ms after conversion


  # SLICE LFP AND TIME ARRAYS APPROPRIATELY ACCORDING TO timeRange
  if timeRange is None:
    timeRange = [0,tt[-1]] # if timeRange is not specified, it takes the entire time range of the recording (ms)
  else:
    dat = dat[:,int(timeRange[0]/dt):int(timeRange[1]/dt)]  # SLICE LFP DATA ACCORDING TO timeRange
    tt = tt[int(timeRange[0]/dt):int(timeRange[1]/dt)]      # DO THE SAME FOR TIMEPOINT ARRAY 

  lfp = dat # set lfp equal to 'dat' so that rest of code is comparable to netpyne lfp.py 
  #lfp = dat * 1000 ## SWITCH LFP UNITS FROM mV to uV !! 

  if filtFreq:  # default is False
    from scipy import signal
    fs = 1000.0/dt # dt is equivalent to sim.cfg.recordStep 
    nyquist = fs/2.0
    if isinstance(filtFreq, list): # bandpass
      Wn = [filtFreq[0]/nyquist, filtFreq[1]/nyquist]
      b, a = signal.butter(filtOrder, Wn, btype='bandpass')
    elif isinstance(filtFreq, Number): # lowpass
      Wn = filtFreq/nyquist
      b, a = signal.butter(filtOrder, Wn)
    for i in range(lfp.shape[1]):
      lfp[:,i] = signal.filtfilt(b, a, lfp[:,i])

  if detrend:   # default is False 
    from scipy import signal
    for i in range(lfp.shape[1]):
      lfp[:,i] = signal.detrend(lfp[:,i])

  if normSignal:  # default is False
    for i in range(lfp.shape[1]):
      offset = min(lfp[:,i])
      if offset <= 0:
        lfp[:,i] += abs(offset)
      lfp[:,i] /= max(lfp[:,i])


  ### ELECTRODES ###
  nrow = lfp.shape[0]  # number of channels in the recording
  if 'all' in electrodes:
    electrodes.remove('all')
    electrodes.extend(list(range(nrow))) 

  ### PLOTTING ### 
  figs = []
  axs = []


  ### TIME SERIES -----------------------------------------
  if 'timeSeries' in plots:
    ydisp = np.absolute(lfp).max() * separation
    print('ydisp = ' + str(ydisp))
    offset = 1.0*ydisp
    if figSize:
      figs.append(plt.figure(figsize=figSize))


    for i,elec in enumerate(electrodes):
      if elec == 'avg':
        lfpPlot = np.mean(lfp,axis=0) ## axis = 1 in netpyne lfp.py, but dims should be flipped for data
        color='k'
        lw = 1.0 

      ## lfpPlot for supgragranular, infragranular, or granular layer(s). 
      elif elec in ['supra', 'infra', 'gran']:
        if dbpath is None:
          print('need path to .csv layer file')
        if fn is None:
          print('need path to .mat data file')
        else:
          s2,g,i1 = getflayers(fn,dbpath=dbpath,abbrev=True)
          print('supra: ' + str(s2))
          print('gran: ' + str(g))
          print('infra: ' + str(i1))

        if elec == 'supra':
          lfpPlot = lfp[s2,:]
          color='r'
          lw = 1.0 
        elif elec == 'infra':
          lfpPlot = lfp[i1,:]
          color='b'
          lw = 1.0 
        elif elec == 'gran':
          lfpPlot = lfp[g,:]
          color='g'
          lw = 1.0 

      elif isinstance(elec, Number):
        lfpPlot = lfp[elec,:] # this is lfp[:,elec] in netpyne lfp.py, but dims should be flipped for data
        color = colors[i%len(colors)]
        lw = 1.0 

      plt.plot(tt[0:len(lfpPlot)], -lfpPlot+(i*ydisp), color=color, linewidth=lw) # chan changed to i
      plt.text(timeRange[0]-0.07*(timeRange[1]-timeRange[0]), (i*ydisp), elec, color=color, ha='center', va='top', fontsize=fontSize, fontweight='bold') # chan changed to i
    
    ax = plt.gca()

    ## Add vertical lines to stimulus onset    ## trigtimes only if non-spontaneous data (click or speech)
    if triglines:
      if trigtimesMS is not None:
        for trig in trigtimesMS:
          if trig >= timeRange[0] and trig <= timeRange[1]:
            #ax.annotate(' ', xy = (trig,24), xytext = (trig,24), arrowprops = dict(facecolor='red', shrink=0.1, headwidth=4,headlength=4),annotation_clip=False)
            ax.axvline(trig,linestyle='dashed') #ymin=,ymax=,


    ## FORMAT PLOT ### 
    ### x axis 
    plt.xlabel('Time (ms)', fontsize=fontSize)
    ### y axis 
    if len(electrodes) > 1:
      plt.text(timeRange[0]-0.14*(timeRange[1]-timeRange[0]), (len(electrodes)*ydisp)/2.0, 'LFP electrode', color='k', ha='left', va='bottom', fontSize=fontSize, rotation=90)
    plt.ylim(-offset, (len(electrodes))*ydisp)
    ax.invert_yaxis()
    ax.axes.yaxis.set_ticks([])
    ### make top, right, left spines of the plot invisible
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ### add title to graph 
    plt.suptitle('LFP Signal')

    ### ADD SCALEBAR ### 
    ##### RUNTIME WARNING BEING GENERATED -- SEE NB 
    round_to_n = lambda x, n, m: int(np.ceil(round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1)) / m)) * m
    scaley = 1000.0  # values in mV but want to convert to uV (??) -- values of lfp are in uV right now.  
    m = 10.0
    sizey = 100/scaley
    #while sizey > 0.25*ydisp:
    if sizey > 0.25 * ydisp: # Change while to if to resolve runtime overflow warning -- see nb for more details 
      try:
        sizey = round_to_n(0.2*ydisp*scaley, 1, m) / scaley
        print('sizey = ' + str(sizey))
      except:
        sizey /= 10.0
      m /= 10.0
    labely = '%.3g $\mu$V'%(sizey*scaley)#)[1:]
    if len(electrodes) > 1:
      add_scalebar(ax,hidey=True, matchy=False, hidex=False, matchx=False, sizex=0, sizey=-sizey, labely=labely, unitsy='$\mu$V', scaley=scaley, loc=3, pad=0.5, borderpad=0.5, sep=3, prop=None, barcolor="black", barwidth=2)
    else:
      add_scalebar(ax, hidey=True, matchy=False, hidex=True, matchx=True, sizex=None, sizey=-sizey, labely=labely, unitsy='$\mu$V', scaley=scaley, unitsx='ms', loc=3, pad=0.5, borderpad=0.5, sep=3, prop=None, barcolor="black", barwidth=2)


    if saveFig:
      if fn is None:
        electrodeString = listToString(electrodes)
        figname = 'NHP_LFP_timeSeries_electrodes_%s.png' % electrodeString  
      else:
        pathParts = fn.split("/")
        for i,pathPart in enumerate(pathParts):
          if '.mat' in pathPart:
            filenameMat = pathParts[i]
            filename = filenameMat[:-4]
        print('Saving LFP timeSeries for file ' + filename)
        figname = 'NHP_LFP_timeSeries_%s.png' % filename
      try:
        plt.savefig(figname, dpi=dpi) #dpi 
      except:
        plt.savefig('NHP_LFP_timeSeries.png',dpi=dpi)


    if showFig:
      plt.show()

  ### SPECTROGRAM -----------------------------------------
  if 'spectrogram' in plots:
    import matplotlib.cm as cm
    numCols = 1 #np.round(nrow/maxPlots) + 1   # maxPlots = 8
    figs.append(plt.figure(figsize=(figSize[0]*numCols, figSize[1])))


    # Morlet wavelet transform method
    if transformMethod == 'morlet':       # transformMethod is 'morlet' by default
      from netpyne.support.morlet import MorletSpec, index2ms

      spec = []
      freqList = None
      if logy:
        freqList = np.logspace(np.log10(minFreq), np.log10(maxFreq), int((maxFreq-minFreq)/stepFreq))

      # TO PLOT LFP SPECTROGRAMS FROM ALL ELECTRODES 
      for i,elec in enumerate(electrodes):
        if elec=='avg':
          lfpPlot = np.mean(lfp,axis=0)   # axis=1 in netpyne lfp.py, but dims should be flipped for data 

        elif elec in ['supra', 'infra', 'gran']:
          if dbpath is None:
            print('need path to .csv layer file')
          if fn is None:
            print('need path to .mat data file')
          else:
            s2,g,i1 = getflayers(fn,dbpath=dbpath,abbrev=True)
            print('supra: ' + str(s2))
            print('gran: ' + str(g))
            print('infra: ' + str(i1))

          if elec == 'supra':
            lfpPlot = lfp[s2,:]
          elif elec == 'infra':
            lfpPlot = lfp[i1,:]
          elif elec == 'gran':
            lfpPlot = lfp[g,:]

        elif isinstance(elec,Number):
          lfpPlot = lfp[elec,:] 
        fs = int(1000.0/dt) # dt is equivalent to sim.cfg.recordStep  (units of dt: ms)
        t_spec = np.linspace(0,index2ms(len(lfpPlot),fs), len(lfpPlot))
        spec.append(MorletSpec(lfpPlot, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq, lfreq=freqList))


      f = freqList if freqList is not None else np.array(range(minFreq, maxFreq+1, stepFreq))   # only used as output for user
      vmin = np.array([s.TFR for s in spec]).min()
      vmax = np.array([s.TFR for s in spec]).max()

      for i,elec in enumerate(electrodes):
        plt.subplot(np.ceil(len(electrodes) / numCols), numCols, i + 1)
        T = timeRange
        F = spec[i].f
        if normSpec:
          spec[i].TFR = spec[i].TFR / vmax
          S = spec[i].TFR
          vc = [0, 1]
        else:
          S = spec[i].TFR
          vc = [vmin, vmax]
      

        plt.imshow(S,extent=(np.amin(T), np.amax(T), np.amin(F), np.amax(F)), origin='lower', interpolation='None', aspect='auto', vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap('viridis'))
        if normSpec:
          plt.colorbar(label='Normalized power')
        else:
          plt.colorbar(label='Power')
        plt.ylabel('Hz')
        plt.tight_layout()
        if len(electrodes) > 1:
          plt.title('Electrode %s' % (str(elec)), fontsize=fontSize - 2)


        ## Add vertical lines to stimulus onset    ## trigtimes only if non-spontaneous data (click or speech)
        ax = plt.gca() 

        if triglines:
          if trigtimesMS is not None:
            for trig in trigtimesMS:
              if trig >= timeRange[0] and trig <= timeRange[1]:
                #ax.annotate(' ', xy = (trig,24), xytext = (trig,24), arrowprops = dict(facecolor='red', shrink=0.1, headwidth=4,headlength=4),annotation_clip=False)
                ax.axvline(trig,linestyle='dashed') #ymin=,ymax=,



    # Skipping --> if transformMethod == 'fft':


    plt.xlabel('time (ms)', fontsize=fontSize)
    plt.suptitle('LFP spectrogram', size=fontSize, fontweight='bold')
    plt.subplots_adjust(bottom=0.08, top=0.90)
    
    if saveFig:
      if fn is None:
        electrodeString = listToString(electrodes)
        figname = 'NHP_LFP_spect_electrodes_%s.png' % electrodeString  
      else:
        pathParts = fn.split("/")
        for i,pathPart in enumerate(pathParts):
          if '.mat' in pathPart:
            filenameMat = pathParts[i]
            filename = filenameMat[:-4]
        print('Saving LFP spectrogram for file ' + filename)
        figname = 'NHP_LFP_spect_%s.png' % filename
      try:
        plt.savefig(figname, dpi=dpi) #dpi 
      except:
        plt.savefig('NHP_LFP_spectrogram.png',dpi=dpi)


    if showFig:
      plt.show()

  ### PSD (Power Spectral Density) -----------------------------------------
  if 'PSD' in plots:
    if overlay:
      figs.append(plt.figure(figsize=figSize))
    else:
      numCols = np.round(len(electrodes) / maxPlots) + 1  # 1  # maxPlots = 8
      figs.append(plt.figure(figsize=(figSize[0]*numCols, figSize[1])))

    allFreqs = []
    allSignal = []

    for i,elec in enumerate(electrodes):
      if elec == 'avg':
        lfpPlot = np.mean(lfp, axis=0)

      elif elec in ['supra', 'infra', 'gran']:
        if dbpath is None:
          print('need path to .csv layer file')
        if fn is None:
          print('need path to .mat data file')
        else:
          s2,g,i1 = getflayers(fn,dbpath=dbpath,abbrev=True)
          print('supra: ' + str(s2))
          print('gran: ' + str(g))
          print('infra: ' + str(i1))

        if elec == 'supra':
          lfpPlot = lfp[s2,:]
        
        elif elec == 'infra':
          lfpPlot = lfp[i1,:]

        elif elec == 'gran':
          lfpPlot = lfp[g,:]

      elif isinstance(elec, Number):
        lfpPlot = lfp[elec,:]

      # Morlet wavelet transform method
      if transformMethod == 'morlet':
        from testMorlet import MorletSpec, index2ms

        Fs = int(1000.0/dt)

        #t_spec = np.linspace(0, index2ms(len(lfpPlot), Fs), len(lfpPlot))
        morletSpec = MorletSpec(lfpPlot, Fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq)
        freqs = F = morletSpec.f
        spec = morletSpec.TFR
        signal = np.mean(spec, 1)
        ylabel = 'Power'

      # FFT transform method ## THIS CONDITION IS UNTESTED 
      elif transformMethod == 'fft':
        Fs = int(1000.0/dt)
        power = mlab.psd(lfpPlot, Fs=Fs, NFFT=NFFT, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=noverlap, pad_to=None, sides='default', scale_by_freq=None)

        if smooth:
          signal = _smooth1d(10*np.log10(power[0]), smooth)
        else:
          signal = 10*np.log10(power[0])
          freqs = power[1]
          ylabel = 'Power (dB/Hz)'


      allFreqs.append(freqs)
      allSignal.append(signal)


    if normPSD:
      vmax = np.max(allSignal)
      for i, s in enumerate(allSignal):
        allSignal[i] = allSignal[i]/vmax

    for i,elec in enumerate(electrodes):
      if not overlay:
        plt.subplot(np.ceil(len(electrodes)/numCols), numCols,i+1)
      if elec == 'avg':
        color = 'k'
      elif elec == 'supra':
        color = 'r'
      elif elec == 'infra':
        color = 'b'
      elif elec == 'gran':
        color = 'g'
      elif isinstance(elec, Number):
        color = colors[i%len(colors)]
      
      freqs = allFreqs[i]
      signal = allSignal[i]
      plt.plot(freqs[freqs<maxFreq], signal[freqs<maxFreq], linewidth=lineWidth, color=color, label='Electrode %s'%(str(elec)))
      plt.xlim([0, maxFreq])
      if len(electrodes) > 1 and not overlay:
        plt.title('Electrode %s'%(str(elec)), fontsize=fontSize)
      plt.ylabel(ylabel, fontsize=fontSize)

    # format plot
    plt.xlabel('Frequency (Hz)', fontsize=fontSize)
    if overlay:
      plt.legend(fontsize=fontSize)
    plt.tight_layout()
    plt.suptitle('LFP Power Spectral Density', fontsize=fontSize, fontweight='bold') # add yaxis in opposite side
    plt.subplots_adjust(bottom=0.08, top=0.92)

    if saveFig:
      if fn is None:
        electrodeString = listToString(electrodes)
        figname = 'NHP_LFP_PSD_electrodes_%s.png' % electrodeString  
      else:
        pathParts = fn.split("/")
        for i,pathPart in enumerate(pathParts):
          if '.mat' in pathPart:
            filenameMat = pathParts[i]
            filename = filenameMat[:-4]
        print('Saving LFP PSD for file ' + filename)
        figname = 'NHP_LFP_PSD_%s.png' % filename
      try:
        plt.savefig(figname, dpi=dpi) #dpi 
      except:
        plt.savefig('NHP_LFP_PSD.png',dpi=dpi)


    if showFig:
      plt.show()

######################################
### FILE PRE-PROCESSING FUNCTIONS #### 
######################################

# Sorts data .mat files by recording region
def sortFiles(pathToData,regions):
  # pathToData -- string -- should go to parent directory with the raw unsorted .mat files
  # regions -- list of string or numbers -- either number code or name code for recording regions of interest (e.g. ['A1' 'MGB'] or [1 7])
  ## ^^ Make it so can either delete or sort the files not recorded in 'regions'

  # (1) Create a list of all the unsorted .mat files 
  ## NOTE: COMBINE THESE LINES? TEST. 
  origDataFiles = [f for f in os.listdir(pathToData) if os.path.isfile(os.path.join(pathToData,f))]
  origDataFiles = [f for f in origDataFiles if '.mat' in f] # list of all the .mat data files to be processed 
  print(origDataFiles)


  # (2) Set up dict to contain A1, MGB, and TRN filenames
  recordingAreaCodes = {1:'A1', 2:'belt', 3:'MGB', 4:'LGN', 5:'Medial Pulvinar', 6:'Pulvinar', 7:'TRN', 8:'Motor Ctx', 9:'Striatum', 10:'SC', 11:'IP', 33:'MGBv'} # All of the area codes -- recording region pairs 
  numCodes = list(recordingAreaCodes.keys()) # [1, 3, 7]
  nameCodes = list(recordingAreaCodes.values()) # ['A1', 'MGB', 'TRN']

  DataFiles = {} 

  for fn in origDataFiles:
    fullFN = pathToData + fn
    #print(fullFN)
    fp = h5py.File(fullFN,'r')
    areaCode = int(fp['params']['filedata']['area'][0][0]) # 1 - A1   # 3 - MGB   # 7 - TRN
    if areaCode in numCodes:
      area = str(recordingAreaCodes[areaCode]) # e.g. 'A1', 'MGB'
      if area not in list(DataFiles.keys()):
        DataFiles[area] = [] 
        DataFiles[area].append(fn)
      else:
        DataFiles[area].append(fn)
    else:
      print('Invalid area code in file %s' % fn)


  # (3) Move files into appropriate subdirectories
  for key in DataFiles.keys():
    for file in DataFiles[key]:
      origFilePath = pathToData + file
      newPath = pathToData + key + '/' # /A1/ etc. 
      newFilePath = newPath + file 
      if os.path.isdir(newPath):
        shutil.move(origFilePath,newFilePath)
      elif not os.path.isdir(newPath):
        os.mkdir(newPath)
        shutil.move(origFilePath,newFilePath)

  return DataFiles # Change this...? 


def moveDataFiles(pathToData,option):  ## deletes or moves irrelevant .mat files
  # pathToData -- path to parent dir with unsorted or unwanted .mat files 
  # option -- delete or move to 'other'

  # list of unsorted files 
  leftoverFiles = [q for q in os.listdir(pathToData) if os.path.isfile(os.path.join(pathToData,q))]
  leftoverFiles = [q for q in leftoverFiles if '.mat' in q]

  for left in leftoverFiles:
    fullLeft = pathToData + left             # full path to leftover .mat file 
    if option is None or option == 'delete':  # DEFAULT IS TO DELETE THE OTHER UNSORTED FILES
      if os.path.isfile(fullLeft):
        print('Deleting ' + left)             # INSTEAD OF DELETING SHOULD I JUST MOVE THE FILE? 
        os.remove(fullLeft)
    elif option == 'move': # MOVE TO 'other' DIRECTORY
      otherDir = pathToData + 'other/'
      otherFilePath = otherDir + left
      if os.path.isdir(otherDir):
        shutil.move(fullLeft,otherFilePath)
      elif not os.path.isdir(otherDir):
        os.mkdir(otherDir)
        shutil.move(fullLeft,otherFilePath)



##################################  
### CSD PROCESSING FUNCTIONS #### 
### ** needs work ** 
##################################

# What functions do I need? 
def plotExpData(pathToData,expCondition,saveFolder,regions):
  # pathToData: string -- path to parent dir containing the .mat files (or the level above that), e.g. '../data/NHPdata/click/contproc/'
  # expCondition: string -- which type of trial? --> 'click' or 'spont' or 'speech'
  # saveFolder: string -- what is the parent dir that you wish to save the .png files in?
  # regions: list of numbers or strings that indicates recording region(s) of interest, e.g. [1, 3, 7], ['A1', 'MGB', 'TRN']

  conditions = ['click', 'spont', 'speech']

  recordingAreaCodes = {1:'A1', 2:'belt', 3:'MGB', 4:'LGN', 5:'Medial Pulvinar', 6:'Pulvinar', 7:'TRN', 8:'Motor Ctx', 9:'Striatum', 10:'SC', 11:'IP', 33:'MGBv'} # All of the area codes -- recording region pairs 
  numCodes = list(recordingAreaCodes.keys())      # e.g. [1, 3, 7]
  nameCodes = list(recordingAreaCodes.values())   # e.g. ['A1', 'MGB', 'TRN']

  # Create directory to save figures 
  if expCondition in conditions: 
    pathToFigs = saveFolder + expCondition + '/' # e.g. '../data/NHPdata/CSD/click/'
    if not os.path.isdir(saveFolder):
      os.mkdir(saveFolder)
    if not os.path.isdir(pathToFigs): # Make subdirs for each condition  
      os.mkdir(pathToFigs)
    for region in regions:
      if region in numCodes:
        areaName = str(recordingAreaCodes[region])
        regionDir = pathToFigs + areaName + '/'  # e.g. '../data/NHPdata/CSD/click/A1/'
        if not os.path.isdir(regionDir):
          os.mkdir(regionDir)
      elif region not in numCodes:
        print('Recording region not recognized.')

  # Find .mat data files for each region specified in 'regions' arg 
  ## MAKE THIS BETTER! 
  for area in nameCodes:
    areaDataDir = pathToData + area + '/'
    if not os.path.isdir(areaDataDir):
      print('No region-sorted .mat files detected for ' + str(area))
    else:
      print('.mat files found for region ' + str(area))
      dataFiles = [f for f in os.listdir(areaDataDir) if os.path.isfile(os.path.join(areaDataDir,f))]
      dataFiles = [f for f in dataFiles if '.mat' in f] # list of all the .mat files in each area's data folder 
    
    for file in dataFiles:
      fullPath = areaDataDir + file # full path to data .mat file 
      [sampr,LFP_data,dt,tt,CSD_data,trigtimes] = loadfile(fn=fullPath, samprds=11*1e3, spacing_um=100)
          # sampr is the sampling rate after downsampling 
          # tt is time array (in seconds)
          # trigtimes is array of stim trigger indices
          # NOTE: make samprds and spacing_um args in this function as well for increased accessibility??? 

      #### PLOT INTERPOLATED CSD COLOR MAP PLOT #### 
      plotCSD(fn=fullPath,dat=CSD_data,tt=tt,timeRange=[1100,1200],showFig=False)
          # NOTE: make timeRange an arg in this function!!! 

      # REMOVE BAD EPOCHS FIRST..?  

      # GET AVERAGE ERP ## 
      ## set epoch params
      swindowms = 0 # start time relative to stimulus 
      ewindowms = 200 # end time of epoch relative to stimulus onset 
        # NOTE: make swindowms and ewindowms args in this function as well?!?!

      # calculate average CSD ERP 
      ttavg,avgCSD = getAvgERP(CSD_data, sampr, trigtimes, swindowms, ewindowms)
      plotAvgCSD(fn=fullPath,dat=avgCSD,tt=ttavg,showFig=False)

      # NOTE: option(s) regarding what is being plotted!!!!! 

    # LAST STEP: move ths saved .png file to the appropriate directory
    # fileName = file[0:-4] # cuts off the .mat 
    # CURRENTLY --> we should be looking / creating files in ooh... analysis. hm. 
    pngFiles = [q for q in os.listdir(areaDataDir) if os.path.isfile(os.path.join(areaDataDir,q))]
    pngFiles = [q for q in pngFiles if '.png' in f]
    


def plot_LFP_PSD(dataFile, LFP_data, electrodes):

    dur_samples = LFP_data.shape[1]
    dur_ms = int(dur_samples/(samprate/1000.))  # ms (eg 348sec = 348,000 ms)
    step = 10000
    tranges = [[t, t+step] for t in range(0, (dur_ms-step), step)]

    allData = []

    for timeRange in tranges:
        #plotLFP(dat=LFP_data, tt=tt, timeRange=[7500,8500], plots=['spectrogram'], electrodes=[2,8,13,18], maxFreq=80, saveFig=True, fn=fullPath, dbpath=dbpath) # fn=fullPath,dbpath = dbpath,  # 16,19 #[4,12]
        out = sim.analysis.plotLFP(**{
                'inputLFP': LFP_data.T,
                'plots': ['PSD'], 
                'electrodes': electrodes, 
                'timeRange': timeRange,
                'maxFreq': 50, 
                'figSize': (8,12), 
                'saveData': False, 
                'saveFig': '../data/NHPdata/spont/spont_LFP_PSD/%s_lfp_psd_%d_%d.png' % (dataFile[:-4], timeRange[0], timeRange[1]),
                'showFig': False})
        allData.append(out[1]['allSignal'])
        
    with open('../data/NHPdata/spont/spont_LFP_PSD/%s_10sec_allData.pkl' % (dataFile[:-4]), 'wb') as f:
        pkl.dump({'allData': allData, 'allFreqs': out[1]['allFreqs']}, f)    




def plot_LFP_PSD_combined(dataFile, norm=False):
    NHP = dataFile[:-4]
    
    with open('../data/NHPdata/spont/spont_LFP_PSD/%s_10sec_allData.pkl' % (NHP), 'rb') as f:
        loadedData = pkl.load(f)
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
        plt.savefig('../data/NHPdata/spont/spont_LFP_PSD/%s_10sec_PSD_combined_norm.png' % (NHP))
    else:
        plt.savefig('../data/NHPdata/spont/spont_LFP_PSD/%s_10sec_PSD_combined.png' % (NHP))




###########################
######## MAIN CODE ########
###########################

if __name__ == '__main__':

    # Parent data directory containing .mat files
    origDataDir = '../data/NHPdata/spont/'   # LOCAL DIR 

    ## Sort these files by recording region 
    # DataFiles = sortFiles(origDataDir, [1, 3, 7]) # path to data .mat files  # recording regions of interest

    ## Delete or move unwanted / unworted .mat data files 
    # moveDataFiles(origDataDir,'move')

    ############### 
    recordingArea = 'A1/' # 'MGB/' 

    test = 1 # set to 1 if testing a particular monkey, 0 if going through all files in data dir
    testFiles = ['1-bu031032017@os_eye06_20.mat', '2-ma031032023@os_eye06_20.mat', '2-rb031032016@os_eye06_20.mat', '2-rb045046026@os_eye06_20.mat']   # CHANGE FILE HERE IF LOOKING AT SPECIFIC MONKEY
    testFiles = ['2-rb031032016@os_eye06_20.mat']

    if test:
        dataFiles = testFiles
    else:
        dataPath = origDataDir + recordingArea
        dataFilesList = os.listdir(dataPath) 
        dataFiles = []
        for file in dataFilesList:
            if '.mat' in file:
                dataFiles.append(file)

    # setup netpyne
    samprate = 11*1e3  # in Hz
    sim.initialize()
    sim.cfg.recordStep = 1000./samprate # in ms


    for dataFile in dataFiles: 
        
        fullPath = origDataDir + dataFile # + recordingArea + dataFile      # Path to data file 

        [sampr,LFP_data,dt,tt,CSD_data,trigtimes] = loadfile(fn=fullPath, samprds=samprate, spacing_um=100)
                # sampr is the sampling rate after downsampling 
                # tt is time array (in seconds)
                # trigtimes is array of stim trigger indices
                # NOTE: make samprds and spacing_um args in this function as well for increased accessibility??? 

        ##### SET PATH TO .csv LAYER FILE ##### 
        dbpath = '../data/NHPdata/spont/21feb02_A1_spont_layers.csv'  # GCP # CHANGE ACCORDING TO MACHINE USED TO RUN 
        
        ##### GET LAYERS FOR OVERLAY #####
        s1low,s1high,s2low,s2high,glow,ghigh,i1low,i1high,i2low,i2high = getflayers(fullPath,dbpath=dbpath,getmid=False,abbrev=False) # fullPath is to data file, dbpath is to .csv layers file 
        lchan = {}
        lchan['S'] = s2high
        lchan['G'] = ghigh
        lchan['I'] = CSD_data.shape[0]-1 #i2high
        print('s2high: ' + str(s2high))

        


        ### Plot batches of CSDs:
        ## Set up time ranges for loop
        tranges = [[x, x+200] for x in range(2000, 3000, 200)] # bring it down to 175-250 if possible
        tranges = [[2800, 3000]]
        for t in tranges:
            plotCSD(fn=fullPath,dat=CSD_data,tt=tt,
                    trigtimes=None,timeRange=[t[0], t[1]],
                    showFig=False, figSize=(6,9), 
                    layerLines=True, layerBounds=lchan, 
                    overlay='LFP', LFP_data=LFP_data, smooth=33,
                    saveFig=dataFile[:-4]+'_CSD_%d-%d' % (t[0], t[1]))

        # -----------------------
        # LFPs
        
        #electrodes = ['avg', list(range(s1low, glow)), list(range(glow, i1low)), list(range(i1low, i2high))]
        
        #plot LFP PSDs
        #plot_LFP_PSD(dataFile, LFP_data, electrodes)
        
        #plot LFP PSD combined    
        #plot_LFP_PSD_combined(dataFile, norm=True)
