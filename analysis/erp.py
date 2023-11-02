from pylab import *
import numpy as np
import scipy.signal as sps
import os
import pandas as pd

def vtoint (vec): return [int(x) for x in vec]
def index2ms (idx, sampr): return 1e3*idx/sampr
def ms2index (ms, sampr): return int(sampr*ms/1e3)

def calPosThresh(dat, sigmathresh):
  #return dat.mean() + sigmathresh * dat.std()
  return 100

def calNegThresh(dat, sigmathresh):
  #return dat.mean() - sigmathresh * dat.std()
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
  ewindowidx = ms2index(ewindowms,sampr)

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

# get the average ERP (dat should be either LFP or CSD)
def getERPOnChan (dat, sampr, chan, trigtimes, swindowms, ewindowms):
  nrow = dat.shape[0]
  tt = np.linspace(swindowms, ewindowms,ms2index(ewindowms - swindowms,sampr))
  swindowidx = ms2index(swindowms,sampr) # could be negative
  ewindowidx = ms2index(ewindowms,sampr)
  lERP = np.zeros((len(trigtimes),len(tt)))
  for i,trigidx in enumerate(trigtimes): # go through stimuli
    sidx = max(0,trigidx+swindowidx)
    eidx = min(dat.shape[1],trigidx+ewindowidx)
    lERP[i,:] = dat[chan, sidx:eidx]
  return tt,lERP

# get the average ERP (dat should be either LFP or CSD)
def getAvgERP (dat, sampr, trigtimes, swindowms, ewindowms):
  nrow = dat.shape[0]
  tt = np.linspace(swindowms, ewindowms,ms2index(ewindowms - swindowms,sampr))
  # print('tt: ' + str(tt))
  swindowidx = ms2index(swindowms,sampr) # could be negative
  ewindowidx = ms2index(ewindowms,sampr)
  avgERP = np.zeros((nrow,len(tt)))
  for chan in range(nrow): # go through channels
    for trigidx in trigtimes: # go through stimuli
      sidx = max(0,trigidx+swindowidx)
      #print('sidx: ' + str(sidx))
      eidx = min(dat.shape[1],trigidx+ewindowidx)
      #print('eidx: ' + str(eidx))
      avgERP[chan,:] += dat[chan, sidx:eidx]
    avgERP[chan,:] /= float(len(trigtimes))
  return tt,avgERP

# draw the average ERP (dat should be either LFP or CSD)
def drawAvgERP (dat, sampr, trigtimes, swindowms, ewindowms, whichchan=None, yl=None, clr=None,lw=1):
  plt.figure(figsize=(9,4))
  ttavg,avgERP = getAvgERP(dat,sampr,trigtimes,swindowms,ewindowms)
  nrow = avgERP.shape[0]
  if isinstance(whichchan, list): # avg across channels
    avgERPChans = mean(avgERP[whichchan, :], 0)
    plot(ttavg,avgERPChans[:],color=clr,linewidth=lw)  
    xlim((swindowms,ewindowms))
    if yl is not None: ylim(yl)
  else:
    for chan in range(nrow): # go through channels
      if whichchan is None:
        subplot(nrow,1,chan+1)
        plot(ttavg,avgERP[chan,:],color=clr,linewidth=lw)
      elif chan==whichchan:
        plot(ttavg,avgERP[chan,:],color=clr,linewidth=lw)
      #xlim((-swindowms,ewindowms))
      xlim((swindowms,ewindowms))
      if yl is not None: ylim(yl)
  
# draw the event related potential (or associated CSD signal), centered around stimulus start (aligned to t=0)
def drawERP (dat, sampr, trigtimes, windowms, whichchan=None, yl=None,clr=None,lw=1):
  if clr is None: clr = 'gray'
  nrow = dat.shape[0]
  tt = np.linspace(-windowms,windowms,ms2index(windowms*2,sampr))
  windowidx = ms2index(windowms,sampr)
  for trigidx in trigtimes: # go through stimuli
    for chan in range(nrow): # go through channels
      sidx = max(0,trigidx-windowidx)
      eidx = min(dat.shape[1],trigidx+windowidx)
      if whichchan is None:
        subplot(nrow,1,chan+1)
        plot(tt,dat[chan, sidx:eidx],color=clr,linewidth=lw)
      elif chan==whichchan:
        plot(tt,dat[chan, sidx:eidx],color=clr,linewidth=lw)
      xlim((-windowms,windowms))
      if yl is not None: ylim(yl)
      #xlabel('Time (ms)')

# normalized cross-correlation between x and y
def normcorr (x, y):
  # Pad shorter array if signals are different lengths
  if x.size > y.size:
    pad_amount = x.size - y.size
    y = np.append(y, np.repeat(0, pad_amount))
  elif y.size > x.size:
    pad_amount = y.size - x.size
    x = np.append(x, np.repeat(0, pad_amount))
  corr = np.correlate(x, y, mode='full')  # scale = 'none'
  lags = np.arange(-(x.size - 1), x.size)
  corr /= np.sqrt(np.dot(x, x) * np.dot(y, y))
  return lags, corr

# x is longer signal; y is short pattern; nsamp is moving window size (in samples) for finding pattern
def windowcorr (x, y, nsamp, verbose=False):
  sz = len(x)
  lsidx,leidx=[],[]
  llag, lc = [],[]
  for sidx in range(0,sz,nsamp):
    lsidx.append(sidx)
    eidx = min(sidx + nsamp, sz-1)
    leidx.append(eidx)
    if verbose: print(sidx,eidx)
    sig = sps.detrend(x[sidx:eidx])
    lags,c = normcorr(sig,y)
    llag.append(lags[int(len(lags)/2):])
    lc.append(c[int(len(lags)/2):])
  return llag, lc, lsidx, leidx

#
def maxnormcorr (x, y):
  lags, corr = normcorr(x,y)
  return max(corr)

#
def findpeakERPtimes (sig, erp, winsz, sampr, dfctr=2, thresh=0.05):
  llag, lc, lsidx, leidx = windowcorr(sig, erp, int(winsz*sampr))
  d = int(dfctr * len(erp)) # minimum distance between peaks and troughs (in samples)
  lpkpos,lpkprop = [],[]
  lT = []
  for i,C in enumerate(lc):
    pkpos, pkprop = sps.find_peaks(C, height = thresh, threshold = None, distance=d)
    lpkpos.append(pkpos)
    lpkprop.append(pkprop)
    for t in pkpos: lT.append(index2ms(lsidx[i] + t, sampr))  
  return {'lT':lT, 'llag':llag, 'lc':lc, 'lsidx':lsidx,'leidx':leidx, 'lpkpos':lpkpos,'lpkprop':lpkprop}

# add ERP score to pdf; ddx must have average s2,g,i1 ERPs
def addERPscore (ddx, lschan, pdf):
  lchan = list(set(pdf['chan']))
  lchan.sort()
  pdf['ERPscore'] = pd.Series([-2 for i in range(len(pdf))], index=pdf.index) # -2 is invalid value
  lERPAvg = [ddx[s] for s in lschan]
  for ERPAvg,chan in zip(lERPAvg,lchan):
    s = pdf[pdf.chan==chan]
    for idx in s.index:
      sig0 = pdf.at[idx,'CSDwvf']
      pdf.loc[idx,'ERPscore'] = maxnormcorr(sig0,ERPAvg)

#      
def getAvgERPInDir (based, stimIntensity, needBBN, needCX, needThal,\
                    swindowms=0, ewindowms=150,
                    dbpath='data/nhpdat/spont/A1/19apr4_A1_spont_LayersForSam.csv',
                    useBIP=False):
  from nhpdat import getflayers, getdownsampr, getorigsampr, getTriggerIDs, closestfile, hasBBNStim,IsCortex,IsThal
  from nhpdat import getStimIntensity, getTriggerTimes
  dd = {}
  for fn in os.listdir(based):
    if not fn.endswith('.mat'): continue
    FN = os.path.join(based,fn)
    if stimIntensity > 0 and getStimIntensity(FN) != stimIntensity: continue
    if needBBN and not hasBBNStim(FN): continue
    if needCX and not IsCortex(FN): continue
    if needThal and not IsThal(FN): continue    
    s2,g,i1=-1,-1,-1; lchan = []
    if IsCortex(FN):
      s2,g,i1=getflayers(closestfile(fn,dbpath=dbpath)[0],abbrev=True)
      if s2 < 0: continue # no layer/channel information
      lchan = [s2,g,i1]
    samprds = getdownsampr(FN)    
    divby = getorigsampr(FN) / samprds
    trigtimes = [int(round(x)) for x in np.array(getTriggerTimes(FN)) / divby] 
    trigIDs = getTriggerIDs(FN)
    if useBIP:
      sampr,dat,dt,tt,CSD,MUA,BIP = loadfile(FN,samprds,getbipolar=True)
    else:
      sampr,dat,dt,tt,CSD,MUA = loadfile(FN,samprds)
    ttrigtimes = [index2ms(t,sampr) for t in trigtimes]
    if useBIP:
      ttavg,avgBIP = getAvgERP(BIP, sampr, trigtimes, swindowms, ewindowms)
      ddf = {'fn':fn,'ttavg':ttavg,'avgBIP':avgBIP,'sampr':sampr}      
    else:
      ttavg,avgCSD = getAvgERP(CSD, sampr, trigtimes, swindowms, ewindowms)
      ddf = {'fn':fn,'ttavg':ttavg,'avgCSD':avgCSD,'sampr':sampr}
    print(fn,lchan)        
    if s2 >= 0:
      if useBIP:
        s2+=1; g+=1; i1+=1;
      ddf['s2']=s2; ddf['g']=g; ddf['i1']=i1        
    else:
      th=int(CSD.shape[0]/2)
      lchan=[th]      
      if useBIP: th+=1
      ddf['th']=th        
    if useBIP:
      ddf['lchan'] = [x+1 for x in lchan]
    else:
      ddf['lchan'] = lchan
    dd[fn] = ddf
  return dd

#
def avgERPOverChan (dd, noiseth=0.75):
  from nhpdat import getflayers, getdownsampr, getorigsampr, getTriggerIDs, closestfile, hasBBNStim,IsCortex,IsThal
  from nhpdat import getStimIntensity, getTriggerTimes  
  ddx = {'s2':[],'g':[],'i1':[],'tt':None}
  for k in dd:
    if getorigsampr('data/nhpdat/bbn/'+k) != 44000.0: continue
    ddf = dd[k]
    lsc = ddf['lchan']
    for idx,c in enumerate([ddf[lsc[0]],ddf[lsc[1]],ddf[lsc[2]]]):
      if max(abs(ddf['avgCSD'][c,:])) > noiseth: continue
      ddx[lsc[idx]].append(ddf['avgCSD'][c,:])
    if ddx['tt'] is None: ddx['tt'] = ddf['ttavg']
  for c in lsc:
    ddx[c+'avg'] = mean(np.array(ddx[c]),axis=0)
    s = std(np.array(ddx[c]),axis=0)
    s /= sqrt(len(ddx[c]))
    ddx[c+'stderr'] = s
  return ddx
  
