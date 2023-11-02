
### Oscillation Event peakF calculations!! 

# IMPORTS # 

import numpy as np
#from bbox import bbox# , p2d
# import bbox 
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion  
from scipy.ndimage.filters import maximum_filter
from morlet import MorletSpec
# from evstats import *  ### from previous attempt(s)
from scipy.stats.stats import pearsonr
from filter import bandpass 
from scipy import ndimage
import pandas as pd


noiseampCSD = 200.0 / 10.0 # amplitude cutoff for CSD noise; was 200 before units fix

###################################################
######### FUNCTIONS FOR getblobsfrompeaks #########
###################################################

# finds boundaries where the image dips below the threshold, starting from x,y and moving left,right,up,down
def findbounds (img,x,y,thresh):
  ysz, xsz = img.shape
  y0 = y
  x0 = x - 1
  # look left
  while True:
    if x0 < 0:
      x0 = 0
      break
    if img[y0][x0] < thresh: break
    x0 -= 1
  left = x0
  # look right
  x0 = x + 1
  while True:
    if x0 >= xsz:
      x0 = xsz - 1
      break
    if img[y0][x0] < thresh: break
    x0 += 1
  right = x0
  # look down
  x0 = x
  y0 = y - 1
  while True:
    if y0 < 0:
      y0 = 0
      break
    if img[y0][x0] < thresh: break
    y0 -= 1  
  bottom = y0
  # look up
  x0 = x
  y0 = y + 1
  while True:
    if y0 >= ysz:
      y0 = ysz - 1
      break
    if img[y0][x0] < thresh: break      
    y0 += 1  
  top = y0
  #print('left,right,top,bottom:',left,right,top,bottom)  
  return left,right,top,bottom

# container for convenience - wavelet info from a single sample
class WaveletInfo:
  def __init__ (self,phs=0.0,idx=0,T=0.0,val=0.0):
    self.phs=phs
    self.idx=idx
    self.T=T
    self.val=val

# # 
# class evblob(bbox):
#   """ event blob class, inherits from bbox
#   """
#   NoneVal = -1e9
#   def __init__ (self):
#     self.avgpowevent=0 # avg val during event but only including suprathreshold pixels
#     self.avgpow=self.avgpoworig=0 # avg val during event bounds
#     self.cmass=evblob.NoneVal
#     self.maxval=self.maxpos=self.maxvalorig=evblob.NoneVal # max spectral amplitude value during event
#     self.minval=evblob.NoneVal # min val during event
#     self.minvalbefore=self.maxvalbefore=self.avgpowbefore=evblob.NoneVal
#     self.minvalafter=self.maxvalafter=self.avgpowafter=evblob.NoneVal
#     self.MUAbefore=self.MUA=self.MUAafter=evblob.NoneVal
#     self.arrMUAbefore=self.arrMUA=self.arrMUAafter=evblob.NoneVal # arrays (across channels) of avg MUA values before,during,after event   
#     self.slicex=self.slicey=evblob.NoneVal
#     self.minF=self.maxF=self.peakF=0 # min,max,peak frequencies
#     self.minT=self.maxT=self.peakT=0 # min,max,peak times
#     self.dur = self.Fspan = self.ncycle = self.dom = self.dombefore = self.domafter = self.Foct = 0
#     self.domevent = 0 # dom during event but only including suprathreshold pixels
#     # Foct is logarithmic frequency span
#     self.hasbefore = self.hasafter = False # whether has before,after period    
#     self.ID = -1
#     self.bbox = bbox()
#     self.windowidx = 0 # window index (from which spectrogram obtained)
#     self.offidx = 0 # offset into time-series (since taking windows)
#     self.duringnoise = 0 # during a noise window?
#     # indicates whether other events of given frequency co-occur (on same channel)
#     self.codelta = self.cotheta = self.coalpha = self.cobeta = self.cogamma = self.cohgamma = self.coother = 0
#     self.band = evblob.NoneVal
#     # correlation between CSD and MUA before,during,after event
#     self.CSDMUACorrbefore = self.CSDMUACorr = self.CSDMUACorrafter = 0.0
#     # a few waveform features:
#     # peak and trough values close to the spectral amplitude peak, used to align signals
#     self.WavePeakVal = self.WavePeakIDX = self.WaveTroughVal = self.WaveTroughIDX = 0
#     self.WaveH = self.WaveW = 0 # wave height and width
#     self.WavePeakT = self.WaveTroughT = 0
#     self.WaveletPeak = WaveletInfo() # for alignment by wavelet peak phase (0)
#     self.WaveletLeftTrough = WaveletInfo() # for alignment by wavelet trough phase (-pi)
#     self.WaveletRightTrough = WaveletInfo() # for alignment by wavelet trough phase (pi)
#     self.WaveletLeftH = self.WaveletLeftW = self.WaveletLeftSlope = 0 # wavelet-based height and width
#     self.WaveletRightH = self.WaveletRightW = self.WaveletRightSlope = 0 # wavelet-based height and width
#     self.WaveletFullW = 0 # full width right minus left trough times
#     self.WaveletFullH = 0 # height offset of the right minus left trough values
#     self.WaveletFullSlope = 0 # slope at full length of troughs (baseline) WaveletFullH/WaveletFullW
#     # lagged coherence - for looking at rhythmicity of signal before,during,after oscillatory event
#     #self.laggedCOH = self.laggedCOHbefore = self.laggedCOHafter = 0
#     self.filtsig = [] # filtered signal
#     self.lfiltpeak = []
#     self.lfilttrough = []
#     self.filtsigcor = 0.0 # correlation between filtered and raw signal
#     # self.oscqual = 0.0
#   def __str__ (self):
#     return str(self.left)+' '+str(self.right)+' '+str(self.top)+' '+str(self.bottom)+' '+str(self.avgpow)+' '+str(self.cmass)+' '+str(self.maxval)+' '+str(self.maxpos)
#   def draw (self,scalex,scaley,offidx=0,offidy=0,bbclr='white',mclr='r',linewidth=3):
#     x0,x1=scalex*(self.left+offidx),scalex*(self.right+offidx)
#     y0,y1=scaley*(self.top+offidy),scaley*(self.bottom+offidy)
#     drline(x0,x0,y0,y1,bbclr,linewidth)
#     drline(x1,x1,y0,y1,bbclr,linewidth)
#     drline(x0,x1,y0,y0,bbclr,linewidth)
#     drline(x0,x1,y1,y1,bbclr,linewidth)
#     plot([scalex*(self.maxpos[1]+offidx)],[scaley*(self.maxpos[0]+offidy)],mclr+'o',markersize=12)

# # get morlet specgrams on windows of dat time series (window size in samples = winsz)
# def getmorletwin (dat,winsz,sampr,freqmin=1.0,freqmax=100.0,freqstep=1.0, noiseamp=noiseampCSD,getphase=False,useloglfreq=False,mspecwidth=7.0):
#   lms = []
#   n,sz = len(dat),len(dat)
#   lnoise = []; lsidx = []; leidx = []
#   if useloglfreq:
#     minstep=0.1
#     loglfreq = getloglfreq(freqmin,freqmax,minstep)  
#   for sidx in range(0,sz,winsz):
#     lsidx.append(sidx)
#     eidx = sidx + winsz
#     if eidx >= sz: eidx = sz - 1
#     leidx.append(eidx)
#     print(sidx,eidx)
#     sig = dat[sidx:eidx]
#     lnoise.append(max(abs(sig)) > noiseamp)
#     if useloglfreq:
#       ms = MorletSpec(sig,sampr,freqmin=freqmin,freqmax=freqmax,freqstep=freqstep,getphase=getphase,lfreq=loglfreq,width=mspecwidth)
#     else:
#       ms = MorletSpec(sig,sampr,freqmin=freqmin,freqmax=freqmax,freqstep=freqstep,getphase=getphase,width=mspecwidth)
#     lms.append(ms)
#   print('exiting getmorletwin-- returning values')
#   return lms,lnoise,lsidx,leidx


# median normalization
def mednorm (dat,byRow=True):
  nrow,ncol = dat.shape[0],dat.shape[1]
  out = np.zeros((nrow,ncol)) # np.zeros <-- zeros 
  if byRow:
    for row in range(nrow):
      med = np.median(dat[row,:]) # np.median <-- median 
      if med != 0.0:
        out[row,:] = dat[row,:] / med
      else:
        out[row,:] = dat[row,:]
  else:
    for col in range(ncol):
      med = median(dat[:,col]) # np.median <-- median 
      if med != 0.0:
        out[:,col] = dat[:,col] / med
      else:
        out[:,col] = dat[:,col]
  return out



# sub average div std normalization
def unitnorm (dat,byRow=True):
  nrow,ncol = dat.shape[0],dat.shape[1]
  out = np.zeros((nrow,ncol)) ## np.zeros <-- zeros 
  if byRow:
    for row in range(nrow):
      avg = np.mean(dat[row,:])
      std = np.std(dat[row,:])
      out[row,:] = dat[row,:] - avg
      if std != 0.0:
        out[row,:] /= std
  else:
    for col in range(ncol):
      avg = np.mean(dat[:,col])
      std = np.std(dat[:,col])
      out[:,col] = dat[:,col] - avg
      if std != 0.0:
        out[:,col] /= std
  return out

######################################################################################################
# extract the event blobs from local maxima image (impk)
def getblobsfrompeaks (imnorm,impk,imorig,medthresh,endfctr,T,F):
  # imnorm is normalized image, lbl is label image obtained from imnorm, imorig is original un-normalized image
  # medthresh is median threshold for significant peaks
  # getblobfeatures returns features of blobs in lbl using imnorm
  lpky,lpkx = np.where(impk) # get the peak coordinates
  lblob = []
  for y,x in zip(lpky,lpkx):
    pkval = imnorm[y][x]
    thresh = max(medthresh, min(medthresh, endfctr * pkval)) # lower value threshold used to find end of event
    #thresh = min(medthresh, endfctr * pkval) # lower value threshold used to find end of event
    #thresh = max(medthresh, endfctr * pkval) # lower value threshold used to find end of event
    #thresh = max(medthresh, (medthresh + endfctr*pkval)/2.0) # threshold used to find end of event 
    left,right,top,bottom = findbounds(imnorm,x,y,thresh)
    #subimg = imnorm[bottom:top+1,left:right+1]
    #thsubimg = subimg > thresh
    #print('L,R,T,B:',left,right,top,bottom,subimg.shape,thsubimg.shape,sum(thsubimg))
    #print('sum(thsubimg)',sum(thsubimg),'amax(subimg)',amax(subimg))    
    b = evblob()
    #b.avgpoworig = ndimage.mean(imorig[bottom:top+1,left:right+1],thsubimg,[1])
    b.maxvalorig = imorig[y][x]
    #b.avgpow = ndimage.mean(subimg,thsubimg,[1])
    b.maxval = pkval
    b.minval = np.amin(imnorm[bottom:top+1,left:right+1]) ## np.amin <-- amin
    b.left = left
    b.right = right
    b.top = top
    b.bottom = bottom
    b.maxpos = (y,x)
    b.minF = F[b.bottom] # get the frequencies
    b.maxF = F[min(b.top,len(F)-1)]
    b.peakF = F[b.maxpos[0]]
    b.band = getband(b.peakF)
    b.minT = T[b.left]
    b.maxT = T[min(b.right,len(T)-1)]
    b.peakT = T[b.maxpos[1]]
    lblob.append(b)
  return lblob

#### FOR 'getband' above in getblobsfrompeaks
from collections import OrderedDict

# frequency band ranges for primate auditory system
def makedbands (useAudGamma = False):
  dbands = OrderedDict()
  gapHz = 1
  dbands['delta'] = [0.5,3.0 + gapHz]
  dbands['theta'] = [4,8 + gapHz]
  dbands['alpha'] = [9,14 + gapHz]
  dbands['beta'] = [15,28 + gapHz]
  if useAudGamma:
    dbands['gamma'] = [29,40 + gapHz] # gamma in aud system has lower max than traditional gamma (30-80 Hz)
    dbands['hgamma'] = [41,200 + gapHz] # considering high gamma anything above 40 Hz
  else:
    dbands['gamma'] = [29,80 + gapHz]
    dbands['hgamma'] = [81,200 + gapHz]    
  return dbands


dbands = makedbands()
lband = list(dbands.keys())

#
def getband (f):
  for k in ['delta','theta','alpha','beta','gamma','hgamma']:
    if f >= dbands[k][0] and f < dbands[k][1]:
      return k
  return 'unknown'


#      
def getFoct (minF, maxF):
  if maxF - minF > 0.0 and minF > 0.0: return np.log(maxF/minF)  # np.log <-- log 
  return 0.0


######################################################################################################
#### NOT NEEDED AT THE MOMENT 
######################################################################################################

# get oscillatory events
# lms is list of windowed morlet spectrograms, lmsnorm is spectrograms normalized by median in each power
# lnoise is whether the window had noise, medthresh is median threshold for significant events,
# lsidx,leidx are starting/ending indices into original time-series, csd is current source density
# on the single chan, MUA is multi-channel multiunit activity, overlapth is threshold for merging
# events when bounding boxes overlap, fctr is fraction of event amplitude to search left/right/up/down
# when terminating events
def getspecevents_norm (lms,lmsnorm,lnoise,medthresh,lsidx,leidx,csd,MUA,chan,sampr,overlapth=0.5,endfctr=0.5,getphase=False):
  llevent = []
  for windowidx,offidx,ms,msn,noise in zip(np.arange(len(lms)),lsidx,lms,lmsnorm,lnoise):   # np.arange <-- arange
    imgpk = detectpeaks(msn) # detect the 2D local maxima
    print('imgpk detected')
    lblob = getblobsfrompeaks(msn,imgpk,ms.TFR,medthresh,endfctr=endfctr,T=ms.t,F=ms.f) # cut out the blobs/events
    print('lblob gotten')
    lblobsig = [blob for blob in lblob if blob.maxval >= medthresh] # take only significant events
    #print('ndups in lblobsig 0 = ', countdups(lblobsig), 'out of ', len(lblobsig))    
    lmergeset,bmerged = getmergesets(lblobsig,overlapth,areaop=min) # determine overlapping events
    lmergedblobs = getmergedblobs(lblobsig,lmergeset,bmerged)
    #print('ndups in lmergedblobs A = ', countdups(lmergedblobs), 'out of ', len(lmergedblobs))
    lmergeset,bmerged = getmergesets(lmergedblobs,1.0,areaop=max) # gets rid of duplicates
    lmergedblobs = getmergedblobs(lmergedblobs,lmergeset,bmerged)
    #print('ndups in lmergedblobs B = ', countdups(lmergedblobs), 'out of ', len(lmergedblobs))
    # get the extra features (before/during/after with MUA,avg,etc.)
    ### getextrafeatures(lmergedblobs,ms,msn,medthresh,csd,MUA,chan,offidx,sampr,endfctr=endfctr,getphase=getphase)
    ### ^^ COMMENTING THIS OUT FOR NOW 
    ### print('extra features gotten')
    ndup = countdups(lmergedblobs)
    if ndup > 0: print('ndup in lmergedblobs = ', ndup, 'out of ', len(lmergedblobs))
    for blob in lmergedblobs: # store offsets for getting to time-series / wavelet spectrograms
      blob.windowidx = windowidx
      blob.offidx = offidx
      blob.duringnoise = noise
    llevent.append(lmergedblobs) # save merged events
    print('one iteration of getspecevents for-loop complete')
  return llevent


### THIS IS NOW IN LOAD.PY ###
def getspecevents_nonNorm (lms,lmsnorm,lnoise,medthresh,lsidx,leidx,csd,MUA,chan,sampr,overlapth=0.5,endfctr=0.5,getphase=False):
  llevent = []
  for windowidx,offidx,ms,msn,noise in zip(np.arange(len(lms)),lsidx,lms,lmsnorm,lnoise):   # np.arange <-- arange
    imgpk = detectpeaks(ms.TFR) # detect the 2D local maxima
    print('imgpk detected')
    lblob = getblobsfrompeaks(ms.TFR,imgpk,ms.TFR,medthresh,endfctr=endfctr,T=ms.t,F=ms.f) # cut out the blobs/events
    print('lblob gotten')
    lblobsig = [blob for blob in lblob if blob.maxval >= medthresh] # take only significant events
    #print('ndups in lblobsig 0 = ', countdups(lblobsig), 'out of ', len(lblobsig))    
    lmergeset,bmerged = getmergesets(lblobsig,overlapth,areaop=min) # determine overlapping events
    lmergedblobs = getmergedblobs(lblobsig,lmergeset,bmerged)
    #print('ndups in lmergedblobs A = ', countdups(lmergedblobs), 'out of ', len(lmergedblobs))
    lmergeset,bmerged = getmergesets(lmergedblobs,1.0,areaop=max) # gets rid of duplicates
    lmergedblobs = getmergedblobs(lmergedblobs,lmergeset,bmerged)
    #print('ndups in lmergedblobs B = ', countdups(lmergedblobs), 'out of ', len(lmergedblobs))
    # get the extra features (before/during/after with MUA,avg,etc.)
    ### getextrafeatures(lmergedblobs,ms,msn,medthresh,csd,MUA,chan,offidx,sampr,endfctr=endfctr,getphase=getphase)
    ### ^^ COMMENTING THIS OUT FOR NOW 
    ### print('extra features gotten')
    ndup = countdups(lmergedblobs)
    if ndup > 0: print('ndup in lmergedblobs = ', ndup, 'out of ', len(lmergedblobs))
    for blob in lmergedblobs: # store offsets for getting to time-series / wavelet spectrograms
      blob.windowidx = windowidx
      blob.offidx = offidx
      blob.duringnoise = noise
    llevent.append(lmergedblobs) # save merged events
    print('one iteration of getspecevents for-loop complete')
  return llevent


def detectpeaks (image):
  """
  Takes an image and detect the peaks usingthe local maximum filter.
  Returns a boolean mask of the peaks (i.e. 1 when
  the pixel's value is the neighborhood maximum, 0 otherwise)
  """
  # from https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array
  # define an 8-connected neighborhood
  neighborhood = generate_binary_structure(2,2)
  #apply the local maximum filter; all pixel of maximal value 
  #in their neighborhood are set to 1
  local_max = maximum_filter(image, footprint=neighborhood)==image
  #local_max is a mask that contains the peaks we are 
  #looking for, but also the background.
  #In order to isolate the peaks we must remove the background from the mask.
  #we create the mask of the background
  background = (image==0)
  #a little technicality: we must erode the background in order to 
  #successfully subtract it form local_max, otherwise a line will 
  #appear along the background border (artifact of the local maximum filter)
  eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
  #we obtain the final mask, containing only peaks, 
  #by removing the background from the local_max mask (xor operation)
  detected_peaks = local_max ^ eroded_background
  return detected_peaks

def countdups (lblob):
  # count duplicate blobs (same left,top,right,bottom)  
  lmergeset,bmerged = getmergesets(lblob,1.0,areaop=max) # determine overlapping events
  return len(lmergeset)

#  ### NOTE: .getintersction?? 
def getmergesets (lblob,prct,areaop=min):
  """ get the merged blobs (bounding boxes)
  lblob is a list of blos (input)
  prct is the threshold for fraction of overlap required to merge two blobs (boxes)
  returns a list of sets of merged blobs and a bool list of whether each original blob was merged
  """                                         
  sz = len(lblob)
  bmerged = [False for i in range(sz)]
  for i,blob in enumerate(lblob): blob.ID = i # make sure ID assigned
  lmergeset = [] # set of merged blobs (boxes)
  for i in range(sz):
    blob0 = lblob[i]
    for j in range(sz):
      if i == j: continue
      blob1 = lblob[j]
      if blob0.band != blob1.band: continue
      # enough overlap between bboxes? 
      if blob0.getintersection(blob1).area() >= prct * areaop(blob0.area(),blob1.area()):
        # merge them
        bmerged[i]=bmerged[j]=True
        found = False
        for k,mergeset in enumerate(lmergeset): # determine if either of these bboxes are in existing mergesets
          if i in mergeset or j in mergeset: # one of the bboxes in an existing mergeset?
            found = True
            if i not in mergeset: mergeset.add(i) # i not already there? add it in
            if j not in mergeset: mergeset.add(j) # j not already there? add it in
        if not found: # did not find either bbox in an existing mergeset? then create a new mergeset
          mergeset = set()
          mergeset.add(i)
          mergeset.add(j)
          lmergeset.append(mergeset)
  return lmergeset, bmerged

#
def getmergedblobs (lblob,lmergeset,bmerged):
  """ create a new list of blobs (boxes) based on lmergeset, and update the new blobs' properties
  """ 
  lblobnew = [] # list of new blobs
  for i,blob in enumerate(lblob):
    if not bmerged[i]: lblobnew.append(blob) # non-merged blobs are copied as is
  for mergeset in lmergeset: # now go through the list of mergesets and create the new blobs
    lblobtmp = [lblob[ID] for ID in mergeset]
    for i,blob in enumerate(lblobtmp):
      if i == 0:
        box = bbox(blob.left,blob.right,blob.bottom,blob.top)
        peakF = blob.peakF
        minF = blob.minF
        maxF = blob.maxF
        minT = blob.minT
        maxT = blob.maxT
        peakT = blob.peakT
        maxpos = blob.maxpos
        maxval = blob.maxval
        minval = blob.minval
      else:
        box = box.getunion(blob)
        minF = min(minF, blob.minF)
        maxF = max(maxF, blob.maxF)
        minT = min(minT, blob.minT)
        maxT = max(maxT, blob.maxT)
        if blob.maxval > maxval:
          peakF = blob.peakF
          peakT = blob.peakT
          maxpos = blob.maxpos
          maxval = blob.maxval
        if blob.minval < minval:
          minval = blob.minval
    blob.left,blob.right,blob.bottom,blob.top = box.left,box.right,box.bottom,box.top
    blob.minF,blob.maxF,blob.peakF,blob.minT,blob.maxT,blob.peakT=minF,maxF,peakF,minT,maxT,peakT
    blob.maxpos,blob.maxval = maxpos,maxval
    blob.minval = minval
    lblobnew.append(blob)
  return lblobnew



# find maximum value, associated index around sig[sidx-winsz:sidx+winsz] as long as within bounds of left,right
def findpeak (sig, sidx, left, right, winsz):
  sz = len(sig)
  maxval = sig[sidx]
  maxidx = sidx
  left=int(left); right=int(right); sidx=int(sidx); winsz=int(winsz)
  SIDX = min(right,sidx+1); EIDX = min(right,sidx+winsz+1)
  for idx in range(SIDX,EIDX,1):
    if sig[idx] > maxval:
      maxval = sig[idx]
      maxidx = idx
  SIDX = max(left,sidx-1); EIDX = max(left,sidx-winsz)
  for idx in range(SIDX,EIDX,-1):
    if sig[idx] > maxval:
      maxval = sig[idx]
      maxidx = idx
  return maxval,maxidx

# find minimum value, associated index around sig[sidx-winsz:sidx+winsz] as long as within bounds of left,right
def findtrough (sig, sidx, left, right, winsz):
  sz = len(sig)
  minval = sig[sidx]
  minidx = sidx
  left=int(left); right=int(right); sidx=int(sidx); winsz=int(winsz)
  SIDX = min(right,sidx+1); EIDX = min(right,sidx+winsz+1)
  for idx in range(SIDX,EIDX,1):
    if sig[idx] < minval:
      minval = sig[idx]
      minidx = idx
  SIDX = max(left,sidx-1); EIDX = max(left,sidx-winsz)
  for idx in range(SIDX,EIDX,-1):
    if sig[idx] < minval:
      minval = sig[idx]
      minidx = idx
  return minval,minidx

# find closest to val, associated index around sig[sidx-winsz:sidx+winsz] as long as within bounds of left,right
def findclosest (sig, sidx, left, right, winsz, val, lookleft=True, lookright=True):
  sz = len(sig)
  closeerr = abs(sig[sidx]-val)  
  closeidx = sidx
  closeval = sig[sidx]
  left=int(left); right=int(right); sidx=int(sidx); winsz=int(winsz)
  if lookright:
    SIDX = min(right,sidx+1); EIDX = min(right,sidx+winsz+1)
    for idx in range(SIDX,EIDX,1):
      err = abs(sig[idx] - val)
      if err < closeerr:
        closeerr = err
        closeval = sig[idx]
        closeidx = idx
        #print(err,closeval,closeidx)
  if lookleft:
    SIDX = max(left,sidx-1); EIDX = max(left,sidx-winsz)
    for idx in range(SIDX,EIDX,-1):
      err = abs(sig[idx] - val)
      if err < closeerr:
        closeerr = err      
        closeval = sig[idx]
        closeidx = idx
        #print(err,closeval,closeidx)      
  return closeval,closeidx

# ### NOTE:  1) mean might be np.mean or ndimage.mean?? 2) bandpass
def getextrafeatures (lblob, ms, img, medthresh, csd, MUA, chan, offidx, sampr, endfctr = 0.5, getphase = True, getfilt = True):
  # get extra features for the event blobs, including:
  # MUA before/after, min/max/avg power before/after
  # ms is the MorletSpec object (contains non-normalized TFR and PHS when getphase==True
  # img is the median normalized spectrogram image; MUA is the multiunit activity (should have same sampling rate)
  # chan is CSD channel (where events detected), note that csd is 1D while MUA is 2D (for now)
  #vec=h.Vector() # for getting the sample entropy
  mua = None
  if MUA is not None: mua = MUA[chan+1,:] # mua on same channel
  for bdx,blob in enumerate(lblob):
    # duration/frequency features
    blob.dur = blob.maxT - blob.minT # duration
    blob.Fspan = blob.maxF - blob.minF # linear frequency span
    blob.ncycle = blob.dur*blob.peakF/1e3 # number of cycles
    blob.Foct = getFoct(blob.minF,blob.maxF)
    ###
    w2 = int(blob.width() / 2.)
    left,right,bottom,top = blob.left,blob.right+1,blob.bottom,blob.top # these are indices into TFR (wavelet spectrogram)
    #print(bdx,left,right,bottom,top,offidx)
    subimg = img[bottom:top+1,left:right+1] # is right+1 correct if already inc'ed above?
    blob.avgpow = np.mean(subimg) # avg power of all pixels in event bounds  ## np.mean <-- mean  OR ndimage.mean?? 
    thresh = max(medthresh, min(medthresh, endfctr * blob.maxval)) # lower value threshold used to find end of event
    #thresh = min(medthresh, endfctr * blob.maxval) # lower value threshold used to find end of event
    #thresh = max(medthresh, endfctr * blob.maxval) # upper value threshold used to find end of event
    #thresh = max(medthresh, (medthresh + endfctr*blob.maxval)/2.0) # upper value threshold used to find end of event 
    thsubimg = subimg >= thresh  # 
    #print(bdx,endfctr,blob.maxval,thresh,subimg.shape,thsubimg.shape,left,right,bottom,top)
    #print(amax(subimg),amin(thsubimg),amax(thsubimg))
    blob.avgpowevent = ndimage.mean(subimg,thsubimg,[1])[0] # avg power of suprathreshold pixels    
    if blob.avgpow>0.0: blob.dom = float(blob.maxval/blob.avgpow) # depth of modulation (all pixels)
    if blob.avgpowevent>0.0: blob.domevent = float(blob.maxval/blob.avgpowevent) # depth of modulation (suprathreshold pixels)
    if mua is not None:
      blob.MUA = mean(mua[left+offidx:right+offidx]) # offset from spectrogram index into original MUA,CSD time-series
      blob.arrMUA = mean(MUA[:,left+offidx:right+offidx],axis=1) # avg MUA from each channel during the event
    if mua is not None and right - left > 1: blob.CSDMUACorr = pearsonr(csd[left+offidx:right+offidx],mua[left+offidx:right+offidx])[0]
    # a few waveform features    
    wvlen2 = (1e3/blob.peakF)/2 # 1/2 wavelength in milliseconds
    wvlensz2 = int(wvlen2*sampr/1e3) # 1/2 wavelength in samples
    blob.WavePeakVal,blob.WavePeakIDX = findpeak(csd, int(blob.maxpos[1])+offidx, left+offidx, right+offidx, wvlensz2)
    blob.WaveTroughVal,blob.WaveTroughIDX = findtrough(csd, int(blob.maxpos[1])+offidx, left+offidx, right+offidx, wvlensz2)
    blob.WavePeakIDX -= offidx; blob.WaveTroughIDX -= offidx; # keep indices within spectrogram image
    blob.WaveH = blob.WavePeakVal - blob.WaveTroughVal
    blob.WavePeakT = 1e3*blob.WavePeakIDX/sampr
    blob.WaveTroughT = 1e3*blob.WaveTroughIDX/sampr
    blob.WaveW = 2 * abs(blob.WavePeakT - blob.WaveTroughT) # should update to use wavelet peak (phase= 0)/trough(phase= -PI) 
    if getphase: # wavelet-based features
      freqIDX = list(ms.f).index(blob.peakF) # index into frequency array
      PHS = ms.PHS[freqIDX,:]
      blob.WaveletPeak.phs,blob.WaveletPeak.idx=findclosest(PHS,int(blob.maxpos[1]),left,right,wvlensz2,0.0)
      blob.WaveletPeak.val = csd[blob.WaveletPeak.idx+offidx] # +offidx for correct index into csd (blob.WaveletPeak.idx is into PHS)
      blob.WaveletPeak.T = 1e3*blob.WaveletPeak.idx/sampr #
      blob.WaveletLeftTrough.phs,blob.WaveletLeftTrough.idx=findclosest(PHS,int(blob.WaveletPeak.idx),left,right,wvlensz2+int(wvlensz2/2),-pi,lookleft=True,lookright=False)
      blob.WaveletLeftTrough.val = csd[blob.WaveletLeftTrough.idx+offidx]# +offidx for correct index into csd (WaveletLeftTrough.idx is into PHS)
      blob.WaveletLeftTrough.T = 1e3*blob.WaveletLeftTrough.idx/sampr #
      blob.WaveletRightTrough.phs,blob.WaveletRightTrough.idx=findclosest(PHS,int(blob.WaveletPeak.idx),left,right,wvlensz2+int(wvlensz2/2),pi,lookleft=False,lookright=True)
      blob.WaveletRightTrough.val = csd[blob.WaveletRightTrough.idx+offidx]# +offidx for correct index into csd (WaveletRightTrough.idx is into PHS)
      blob.WaveletRightTrough.T = 1e3*blob.WaveletRightTrough.idx/sampr #      
      blob.WaveletLeftH = blob.WaveletPeak.val - blob.WaveletLeftTrough.val
      blob.WaveletLeftW = blob.WaveletPeak.T - blob.WaveletLeftTrough.T
      if blob.WaveletLeftW != 0.0: blob.WaveletLeftSlope = blob.WaveletLeftH / blob.WaveletLeftW      
      blob.WaveletRightH = blob.WaveletPeak.val - blob.WaveletRightTrough.val
      blob.WaveletRightW = blob.WaveletRightTrough.T - blob.WaveletPeak.T
      if blob.WaveletRightW != 0.0: blob.WaveletRightSlope = blob.WaveletRightH / blob.WaveletRightW
      blob.WaveletFullH = blob.WaveletRightTrough.val - blob.WaveletLeftTrough.val
      blob.WaveletFullW = blob.WaveletRightTrough.T - blob.WaveletLeftTrough.T
      if blob.WaveletFullW != 0.0: blob.WaveletFullSlope = blob.WaveletFullH / blob.WaveletFullW
    if getfilt:
      padsz = int(sampr*0.2)
      x0 = left+offidx
      x1 = right+offidx-1
      x0p = max(0,x0-padsz)
      x1p = min(len(csd),x1+padsz)
      fsig = bandpass(csd[x0p:x1p], blob.minF, blob.maxF, sampr, zerophase=True)
      #print(x0,x1,x0p,x1p,len(fsig))
      blob.filtsig = fsig[x0-x0p:x0-x0p+x1-x0]
      if x1-x0>1: blob.filtsigcor = pearsonr(blob.filtsig,csd[x0:x1])[0]
    # look at values in period before event
    idx0 = max(0,blob.left - wvlensz2) #max(0,blob.left - w2)
    idx1 = blob.left
    if idx1 > idx0 + 1: # any period before?
      subimg = img[blob.bottom:blob.top+1,idx0:idx1]
      blob.minvalbefore = amin(subimg)
      blob.maxvalbefore = amax(subimg)
      blob.avgpowbefore = mean(subimg)
      if blob.avgpowbefore>0.0: blob.dombefore = float(blob.maxvalbefore/blob.avgpowbefore)
      idx0 += offidx; idx1 += offidx # offset from spectrogram index into original MUA,CSD time-series
      if mua is not None:
        blob.MUAbefore = mean(mua[idx0:idx1]) # offset from spectrogram index into original MUA,CSD time-series
        blob.arrMUAbefore = mean(MUA[:,idx0:idx1],axis=1) # avg MUA from each channel before the event
      if mua is not None and idx1-idx0>1:
        blob.CSDMUACorrbefore = pearsonr(csd[idx0:idx1],mua[idx0:idx1])[0]
      blob.hasbefore = True      
      #idx0 = max(0,blob.left - wvlensz2*2) 
      #idx1 = blob.left    
      #blob.laggedCOHbefore = lagged_coherence(csd[idx0:idx1], (blob.minF,blob.maxF), sampr)#, n_cycles=blob.ncycle) # quantifies rhythmicity
      #if isnan(blob.laggedCOHbefore): blob.laggedCOHbefore = -1.0
    else:
      blob.hasbefore = False
    # look at values in period after event
    idx0 = blob.right+1
    idx1 = min(idx0 + wvlensz2,img.shape[1]) # min(idx0 + w2,img.shape[1])
    if idx1 > idx0 + 1: # any period after?
      subimg = img[blob.bottom:blob.top+1,idx0:idx1]
      blob.minvalafter = amin(subimg)
      blob.maxvalafter = amax(subimg)
      blob.avgpowafter = mean(subimg)
      if blob.avgpowafter>0.0: blob.domafter = float(blob.maxvalafter/blob.avgpowafter)
      idx0 += offidx; idx1 += offidx # offset from spectrogram index into original MUA,CSD time-series
      if mua is not None:
        blob.MUAafter = mean(mua[idx0:idx1])
        blob.arrMUAafter = mean(MUA[:,idx0:idx1],axis=1) # avg MUA from each channel after the event
      if mua is not None and idx1-idx0>1:
        blob.CSDMUACorrafter = pearsonr(csd[idx0:idx1],mua[idx0:idx1])[0]
      blob.hasafter = True      
      #idx1 = min(idx0 + wvlensz2*2,img.shape[1]) # min(idx0 + w2,img.shape[1])
      #blob.laggedCOHafter = lagged_coherence(csd[idx0:idx1], (blob.minF,blob.maxF), sampr)#, n_cycles=blob.ncycle) # quantifies rhythmicity
      #if isnan(blob.laggedCOHafter): blob.laggedCOHafter = -1.0
    else:
      blob.hasafter = False
  getcoband(lblob) # get band of events co-occuring on same channel



## getCV2, getLV, getFF?  -> not causing any problems for now. AH -- these come from evstats.py in a1dat

#
def getIEIstatsbyBand (dat,winsz,sampr,freqmin,freqmax,freqstep,medthresh,lchan,MUA,overlapth=0.5,getphase=True,savespec=False,useDynThresh=False,threshfctr=2.0,useloglfreq=False,mspecwidth=7.0,noiseamp=noiseampCSD,endfctr=0.5,normop=mednorm):
  # get the interevent statistics split up by frequency band
  dout = {'sampr':sampr,'medthresh':medthresh,'winsz':winsz,'freqmin':freqmin,'freqmax':freqmax,'freqstep':freqstep,'overlapth':overlapth}
  dout['threshfctr'] = threshfctr; dout['useDynThresh']=useDynThresh; dout['mspecwidth'] = mspecwidth; dout['noiseamp']=noiseamp
  dout['endfctr'] = endfctr
  for chan in lchan:
    dout[chan] = doutC = {'delta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'theta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'alpha':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'beta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'gamma':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'hgamma':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'lnoise':[]}
    print('up to channel', chan,'getphase:',getphase)
    if dat.shape[0] > dat.shape[1]:
      sig = dat[:,chan] # signal (either CSD or LFP)
      lms,lnoise,lsidx,leidx = getmorletwin(dat[:,chan],int(winsz*sampr),sampr,freqmin=freqmin,freqmax=freqmax,freqstep=freqstep,getphase=getphase,useloglfreq=useloglfreq,mspecwidth=mspecwidth,noiseamp=noiseamp)
    else:
      sig = dat[chan,:] # signal (either CSD or LFP)
      lms,lnoise,lsidx,leidx = getmorletwin(dat[chan,:],int(winsz*sampr),sampr,freqmin=freqmin,freqmax=freqmax,freqstep=freqstep,getphase=getphase,useloglfreq=useloglfreq,mspecwidth=mspecwidth,noiseamp=noiseamp)
    print('completed morlet in getIEIstatsbyBand')
    if 'lsidx' not in dout: dout['lsidx'] = lsidx # save starting indices into original data array
    if 'leidx' not in dout: dout['leidx'] = leidx # save ending indices into original data array
    lmsnorm = [normop(ms.TFR) for ms in lms] # normalize wavelet specgram by median (when normop==mednorm) or unitnorm (sub avg div std)
    print('done with lmsnorm line')
    if useDynThresh: # using dynamic threshold?
      evthresh = getDynamicThresh(lmsnorm, lnoise, threshfctr, medthresh)
      print('useDynThresh=True, evthresh=',evthresh)
    else: #  otherwise use the default medthresh
      evthresh = medthresh
    print('evthresh = ' + str(evthresh))
    doutC['evthresh'] = evthresh # save the threshold used    
    print('doutC line completed')
    specsamp = lms[0].TFR.shape[1] # number of samples in spectrogram time axis
    print('specsamp line completed')
    specdur = specsamp / sampr # spectrogram duration in seconds
    print('specdur = ' + str(specdur))
    if 'specsamp' not in dout: dout['specsamp'] = specsamp
    if 'specdur' not in dout: dout['specdur'] = specdur    
    print ('2 if statements on specsamp and specsdur completed')
    llevent = getspecevents(lms,lmsnorm,lnoise,evthresh,lsidx,leidx,sig,MUA,chan,sampr,overlapth=overlapth,getphase=getphase,endfctr=endfctr) # get the spectral events
    print('completed llevent getspecevents')
    scalex = 1e3*specdur/specsamp # to scale indices to times
    print('scalex = ' + str(scalex))
    if 'scalex' not in dout: dout['scalex'] = scalex
    doutC['lnoise'] = lnoise # this is per channel - diff noise on each channel
    myt = 0
    for levent,msn,ms in zip(llevent,lmsnorm,lms):
      print(myt)
      """ do not skip noise so can look at noise event waveforms in eventviewer; can always filter out noise from dframe
      if lnoise[myt]: # skip noise
        myt+=1
        continue      
      """
      for band in dbands.keys(): # check events by band
        lband = getblobinrange(levent,dbands[band][0],dbands[band][1])
        count = len(lband)
        doutC[band]['Count'].append(count)
        doutC[band]['levent'].append(lband)
        if count > 2:
          lbandIEI = getblobIEI(lband,scalex)
          cv = getCV2(lbandIEI)
          doutC[band]['CV'].append(cv)
          doutC[band]['IEI'].append(lbandIEI)
        else:
          doutC[band]['IEI'].append([])
        if count > 3:
          lv = getLV(lbandIEI)
          doutC[band]['LV'].append(lv)
          print(band,len(lband),lv,cv)
      myt+=1
      for band in dbands.keys(): doutC[band]['FF'] = getFF(doutC[band]['Count'])
    if savespec:
      for MS,MSN in zip(lms,lmsnorm): MS.TFR = MSN # do not save lmsnorm separately, just copy it over to lms
      doutC['lms'] = lms
    else:
      del lms,lmsnorm # cleanup memory
      gc.collect()
  dout['lchan'] = lchan
  return dout


#### THIS FUNCTION IS NOW IN LOAD.PY as getIEIstatsbyBand_nonNorm ####
# WORKS WITH NON-NORMALIZED SPECTROGRAM DATA!!!! 
def getIEIstatsbyBand2 (inputData,winsz,sampr,freqmin,freqmax,freqstep,medthresh,lchan,MUA,overlapth=0.5,
  getphase=True,savespec=True,threshfctr=2.0,useloglfreq=False,mspecwidth=7.0,noiseamp=noiseampCSD,
  endfctr=0.5,normop=mednorm):
  ### inputData --> differs from 'dat' in that it already has the channel segmented out; SHOULD I CHANGE THIS AT SOME POINT? 
  # get the interevent statistics split up by frequency band
  dout = {'sampr':sampr,'medthresh':medthresh,'winsz':winsz,'freqmin':freqmin,'freqmax':freqmax,'freqstep':freqstep,'overlapth':overlapth}
  dout['threshfctr'] = threshfctr; dout['mspecwidth'] = mspecwidth; dout['noiseamp']=noiseamp
  dout['endfctr'] = endfctr

  # for chan in lchan:
  #   dout[chan] = doutC = {'delta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
  #                         'theta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
  #                         'alpha':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
  #                         'beta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
  #                         'gamma':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
  #                         'hgamma':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
  #                         'lnoise':[]}
  #   print('up to channel', chan,'getphase:',getphase)
  #   if dat.shape[0] > dat.shape[1]:
  #     sig = dat[:,chan] # signal (either CSD or LFP)
  #     lms,lnoise,lsidx,leidx = getmorletwin(dat[:,chan],int(winsz*sampr),sampr,freqmin=freqmin,freqmax=freqmax,freqstep=freqstep,getphase=getphase,useloglfreq=useloglfreq,mspecwidth=mspecwidth,noiseamp=noiseamp)
  #   else:
  #     sig = dat[chan,:] # signal (either CSD or LFP)
  #     lms,lnoise,lsidx,leidx = getmorletwin(dat[chan,:],int(winsz*sampr),sampr,freqmin=freqmin,freqmax=freqmax,freqstep=freqstep,getphase=getphase,useloglfreq=useloglfreq,mspecwidth=mspecwidth,noiseamp=noiseamp)
  #   print('completed morlet in getIEIstatsbyBand')

  for chan in lchan:
    dout[chan] = doutC = {'delta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'theta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'alpha':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'beta':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'gamma':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'hgamma':{'LV':[],'CV':[],'Count':[],'FF':None,'levent':[],'IEI':[]},
                          'lnoise':[]}

    # Get morlet specgrams on windows of dat time series (window size in samples = winsz)
    lms,lnoise,lsidx,leidx = getmorletwin(inputData,int(winsz*sampr),sampr,freqmin=freqmin,
      freqmax=freqmax,freqstep=freqstep,getphase=getphase,useloglfreq=useloglfreq,mspecwidth=mspecwidth,
      noiseamp=noiseamp)

    if 'lsidx' not in dout: dout['lsidx'] = lsidx # save starting indices into original data array
    if 'leidx' not in dout: dout['leidx'] = leidx # save ending indices into original data array
    
    lmsnorm = [normop(ms.TFR) for ms in lms] # normalize wavelet specgram by median (when normop==mednorm) or unitnorm (sub avg div std)
    print('done with lmsnorm line')


    ## evthresh 
    evthresh = medthresh
    print('evthresh = ' + str(evthresh))
    doutC['evthresh'] = evthresh # save the threshold used  

    print('doutC line completed')

    specsamp = lms[0].TFR.shape[1] # number of samples in spectrogram time axis
    print('specsamp line completed')
    specdur = specsamp / sampr # spectrogram duration in seconds
    print('specdur = ' + str(specdur))
    if 'specsamp' not in dout: dout['specsamp'] = specsamp
    if 'specdur' not in dout: dout['specdur'] = specdur    
    print ('2 if statements on specsamp and specsdur completed')

    ### GET SPEC EVENTS 
    # for normalized spectrogram   ### sig --> inputData 
    llevent_norm = getspecevents_norm(lms,lmsnorm,lnoise,evthresh,lsidx,leidx,inputData,MUA,chan,sampr,overlapth=overlapth,getphase=getphase,endfctr=endfctr) # get the spectral events
    print('completed llevent_norm getspecevents')

    # for non-normalized spectrogram  ### sig --> inputData 
    llevent = getspecevents_nonNorm(lms,lmsnorm,lnoise,evthresh,lsidx,leidx,inputData,MUA,chan,sampr,overlapth=overlapth,getphase=getphase,endfctr=endfctr)
    print('completed llevent getspecevents')

    ### SCALEX 
    scalex = 1e3*specdur/specsamp # to scale indices to times
    print('scalex = ' + str(scalex))
    if 'scalex' not in dout: dout['scalex'] = scalex

    ### NOISE 
    doutC['lnoise'] = lnoise # this is per channel - diff noise on each channel
    myt = 0
    for levent,msn,ms in zip(llevent,lmsnorm,lms):
      print(myt)
      """ do not skip noise so can look at noise event waveforms in eventviewer; can always filter out noise from dframe
      if lnoise[myt]: # skip noise
        myt+=1
        continue      
      """
      for band in dbands.keys(): # check events by band
        lband = getblobinrange(levent,dbands[band][0],dbands[band][1])
        count = len(lband)
        doutC[band]['Count'].append(count)
        doutC[band]['levent'].append(lband)
    #     if count > 2:
    #       lbandIEI = getblobIEI(lband,scalex)
    #       cv = getCV2(lbandIEI)
    #       doutC[band]['CV'].append(cv)
    #       doutC[band]['IEI'].append(lbandIEI)
    #     else:
    #       doutC[band]['IEI'].append([])
    #     if count > 3:
    #       lv = getLV(lbandIEI)
    #       doutC[band]['LV'].append(lv)
    #       print(band,len(lband),lv,cv)
    #   myt+=1
    #   for band in dbands.keys(): doutC[band]['FF'] = getFF(doutC[band]['Count'])
    if savespec:
      for MS,MSN in zip(lms,lmsnorm): MS.TFR = MSN # do not save lmsnorm separately, just copy it over to lms
      doutC['lms'] = lms
    else:
      del lms,lmsnorm # cleanup memory
      gc.collect()
  dout['lchan'] = lchan
  return dout


#  getcyclekeys(), getalignoffset, index2ms, getcyclefeatures, addOSCscore
def GetDFrame (dout,sampr,CSD, MUA, alignby = 'bywaveletpeak',haveMUA=True):
  totsize = 0 # total number of events
  for chan in dout['lchan']:
    for band in lband:
      for levent in dout[chan][band]['levent']:
        totsize += len(levent)
  print('total number of events = ' + str(totsize))
  row_list = []
  columns=['chan','dur','maxvalbefore','maxval','maxvalafter','ncycle',\
           'dom','dombefore','domafter','domevent',\
           'MUA','MUAbefore','MUAafter','avgpowbefore','avgpow','avgpowafter','avgpowevent',\
           'minvalbefore','minval','minvalafter','hasbefore','hasafter',\
           'band','windowidx','offidx','duringnoise',\
           'minF','maxF','peakF','Fspan','Foct',\
           'minT','maxT','peakT','left','right','bottom','top','maxpos',\
           'codelta','cotheta','coalpha','cobeta','cogamma','cohgamma','coother',\
           'CSDMUACorr','CSDMUACorrbefore','CSDMUACorrafter',\
           'WavePeakVal','WavePeakIDX','WaveTroughVal','WaveTroughIDX','WaveH','WaveW','WavePeakT','WaveTroughT',\
           'WaveletPeakPhase','WaveletPeakVal','WaveletPeakIDX','WaveletPeakT',\
           'WaveletLeftTroughPhase','WaveletLeftTroughVal','WaveletLeftTroughIDX','WaveletLeftTroughT',\
           'WaveletRightTroughPhase','WaveletRightTroughVal','WaveletRightTroughIDX','WaveletRightTroughT',\
           'WaveletLeftH','WaveletLeftW','WaveletLeftSlope',\
           'WaveletRightH','WaveletRightW','WaveletRightSlope',\
           'WaveletFullH','WaveletFullW','WaveletFullSlope',\
           'absPeakT',\
           'absminT',\
           'absmaxT',\
           'absWaveletLeftTroughT',\
           'absWaveletRightTroughT',\
           'absWaveletPeakT',\
           'filtsig','filtsigcor',\
           'MUARatDOB','MUARatDOA','arrMUAbefore','arrMUA','arrMUAafter',\
           'RLWidthRat','RLHeightRat','RLSlopeRat',
           'CSDwvf','MUAwvf','alignoffset','siglen']
  lcyckeys = getcyclekeys()
  print('lcyckeys completed')
  for k in lcyckeys: columns.append('cyc_'+k)
  print('columns.append completed')
  allevents = []
  dchanevents = {}
  MUAwvf = None # MUA waveform (if MUA provided)
  for chan in dout['lchan']:
    for band in lband:
      for levent in dout[chan][band]['levent']:
        print('levent loop starting')
        print('band: ' + str(band)) 
        print('chan: ' + str(chan)) #+ ', levent: ' + print(levent))
        for ev in levent:
          # more featurs - ratio of mua, during event over mua before,after
          MUARatDOB = MUARatDOA = arrMUAbefore = arrMUA = arrMUAafter = 0
          # right div by left width, height, slope ratios
          RLWidthRat = RLHeightRat = RLSlopeRat = 0
          #print('RLWidthRat = ' + str(RLWidthRat))
          if haveMUA:
            if ev.hasbefore and ev.hasafter:
              if ev.MUAbefore > 0: MUARatDOB = ev.MUA / ev.MUAbefore
              if ev.MUAafter > 0: MUARatDOA = ev.MUA / ev.MUAafter
            arrMUAbefore = ev.arrMUAbefore
            arrMUA = ev.arrMUA
            arrMUAafter = ev.arrMUAafter
          # ratio of wavelet left,right widths, heights, slopes
          if ev.WaveletLeftW > 0.: RLWidthRat = ev.WaveletRightW / ev.WaveletLeftW
          if ev.WaveletLeftH > 0.: RLHeightRat = ev.WaveletRightH / ev.WaveletLeftH
          if ev.WaveletLeftSlope != 0.: RLSlopeRat = ev.WaveletRightSlope / ev.WaveletLeftSlope
          ######################################################################################### 
          # get the waveforms for storage in the dataframe
          left,right = int(ev.left+ev.offidx), int(ev.right+ev.offidx) # offidx to get back into the original time-series
          alignoffset = getalignoffset(ev, alignby)  
          #print('align offset completed')        
          CSDwvf = CSD[chan,left:right] # CSD waveform
          if haveMUA: MUAwvf = MUA[chan+1,left:right] # MUA waveform
          siglen = len(CSDwvf) # signal length
          #print('siglen completed')
          #########################################################################################
          # vals is a list of values for each event
          vals = [chan,ev.dur,ev.maxvalbefore,ev.maxval,ev.maxvalafter,ev.ncycle,\
                  ev.dom,ev.dombefore,ev.domafter,ev.domevent,\
                  ev.MUA,ev.MUAbefore,ev.MUAafter,ev.avgpowbefore,ev.avgpow,ev.avgpowafter,ev.avgpowevent,\
                  ev.minvalbefore,ev.minval,ev.minvalafter,int(ev.hasbefore),int(ev.hasafter),\
                  band,ev.windowidx,ev.offidx,ev.duringnoise,\
                  ev.minF,ev.maxF,ev.peakF,ev.Fspan,ev.Foct,\
                  ev.minT,ev.maxT,ev.peakT,ev.left,ev.right,ev.bottom,ev.top,ev.maxpos[1],\
                  ev.codelta,ev.cotheta,ev.coalpha,ev.cobeta,ev.cogamma,ev.cohgamma,ev.coother,\
                  ev.CSDMUACorr,ev.CSDMUACorrbefore,ev.CSDMUACorrafter,\
                  ev.WavePeakVal,ev.WavePeakIDX,ev.WaveTroughVal,ev.WaveTroughIDX,ev.WaveH,ev.WaveW,ev.WavePeakT,ev.WaveTroughT,\
                  ev.WaveletPeak.phs,ev.WaveletPeak.val,ev.WaveletPeak.idx,ev.WaveletPeak.T,\
                  ev.WaveletLeftTrough.phs,ev.WaveletLeftTrough.val,ev.WaveletLeftTrough.idx,ev.WaveletLeftTrough.T,\
                  ev.WaveletRightTrough.phs,ev.WaveletRightTrough.val,ev.WaveletRightTrough.idx,ev.WaveletRightTrough.T,\
                  ev.WaveletLeftH,ev.WaveletLeftW,ev.WaveletLeftSlope,\
                  ev.WaveletRightH,ev.WaveletRightW,ev.WaveletRightSlope,\
                  ev.WaveletFullH,ev.WaveletFullW,ev.WaveletFullSlope,\
                  ev.peakT+index2ms(ev.offidx,sampr),\
                  ev.minT+index2ms(ev.offidx,sampr),\
                  ev.maxT+index2ms(ev.offidx,sampr),\
                  ev.WaveletLeftTrough.T+index2ms(ev.offidx,sampr),\
                  ev.WaveletRightTrough.T+index2ms(ev.offidx,sampr),\
                  ev.WaveletPeak.T+index2ms(ev.offidx,sampr),\
                  ev.filtsig,ev.filtsigcor,\
                  MUARatDOB,MUARatDOA,arrMUAbefore,arrMUA,arrMUAafter,\
                  RLWidthRat,RLHeightRat,RLSlopeRat,
                  CSDwvf, MUAwvf, alignoffset, siglen]
          ######################################################################################### 
          # get the cycle features for storage in the dataframe
          dprop = getcyclefeatures(ev.filtsig, sampr, 1.5 * ev.maxF)
          #print('dprop completed')
          for k in lcyckeys: vals.append(dprop[k])
          ######################################################################################### 
          # based on https://stackoverflow.com/questions/10715965/add-one-row-to-pandas-dataframe
          row_list.append(dict((c,v) for c,v in zip(columns,vals)))
          allevents.append(ev)
  # now create the final dataframe
  pdf = pd.DataFrame(row_list, index=np.arange(0,totsize), columns=columns)
  pdf = pdf.sort_values('absPeakT') # sort by absPeakT; index will be out of order, but will correspond to dout order
  addOSCscore(pdf) # add oscillation score
  return pdf        


# get interevent interval distribution
def getblobIEI (lblob,scalex=1.0):
  liei = []
  newlist = sorted(lblob, key=lambda x: x.left)
  for i in range(1,len(newlist),1):
    liei.append((newlist[i].left-newlist[i-1].right)*scalex)
  return liei

# get event blobs in (inclusive for lower bound, strictly less than for upper bound) range of minf,maxf
def getblobinrange (lblobf, minF,maxF): return [blob for blob in lblobf if blob.peakF >= minF and blob.peakF < maxF]


