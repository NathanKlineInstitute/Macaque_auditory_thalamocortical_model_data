
import pandas as pd
import pickle
import numpy as np
import matplotlib.pyplot as plt
from morlet import MorletSpec
import os 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker  ## for colorbar 
from pylab import *

import IPython as ipy

fontsize = 16 # fig 6
fontsize = 10 # fig 7

def readChannelFiles(dataPath, subjectName, frequencyBand, chanNumber):
	# dataPath: str; to base directory with osc event files 
	# subjectName: str; Filename of monkey or sim (e.g. 1-bu031032017@os_eye06_20, A1_v34_batch27_v34_batch27_0_0)
	# frequencyBand: 'str; delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma'
	# chanNumber: int; relevant channel number

	### dfs data frame (.pkl)
	dfsFileName = subjectName + '_chan' + str(chanNumber) + '_' + frequencyBand + '.pkl' 
	dfsFullPath = dataPath + '/' + dfsFileName
	dfs = pd.read_pickle(dfsFullPath)

	return dfs 

def readSubjectFiles(dataPath, subjectName, dlms=True):
	### df data frame (.pkl):
	dfFullPath = dataPath + subjectName + '_df.pkl'
	df = pd.read_pickle(dfFullPath)

	### dlms (.pkl): 
	if dlms:
		dlmsFullPath = dataPath + subjectName + '_dlms.pkl'
		dlms_file = open(dlmsFullPath, 'rb')
		dlms = pickle.load(dlms_file)
		dlms_file.close()


	## read allDataDict as well (contains CSD, dt, tt, sampr, dat, timeRange)
	allDataFullPath = dataPath + subjectName + '_allData.pkl'
	allData_file = open(allDataFullPath, 'rb')
	allData = pickle.load(allData_file)
	allData_file.close()

	CSD = allData['CSD']
	dt = allData['dt']
	tt = allData['tt']
	sampr = allData['sampr']
	dat = allData['dat']
	timeRange = allData['timeRange']

	if dlms:
		return df, dlms, allData, CSD, dt, tt, sampr, dat, timeRange
	else:
		return df, allData, CSD, dt, tt, sampr, dat, timeRange

# get all major event properties, used for drawing the event or other...
def geteventprop (dframe,evidx,align):
	evidx=int(evidx)
	dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,filtsigcor,Foct = [dframe.at[evidx,c] for c in ['dur','chan','hasbefore','hasafter','windowidx','offidx','left','right','minT','maxT','peakT','minF','maxF','peakF','avgpowevent','ncycle','WavePeakT','WaveTroughT','WaveletPeakT','WaveletLeftTroughT','WaveletRightTroughT','filtsigcor','Foct']]
	if 'cyc_npeak' in dframe.columns:
		cycnpeak = dframe.at[evidx,'cyc_npeak']
	else:
		cycnpeak = -1
	if 'ERPscore' in dframe.columns:
		ERPscore = dframe.at[evidx,'ERPscore']
	else:
		ERPscore = -2
	if False and 'OSCscore' in dframe.columns:
		OSCscore = dframe.at[evidx,'OSCscore']
	else:
		OSCscore = -2
	band=dframe.at[evidx,'band']
	w2 = int((right-left+1)/2.)
	left=int(left+offidx); right=int(right+offidx);
	alignoffset = 0 # offset to align waveforms to 0, only used when specified as below
	if align == 'byspecpeak':
		alignoffset = -peakT
	elif align == 'bywavepeak':
		alignoffset = -WavePeakT
	elif align == 'bywavetrough':
		alignoffset = -WaveTroughT
	elif align == 'bywaveletpeak':
		alignoffset = -WaveletPeakT
	elif align == 'bywaveletlefttrough':
		alignoffset = -WaveletLeftTroughT
	elif align == 'bywaveletrighttrough':
		alignoffset = -WaveletRightTroughT
	#print('align:',peakT,align,alignoffset)
	return dur,int(chan),hasbefore,hasafter,int(windowidx),offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT ,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore

def getStats(df, evidx,align='bywaveletpeak',verbose=False):
	dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore = geteventprop(df, evidx, align)   #= self.getallprop(evidx,align) 

	if verbose:
		return dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore
	else:
		return dur,peakF,ncycle,band

# draw a line
def drline (x0,x1,y0,y1,clr,w,ax=None):
	if ax is None: ax=gca()
	ax.plot([x0,x1],[y0,y1],clr,linewidth=w)

# draw a box
def drbox (x0,x1,y0,y1,clr,w,ax=None):
	drline(x0,x0,y0,y1,clr,w,ax)
	drline(x1,x1,y0,y1,clr,w,ax)
	drline(x0,x1,y0,y0,clr,w,ax)
	drline(x0,x1,y1,y1,clr,w,ax)

# viewer for oscillatory events
class eventviewer():
	def __init__ (self,dframe,CSD,MUA,tt,sampr,winsz,dlms,useloglfreq=False):
		self.fig = plt.figure()
		self.dframe = dframe
		self.CSD=CSD
		self.MUA=MUA
		self.tt = tt
		self.sampr=sampr
		self.winsz=winsz
		self.dlms=dlms
		self.avgCSD = self.avgMUA = self.avgSPEC = None
		self.ntrial = 0
		if self.MUA is not None:
			self.nrow = 3
		else:
			self.nrow = 2
		self.useloglfreq = useloglfreq
		self.setupax()
		self.specrange = None
	def setupax (self):
		# setup axes
		self.lax = [self.fig.add_subplot(self.nrow,1,i+1) for i in range(self.nrow)]
		self.lax[-1].set_xlabel('Time (ms)',fontsize=fontsize)#12)
		self.lax[-1].tick_params(labelsize=fontsize)	#10)
		if self.MUA is not None:
			self.lax[-1].set_ylabel('MUA') # MUA is filtered, rectified version of the LFP, so its units should be mV?
			self.lax[-2].set_ylabel(r'CSD ($mV/mm^2$)',fontsize=fontsize)
		else:
			self.lax[-1].set_ylabel('CSD Amplitude (' + r'$\frac{mV}{mm^2}$' + ')', fontsize=fontsize) #(r'CSD ($mV/mm^2$)',fontsize=12)
		self.lax[0].set_ylabel('Frequency test (Hz)',fontsize=fontsize);		#12);
		self.lax[0].tick_params(labelsize=fontsize)#10)
		
	def clf (self):
		# clear figure
		self.fig.clf()
		self.setupax()
	def clear (self): self.clf()
	def draw (self, evidx, align='bywaveletpeak', ylspec=None, clr=None, lw=1, drawfilt=True, filtclr='b', lwfilt=3, lwbox=3, verbose=True):
		# draw an event
		dframe,CSD,MUA,sampr,winsz,dlms,fig = self.dframe,self.CSD,self.MUA,self.sampr,self.winsz,self.dlms,self.fig
		gdx = 0
		lclr = ['r','g','b','c','m','y','k']
		evidx=int(evidx)
		if clr is None: clr = lclr[evidx%len(lclr)]
		dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore = geteventprop(dframe, evidx, align)	# self.getallprop(evidx,align)  
		print('windowidx: ' + str(windowidx)) ## testing dlms 
		### MAKING w2 SMALLER FOR LATER ON ###
		w2 = int(w2*0.6) 
		ax = self.lax[gdx]
		MS = dlms[chan][windowidx] # used [0] for delta since no dlms 
		if self.specrange is not None:
			vmin,vmax=self.specrange
		else:
			vmin,vmax=amin(MS.TFR),amax(MS.TFR)
		if self.useloglfreq:
			global lfreq
			img=ax.imshow(MS.TFR,extent=(MS.t[0]+alignoffset,MS.t[-1]+alignoffset,0,len(lfreq)-1),origin='lower',interpolation='None',aspect='auto',cmap=plt.get_cmap('jet'),vmin=vmin,vmax=vmax);
		else:
			img=ax.imshow(MS.TFR,extent=(MS.t[0]+alignoffset,MS.t[-1]+alignoffset,MS.f[0],MS.f[-1]),origin='lower',interpolation='None',aspect='auto',cmap=plt.get_cmap('jet'),vmin=vmin,vmax=vmax);
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='3%', pad=0.2)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
		#fmt.set_powerlimits((0,0))
		cbar=plt.colorbar(img, cax = cax, orientation='vertical', format=fmt)#, ticks=[0,5,10,15,20])		#label='Power'
		#cbar.set_label('Power (' + r'$\frac{mVˆ2}{mm^4}$' + ', size=fontsize)
		cbar.set_label('Power', size=fontsize) 
		cbar.ax.tick_params(labelsize=fontsize)
		# cbar.ax.set_yticklabels([0,5,10,15,20])  ## CUSTOM!!! 
		
		drbox(minT+alignoffset,maxT+alignoffset,minF,maxF,'r',lwbox,ax)  
		if ylspec is not None: ax.set_ylim(ylspec)
		axtstr = 'Peak Freq: ' + str(round(peakF,2)) + ' Hz'   #, channel: ' + str(chan)

		# print(axtstr) # print the info
		if verbose: ax.set_title(axtstr, fontsize=fontsize)
		gdx += 1
		#####################################################      
		#                     PLOT BEFORE
		if hasbefore:
			idx0 = max(0,left - w2)
			idx1 = left    
			# print('idx1: ' + str(idx1))  ## Adding these lines to strategize about tightening x axis 
			sig = CSD[chan,idx0:idx1]
			beforeT = (maxT-minT) * (idx1 - idx0) / (right - left + 1)
			# print('beforeT: ' + str(beforeT)) ## Adding these lines to strategize about tightening x axis 
			tt = np.linspace(minT-beforeT,minT,len(sig)) + alignoffset
			ax = self.lax[gdx+0]
			ax.plot(tt,sig,'k',linewidth=lw)
			if MUA is not None:
				ax = self.lax[gdx+1]
				ax.plot(tt,MUA[chan+1,idx0:idx1],'k',linewidth=lw)
			# print(tt[0],tt[-1])

		#####################################################      
		#                     PLOT DURING
		sig = CSD[chan,left:right]
		tt = np.linspace(minT,maxT,len(sig)) + alignoffset
		ax = self.lax[gdx+0]
		axtstr = 'Duration: '+str(round(dur,1))+' ms, '
		# if cycnpeak > -1: axtstr += str(int(cycnpeak)) + ' peaks, '
		axtstr += str(round(ncycle,1))+' cycles' #, '
		#axtstr += 'filtsigcor:'+str(round(filtsigcor,2))
		if ERPscore > -2: axtstr += ', ERPscore:'+str(round(ERPscore,2))
		if OSCscore > -2: axtstr += ', OSCscore:'+str(round(OSCscore,2))
		# print(axtstr) # print the info
		if verbose: ax.set_title(axtstr, fontsize=fontsize)#14)
		ax.plot(tt,CSD[chan,left:right],clr,linewidth=lw)
		if drawfilt:
			fsig = np.array(dframe.at[evidx,'filtsig'])
			offY = np.mean(CSD[chan,left:right]) - np.mean(fsig)
			ax.plot(tt,fsig+offY,filtclr,linewidth=lwfilt)
		if MUA is not None:
			ax = self.lax[gdx+1]
			ax.plot(tt,MUA[chan+1,left:right],clr,linewidth=lw)

		#####################################################
		#                     PLOT AFTER
		if hasafter:
			idx0 = int(right)
			idx1 = min(idx0 + w2,max(CSD.shape[0],CSD.shape[1]))
			sig = CSD[chan,idx0:idx1]
			afterT = (maxT-minT) * (idx1 - idx0) / (right - left + 1)
			tt = np.linspace(maxT,maxT+afterT,len(sig)) + alignoffset # tt = linspace(maxT,maxT+beforeT,len(sig)) + alignoffset 
			ax = self.lax[gdx+0]
			ax.plot(tt,sig,'k',linewidth=lw)
			if MUA is not None:
				ax = self.lax[gdx+1]
				ax.plot(tt,MUA[chan+1,idx0:idx1],'k',linewidth=lw)
			# print(tt[0],tt[-1])
		### reset xlim on top & bottom panels 
		# xl = ax.get_xlim()
		if hasbefore and hasafter:
			xl = (minT-beforeT + alignoffset, maxT+beforeT + alignoffset) # xl = (minT-beforeT + alignoffset, maxT+afterT + alignoffset) 
		else:
			print('NO HASBEFORE OR HASAFTER')
			xl = ax.get_xlim() # xl = (minT + alignoffset, maxT + alignoffset) 
		### reset xlim on spectrogram 
		ax = self.lax[0]
		ax.set_xlim((xl))
		### reset xlim on bottom panel (CSD waveform)
		ax1 = self.lax[1]
		ax1.set_xlim((xl))
		####### These lines will insert an invisible colorbar on the time series plot in order to match the size of the spectrogram
		divider1=make_axes_locatable(ax1)
		cax1=divider1.append_axes('right',size='3%', pad=0.2)
		cax1.axis('off')

		####
		plt.tight_layout()

def plotWavelets(dfs, df, dat, tt, sampr, dlms, subjectName, chanNumber, frequencyBand, eventIdx=None, specrange=None, ylspec=None, saveFig=0):
	# dfs; pandas dataframe read by readWaveletFile
	# dat; CSD data from getData
	# tt; time array, from getData
	# sampr; sampling rate from getData
	# dlms; dictionary, from readWaveletFile
	# df; pandas dataframe read by readWaveletFile
	# specrange; list --> for color contrasts etc 
	# ylspec; list --> spectrogram range for top panel (e.g. (1,50))
	# subjectName; str 
	# chanNumber; int
	# frequencyBand; str 

	winsz = 10

	if specrange is None:
		specrange = (0,30) # (30,80)

	if ylspec is None:
		if frequencyBand == 'gamma':
			ylspec = (30,85)
		else:
			ylspec = (1,50) # (1,10) # (30,85)


	if eventIdx is None:  
		if len(dfs) > 0:
			for idx in range(len(dfs.index)):
				evv = eventviewer(df,dat,None,tt,sampr,winsz,dlms) # event viewer to examine the output
				evv.specrange = specrange  
				evv.draw(dfs.index[idx],clr='red',drawfilt=True,filtclr='b',ylspec=ylspec) 

			print('dfs has no wavelets stored!')
	else:  			### PLOT AN INDIVIDUAL WAVELET 
		evv = eventviewer(df,dat,None,tt,sampr,winsz,dlms)
		evv.specrange = specrange  
		evv.draw(eventIdx,clr='red',drawfilt=True,filtclr='b',ylspec=ylspec)
		## SAVE FIGURE 
		if saveFig:
			figname = 'figs/fig6_' + frequencyBand + '_' + subjectName #+ '.png'
			savefig(figname + '.png', dpi=600)
		else:
			plt.show()


def plotOscEvents(oscEventsInfo, dataPaths, frequencyBands, eventTypes=['sim', 'nhp'], saveFig=0):
	### oscEventsInfo: dict w/ info about the oscillation events to plot
	### dataPaths: dict w/ paths to sim & nhp fig6_data
	### frequencyBands: list w/ desired frequency bands 
	### eventTypes: list w/ desired eventType(s) to plot 

	# eventTypes = ['sim', 'nhp']

	for band in frequencyBands:
		for eventType in eventTypes:
			## Check for dlms, df, and allData files
			dataPath = dataPaths[eventType]
			filenamePrefix = oscEventsInfo[band][eventType]['subjName']
			dlmsFile = filenamePrefix + '_dlms.pkl'
			dfFile = filenamePrefix + '_df.pkl'
			allDataFile = filenamePrefix + '_allData.pkl'
			allFiles = os.listdir(dataPath)
			if dlmsFile not in allFiles or dfFile not in allFiles or allDataFile not in allFiles:
				subjectName = oscEventsInfo[band][eventType]['subjName']
				if dlmsFile not in allFiles:
					print(band + ' ' + eventType + ' requires ' + subjectName + '_dlms.pkl file')
				if dfFile not in allFiles:
					print(band + ' ' + eventType + ' requires ' + subjectName + '_df.pkl file')
				if allDataFile not in allFiles:
					print(band + ' ' + eventType + ' requires ' + subjectName + '_allData.pkl file')
				continue
			else:
				## Plot osc event 
				print('Plotting ' + band + ' ' + eventType + ' oscillation event')
				### get data 
				subjectName = oscEventsInfo[band][eventType]['subjName']
				chan = oscEventsInfo[band][eventType]['chan']
				dfs = readChannelFiles(dataPath, subjectName, band, chan)
				# df, dlms, allData, CSD, dt, tt, sampr, dat, timeRange = readSubjectFiles(simDataPath, subjectName, dlms=False)
				df, dlms, allData, CSD, dt, tt, sampr, dat, timeRange = readSubjectFiles(dataPath, subjectName, dlms=True)

				#### Get timing data for x-axis
				## NOTE: for the line below, can use either df or dfs, in this context. 
				eventIdx = oscEventsInfo[band][eventType]['eventIdx']
				dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore = getStats(dfs, evidx=eventIdx,align='bywaveletpeak',verbose=True)
				print('minf: ' + str(minF))
				print('maxf: ' + str(maxF))

				## Resize w2 to match the load.py calculation for the osc event plotting (in def draw() in class eventviewer)
				w2 = int(w2*0.6)
				# print('w2: ' + str(w2))

				## Calculate beforeT
				idx0_before = max(0,left - w2)
				idx1_before = left 
				beforeT = (maxT-minT) * (idx1_before - idx0_before) / (right - left + 1)
				# print('beforeT: ' + str(beforeT))

				## Calculate afterT 
				idx0_after = int(right)
				idx1_after = min(idx0_after + w2,max(CSD.shape[0],CSD.shape[1]))
				afterT = (maxT-minT) * (idx1_after - idx0_after) / (right - left + 1)
				# print('afterT: ' + str(afterT))

				## Calculate tt for before: 
				sig_before = CSD[chan,idx0_before:idx1_before]
				tt_before = np.linspace(minT-beforeT,minT,len(sig_before)) + alignoffset

				## Calculate tt for during:
				sig_during = CSD[chan,left:right]
				tt_during = np.linspace(minT,maxT,len(sig_during)) + alignoffset

				## Calculate tt for after:
				sig_after = CSD[chan,idx0_after:idx1_after]
				tt_after = np.linspace(maxT,maxT+afterT,len(sig_after)) + alignoffset

				specrange = oscEventsInfo[band][eventType]['specrange']
				ylspec = oscEventsInfo[band][eventType]['ylspec']

				plotWavelets(dfs, df, dat, tt, sampr, dlms, subjectName, chan, band, eventIdx, specrange=specrange, ylspec=ylspec, saveFig=saveFig)


