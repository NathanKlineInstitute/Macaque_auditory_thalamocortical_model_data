
### IMPORTS ### 
from netpyne import sim
# from netpyne.analysis import csd
import csd 
import numpy as np
from matplotlib import pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker  ## for colorbar 
from morlet import MorletSpec, index2ms
import pandas as pd
from numbers import Number
import seaborn as sns 

def getWaveletInfo(freqBand, based, verbose=0): 
	## freqBand: str  --> e.g. 'delta', 'alpha', 'theta'
	## based: str --> path to directory with the .pkl data files 
	## verbose: bool --> if 0, default to only putting out timeRange and dataFile, if 1 --> include channel as well 

	waveletInfo = {
	'delta': {'dataFile': 'v34_batch57_3_4/A1_v34_batch65_v34_batch65_0_0_data.pkl', 'timeRange': [1480, 2520], 'channel': 14},
	'beta': {'dataFile': 'v34_batch57_3_2/A1_v34_batch67_v34_batch67_0_0_data.pkl',	'timeRange': [456, 572], 'channel': 14}, 
	'alpha': {'dataFile': 'v34_batch57_3_2/A1_v34_batch67_v34_batch67_0_0_data.pkl', 'timeRange': [3111, 3325], 'channel': 9}, 
	'theta': {'dataFile': 'v34_batch57_3_3/A1_v34_batch67_v34_batch67_1_1_data.pkl', 'timeRange': [2785, 3350], 'channel': 8}}

	timeRange = waveletInfo[freqBand]['timeRange']
	dataFileNoPath = waveletInfo[freqBand]['dataFile']
	dataFile = based + dataFileNoPath
	channel = waveletInfo[freqBand]['channel']

	if verbose:
		return timeRange, dataFile, channel
	else:
		return timeRange, dataFile

## CSD panels 
def getCSDdata(dataFile=None, outputType=['timeSeries', 'spectrogram'], oscEventInfo=None, dt=None, sampr=None, pop=None, spacing_um=100, minFreq=1, maxFreq=110, stepFreq=0.25):
	#### Outputs a dict with CSD and other relevant data for plotting! 
	## dataFile: str     					--> .pkl file with recorded simulation 
	## outputType: list of strings 			--> options are 'timeSeries' +/- 'spectrogram' --> OR could be empty, if want csdData from all electrodes!! 
	## oscEventInfo: dict 					---> dict w/ chan, left, right, minT, maxT, alignoffset
	## dt: time step of the simulation 		--> (usually --> sim.cfg.recordStep)
	## sampr: sampling rate (Hz) 			--> (usually --> 1/(dt/1000))
	## pop: str or list 					--> e.g. 'ITS4' or ['ITS4']
	## spacing_um: int 						--> 100 by DEFAULT (spacing between electrodes in MICRONS)
	## minFreq: float / int 				--> DEFAULT: 1 Hz  
	## maxFreq: float / int 				--> DEFAULT: 110 Hz 
	## stepFreq: float / int 				--> DEFAULT: 0.25 Hz 


	## load .pkl simulation file 
	if dataFile:
		sim.load(dataFile, instantiate=False)
	else:
		print('No dataFile; will use data from dataFile already loaded elsewhere!')

	## Determine timestep, sampling rate, and electrode spacing 
	dt = sim.cfg.recordStep#/1000.0	# This is in Seconds 
	sampr = 1.0 / (dt/1000.0) 				# Hz 
	spacing_um = spacing_um			# 100 um by default # 


	## Get LFP data   
	if pop is None:
		lfpData = sim.allSimData['LFP']
	else:
		if type(pop) is list:
			pop = pop[0]
		lfpData = sim.allSimData['LFPPops'][pop]

	## Set up dictionary for output data 
	outputData = {}


	#### ALL CSD DATA -- ALL CHANS, ALL TIMEPOINTS!!! 
	csdData = csd.getCSD(LFP_input_data=lfpData, dt=dt, sampr=sampr, spacing_um=spacing_um, vaknin=True)
	tt = np.linspace(0, sim.cfg.duration, len(csdData[1])) 
	# Input full CSD data (and time array) into outputData dict 
	outputData.update({'csdData': csdData})
	outputData.update({'tt': tt}) ## Keep in mind this is in milliseconds! 


	## Extract oscillation event info 
	if oscEventInfo is not None:
		# Extract chan, left, right, minT, maxT, alignoffset, w2 
		chan = oscEventInfo['chan']
		left = oscEventInfo['left']
		right = oscEventInfo['right']
		minT = oscEventInfo['minT']
		maxT = oscEventInfo['maxT']
		alignoffset = oscEventInfo['alignoffset']
		w2 = oscEventInfo['w2']
		## Calculate idx0 and idx1 for before, and beforeT
		idx0_before = max(0,left - w2)
		idx1_before = left 
		beforeT = (maxT-minT) * (idx1_before - idx0_before) / (right - left + 1)
		outputData.update({'beforeT': beforeT})
		# Calculate idx0 and idx1 for after, and afterT
		idx0_after = int(right)
		idx1_after = min(idx0_after + w2,max(csdData.shape[0],csdData.shape[1]))
		afterT = (maxT-minT) * (idx1_after - idx0_after) / (right - left + 1)
		outputData.update({'afterT': afterT})
	else:
		chan=None
		# print('No oscillation event data detected!')


		# timeSeries --------------------------------------------
	if 'timeSeries' in outputType:
		# print('Returning timeSeries data')

		##############################
		#### BEFORE THE OSC EVENT ####
		# (1) Calculate CSD before osc event
		csdBefore = csdData[chan,idx0_before:idx1_before]
		# (2) Calculate timepoint data 
		tt_before = np.linspace(minT-beforeT,minT,len(csdBefore)) + alignoffset
		# (3) Input time and CSD data for before osc event into outputData dict 
		outputData.update({'tt_before': tt_before})
		outputData.update({'csdBefore': csdBefore})

		##############################
		#### DURING THE OSC EVENT ####
		# (1) Calculate CSD during osc event
		csdDuring = csdData[chan,left:right]
		# (2) Calculate timepoint data 
		tt_during = np.linspace(minT,maxT,len(csdDuring)) + alignoffset
		outputData.update({'tt_during': tt_during})
		outputData.update({'csdDuring': csdDuring})

		#############################
		#### AFTER THE OSC EVENT #### 
		# (1) Calculate CSD after osc event
		csdAfter = csdData[chan,idx0_after:idx1_after]
		# (2) Calculate timepoint data 
		tt_after = np.linspace(maxT,maxT+afterT,len(csdAfter)) + alignoffset
		# (3) Input time and CSD data for after osc event into outputData dict 
		outputData.update({'tt_after': tt_after})
		outputData.update({'csdAfter': csdAfter})


		#################################################
		#### BEFORE, DURING, AND AFTER THE OSC EVENT #### 
		# (1) Calculate CSD after osc event
		csdOscChan_plusTimeBuffer = csdData[chan,idx0_before:idx1_after]
		# (2) Calculate timepoint data 
		tt_plusTimeBuffer = np.linspace(minT-beforeT,maxT+afterT,len(csdOscChan_plusTimeBuffer)) + alignoffset
		# (3) Input time and CSD data for after osc event into outputData dict 
		outputData.update({'tt_plusTimeBuffer': tt_plusTimeBuffer})
		outputData.update({'csdOscChan_plusTimeBuffer': csdOscChan_plusTimeBuffer})


		# Update outputData with the channel info as well
		outputData.update({'chan': chan})

		# Get xlim for plotting
		xl = (minT-beforeT + alignoffset, maxT+afterT + alignoffset)
		outputData.update({'xl': xl})


	# spectrogram -------------------------------------------
	if 'spectrogram' in outputType:
		# print('Returning spectrogram data')

		spec = []
		freqList = None

		## Spectrogram Data Calculations ## 
		fs = sampr 

		##############################
		#### DURING THE OSC EVENT #### 
		csdDuring = csdData[chan,left:right]
		specDuring = []
		specDuring.append(MorletSpec(csdDuring, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq, getphase=True,width=7.0))
		## vmin, vmax 
		vminDuring = np.array([s.TFR for s in specDuring]).min()
		vmaxDuring = np.array([s.TFR for s in specDuring]).max()
		vcDuring = [vminDuring, vmaxDuring]
		## T 
		T_during = [minT + alignoffset, maxT + alignoffset]
		## F, S 
		F_during = specDuring[0].f
		S_during = specDuring[0].TFR
		## outputData update 
		outputData.update({'T_during': T_during, 'F_during': F_during, 'S_during': S_during, 'vcDuring': vcDuring})
		outputData.update({'specDuring': specDuring})


		###############################################
		#### DURING THE OSC EVENT + BEFORE & AFTER ####  
		csdFull = csdData[chan, idx0_before:idx1_after] 
		specFull = []
		specFull.append(MorletSpec(csdFull, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq, getphase=True,width=7.0))
		## vmin, vmax 
		vminFull = np.array([s.TFR for s in specFull]).min()
		vmaxFull = np.array([s.TFR for s in specFull]).max()
		vcFull = [vminFull, vmaxFull]
		## T 
		T_full = [(minT - beforeT) + alignoffset, (maxT + afterT) + alignoffset]
		## F, S
		F_full = specFull[0].f
		S_full = specFull[0].TFR
		## outputData update 
		outputData.update({'T_full': T_full, 'F_full': F_full, 'S_full': S_full, 'vcFull': vcFull})
		outputData.update({'specFull': specFull})


		## Get frequency list 
		f = freqList if freqList is not None else np.arange(minFreq, maxFreq+stepFreq, stepFreq)
		outputData.update({'freqs': f[f<=maxFreq]})

		# Update outputData with the channel info as well
		outputData.update({'chan': chan})

	return outputData
def plotCombinedCSD(pop, electrode, colorDict, timeSeriesDict=None, spectDict=None, alignoffset=None, vmaxContrast=None, colorMap='jet', figSize=(10,7), minFreq=None, maxFreq=None, plotTypes=['timeSeries', 'spectrogram'], hasBefore=1, hasAfter=1, savePath=None, saveFig=True):
	### pop: list or str 				--> relevant population to plot data for 
	### electrode: int 					--> electrode at which to plot the CSD data 
	### colorDict: dict 				--> corresponds pop to color 
	### timeSeriesDict: dict 			--> output of getCSDdata (with outputType=['timeSeries'])
	### spectDict: dict 				--> output of getCSDdata (with outputType=['spectrogram'])
	### alignoffset: float 				!! --> Can be gotten from plotWavelets.py in a1dat [TO DO: ELABORATE ON THIS / LINK CALCULATIONS]
	### vmaxContrast: float or int 		--> Denominator; This will help with color contrast if desired!, e.g. 1.5 or 3 (DEFAULT: None)
	### colorMap: str 					--> DEFAULT: jet; cmap for ax.imshow lines --> Options are currently 'jet' or 'viridis'
	### figSize: tuple 					--> DEFAULT: (10,7)
	### minFreq: int 					--> DEFAULT: None
	### maxFreq: int 					--> DEFAULT: None
	### plotTypes: list 				--> DEFAULT: ['timeSeries', 'spectrogram']
	### hasBefore: bool 				--> DEFAULT: 1 (before osc event data will be plotted)
	### hasAfter: bool  				--> DEFAULT: 1 (after osc event data will be plotted)
	### savePath: str    				--> Path to directory where fig should be saved; DEFAULT: '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/popContribFigs/'
	### saveFig: bool 					--> DEFAULT: True

	# Get relevant pop
	if type(pop) is str:
		popToPlot = pop
	elif type(pop) is list:
		popToPlot = pop[0]
	elif pop is None:
		popToPlot = None


	# Create figure 
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figSize)

	# Set font sizes
	labelFontSize = 20#12  ## NOTE: spike has this as 12, lfp plotting has this as 15 
	titleFontSize = 25#20

	# Set electrode variable, for plot title(s)
	if type(electrode) is list:
		electrode = electrode[0]


	#### SPECTROGRAM CALCULATIONS ####-------------------------------------------------------
	if 'spectrogram' in plotTypes:

		if hasBefore and hasAfter:
			## S, F, T
			S = spectDict['S_full']
			F = spectDict['F_full']
			T = spectDict['T_full']
			## vmin, vmax 	-- adjust for color contrast purposes, if desired 
			vc = spectDict['vcFull']


		else:
			## S, F, T
			S = spectDict['S_during']
			F = spectDict['F_during']
			T = spectDict['T_during']
			## vmin, vmax 	-- adjust for color contrast purposes, if desired 
			vc = spectDict['vcDuring']


		orig_vmin = vc[0]
		orig_vmax = vc[1]

		if vmaxContrast is None:
			vmin = orig_vmin
			vmax = orig_vmax
		else:
			vmin = orig_vmin
			vmax = orig_vmax / vmaxContrast
			vc = [vmin, vmax]


		## Set up minFreq and maxFreq for spectrogram
		if minFreq is None:
			minFreq = np.amin(F) 	# 1
		if maxFreq is None:
			maxFreq = np.amax(F)	# 100

		## Set up imshowSignal
		imshowSignal = S 


	#### TIME SERIES CALCULATIONS ####-------------------------------------------------------
	if 'timeSeries' in plotTypes:
		# time (x-axis)
		tt_during = timeSeriesDict['tt_during']

		# CSD (y-axis)
		csdDuring = timeSeriesDict['csdDuring']

		if hasBefore:
			# time (x-axis)
			tt_before = timeSeriesDict['tt_before']
			# CSD (y-axis)
			csdBefore = timeSeriesDict['csdBefore']

		if hasAfter:
			# time (x-axis)
			tt_after = timeSeriesDict['tt_after']
			# CSD (y-axis)
			csdAfter = timeSeriesDict['csdAfter']

	#### PLOTTING ####-------------------------------------------------------
	if 'spectrogram' in plotTypes and 'timeSeries' in plotTypes:

		### PLOT SPECTROGRAM ### 
		# spectrogram title
		if popToPlot is None:
			spectTitle = 'CSD Spectrogram, channel ' + str(spectDict['chan'])
		else:
			spectTitle = 'CSD Spectrogram for ' + popToPlot + ', channel ' + str(spectDict['chan'])
		# plot and format 
		ax1 = plt.subplot(2, 1, 1)
		img = ax1.imshow(imshowSignal, extent=(np.amin(T), np.amax(T), minFreq, maxFreq), origin='lower', interpolation='None', aspect='auto', 
			vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap(colorMap))  # minFreq --> np.amin(F)  	# maxFreq --> np.amax(F)	# imshowSignal --> S
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad=0.2)
		cax1.tick_params(labelsize=15) #12)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar scientific notation
		fmt.set_powerlimits((0,0))
		cbar=plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		cbar.ax.set_ylabel('Power', size=15)
		ax1.set_title(spectTitle, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=T[0], right=T[1]) 			# ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq)				# Uncomment this if using the commented-out ax1.imshow (with S, and with np.amin(F) etc.)
		## Set the spectrogram tick params (font size) ## 
		ax1.tick_params(axis='both', labelsize=15) #12)


		### PLOT TIMESERIES ###
		# timeSeries title
		if popToPlot is None:
			timeSeriesTitle = 'CSD Signal, channel ' + str(timeSeriesDict['chan'])
		else:
			timeSeriesTitle = 'CSD Signal for ' + popToPlot + ', channel ' + str(timeSeriesDict['chan'])

		# y-axis label
		timeSeriesYAxis = 'CSD Amplitude (' + r'$\frac{mV}{mm^2}$' + ')'

		# Format timeSeries plot 
		lw = 1.0
		ax2 = plt.subplot(2, 1, 2)
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes('right', size='3%', pad=0.2)
		cax2.axis('off')
		# Plot CSD before Osc Event
		if hasBefore:  ### NOTE: combined with hasAfter I think
			ax2.plot(tt_before, csdBefore, color='k', linewidth=1.0)
		# Plot CSD during Osc Event
		if popToPlot is None:
			ax2.plot(tt_during, csdDuring, color='r', linewidth=lw) 		# ax2.plot(t, csdTimeSeries, color=colorDict[popToPlot], linewidth=lw) # 	# ax2.plot(t[0:len(lfpPlot)], lfpPlot, color=colorDict[popToPlot], linewidth=lw)
		else:
			ax2.plot(tt_during, csdDuring, color=colorDict[popToPlot], linewidth=lw) 		# ax2.plot(t, csdTimeSeries, color=colorDict[popToPlot], linewidth=lw) # 	# ax2.plot(t[0:len(lfpPlot)], lfpPlot, color=colorDict[popToPlot], linewidth=lw)
		# Plot CSD after Osc Event 
		if hasAfter:
			ax2.plot(tt_after, csdAfter, color='k', linewidth=1.0)
		# Set title and x-axis label
		ax2.set_title(timeSeriesTitle, fontsize=titleFontSize)
		ax2.set_xlabel('Time (ms)', fontsize=labelFontSize)
		# Set x-axis limits 
		if hasBefore and hasAfter:
			xl = timeSeriesDict['xl']
			ax2.set_xlim((xl))
		else:
			ax2.set_xlim(left=tt_during[0], right=tt_during[-1]) 		# ax2.set_xlim(left=t[0], right=t[-1]) 			# ax2.set_xlim(left=timeRange[0], right=timeRange[1])
		# Set y-axis label 
		ax2.set_ylabel(timeSeriesYAxis, fontsize=labelFontSize)
		# Set the timeSeries tick params (font size)
		ax2.tick_params(axis='both', labelsize=15) #12)

		# For potential saving 
		if popToPlot is None:
			figFilename = 'CombinedCSD_chan_' + str(timeSeriesDict['chan']) + '.png'
		else:
			figFilename = popToPlot + '_combinedCSD_chan_' + str(timeSeriesDict['chan']) + '.png' # figFilename = popToPlot + '_combinedCSD_chan_' + str(timeSeriesDict['chan']) + '.png'


	elif 'spectrogram' in plotTypes and 'timeSeries' not in plotTypes:
		### PLOT SPECTROGRAM ### 
		# spectrogram title
		if popToPlot is None:
			spectTitle = 'CSD Spectrogram, channel ' + str(spectDict['chan'])
		else:
			spectTitle = 'CSD Spectrogram for ' + popToPlot + ', channel ' + str(spectDict['chan'])  	# spectTitle = 'CSD Spectrogram for ' + popToPlot + ', channel ' + str(spectDict['chan'])
		# plot and format 
		ax1 = plt.subplot(1, 1, 1)
		# img = ax1.imshow(S, extent=(np.amin(T), np.amax(T), np.amin(F), np.amax(F)), origin='lower', interpolation='None', aspect='auto', 
		# 	vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap(colorMap)) ## ANSWER: NO --> NOTE: instead of np.amax(F) should I be doing maxFreq? (similarly for minFreq // np.amin(F))
		img = ax1.imshow(imshowSignal, extent=(np.amin(T), np.amax(T), minFreq, maxFreq), origin='lower', interpolation='None', aspect='auto', 
			vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap(colorMap)) 
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad=0.2)
		cax1.tick_params(labelsize=15) #12)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar scientific notation
		fmt.set_powerlimits((0,0))
		cbar=plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		cbar.ax.set_ylabel('Power', size=15)
		ax1.set_title(spectTitle, fontsize=titleFontSize)
		ax1.set_xlabel('Time (ms)', fontsize=labelFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=T[0], right=T[1]) 			# ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq) 				# Uncomment this if using the commented-out ax1.imshow (with S, and with np.amin(F) etc.)
		## Set the spectrogram tick params (font size)
		ax1.tick_params(axis='both', labelsize=15) #12)

		# For potential saving 
		if popToPlot is None:
			figFilename = 'CSD_spectrogram_chan_' + str(timeSeriesDict['chan']) + '.png'
		else:
			figFilename = popToPlot + '_CSD_spectrogram_chan_' + str(timeSeriesDict['chan']) + '.png'	# figFilename = popToPlot + '_CSD_spectrogram_chan_' + str(spectDict['chan']) + '.png'


	elif 'spectrogram' not in plotTypes and 'timeSeries' in plotTypes:

		### PLOT TIME SERIES ### 
		# timeSeries title
		if popToPlot is None:
			timeSeriesTitle = 'CSD Signal, channel ' + str(timeSeriesDict['chan'])
		else:
			timeSeriesTitle = 'CSD Signal for ' + popToPlot + ', channel ' + str(timeSeriesDict['chan'])	# timeSeriesTitle = 'CSD Signal for ' + popToPlot + ', channel ' + str(timeSeriesDict['chan'])

		# y-axis label
		timeSeriesYAxis = 'CSD Amplitude (' + r'$\frac{mV}{mm^2}$' + ')'

		# Format timeSeries plot 
		lw = 1.0
		ax2 = plt.subplot(1, 1, 1)
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes('right', size='3%', pad=0.2)
		cax2.axis('off')
		# Plot CSD before Osc Event  
		if hasBefore: 
			ax2.plot(tt_before, csdBefore, color='r', linewidth=1.0)
		# Plot CSD during Osc Event 			# 		ax2.plot(t, csdTimeSeries, color=colorDict[popToPlot], linewidth=lw) # 	# ax2.plot(t[0:len(lfpPlot)], lfpPlot, color=colorDict[popToPlot], linewidth=lw)
		if popToPlot is None:
			ax2.plot(tt_during, csdDuring, color='k', linewidth=lw) # 	# ax2.plot(t[0:len(lfpPlot)], lfpPlot, color=colorDict[popToPlot], linewidth=lw)
		else:
			ax2.plot(tt_during, csdDuring, color=colorDict[popToPlot], linewidth=lw) # 	# ax2.plot(t[0:len(lfpPlot)], lfpPlot, color=colorDict[popToPlot], linewidth=lw)
		# Plot CSD after Osc Event 
		if hasAfter:
			ax2.plot(tt_after, csdAfter, color='k', linewidth=1.0)
		# Set title and x-axis label
		ax2.set_title(timeSeriesTitle, fontsize=titleFontSize)
		ax2.set_xlabel('Time (ms)', fontsize=labelFontSize)
		# Set x-axis limits 
		if hasBefore and hasAfter:
			xl = timeSeriesDict['xl']
			ax2.set_xlim((xl))
		elif hasBefore and not hasAfter:
			ax2.set_xlim(left=tt_before[0], right=tt_during[-1])
		elif not hasBefore and hasAfter:
			ax2.set_xlim(left=tt_during[0], right=tt_after[-1])
		else:
			ax2.set_xlim(left=tt_during[0], right=tt_during[-1])	# ax2.set_xlim(left=t[0], right=t[-1]) 	# ax2.set_xlim(left=timeRange[0], right=timeRange[1])
		# Set y-axis label 
		ax2.set_ylabel(timeSeriesYAxis, fontsize=labelFontSize)
		## Set the timeSeries tick params (font size)
		ax2.tick_params(axis='both', labelsize=15) #12)

		# For potential saving 
		if popToPlot is None:
			figFilename = 'CSD_timeSeries_chan_' + str(timeSeriesDict['chan']) + '.png'
		else:
			figFilename = popToPlot + '_CSD_timeSeries_chan_' + str(timeSeriesDict['chan']) + '.png'
 
	plt.tight_layout()
def plotCSDPanels(based):
	### based: str 		--> path to directory with sim files 
	colorDictCustom = {}
	colorDictCustom['ITP4'] = 'magenta'
	colorDictCustom['ITS4'] = 'dodgerblue'
	colorDictCustom['IT5A'] = 'lime'

	timeRange, dataFile, channel = getWaveletInfo('theta', based=based, verbose=1)

	### OSC EVENT INFO DICTS !!
	thetaOscEventInfo = {'chan': 8, 'minT': 2785.22321038684, 
					'maxT': 3347.9278996316607, 'alignoffset':-3086.95, 'left': 55704, 'right':66958,
					'w2': 3376} 

	### plot CSD spect + timeSeries for theta osc event, chan 8
	electrode=channel
	includePops=[None, 'ITS4', 'ITP4', 'IT5A']
	minFreq = 1
	maxFreq = 12
	stepFreq = 0.25 
	for pop in includePops:
		timeSeriesDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=thetaOscEventInfo, pop=pop, minFreq=minFreq, maxFreq=maxFreq, stepFreq=stepFreq)
		spectDict = getCSDdata(dataFile=dataFile, outputType=['spectrogram'], oscEventInfo=thetaOscEventInfo, pop=pop, minFreq=minFreq, maxFreq=maxFreq, stepFreq=stepFreq)

		plotCombinedCSD(timeSeriesDict=timeSeriesDict, spectDict=spectDict, colorDict=colorDictCustom, pop=pop, electrode=electrode, 
			minFreq=1, maxFreq=maxFreq, vmaxContrast=None, colorMap='jet', figSize=(10,7), plotTypes=['timeSeries', 'spectrogram'], 
			hasBefore=1, hasAfter=1, saveFig=False) 

	plt.show()

## CSD heatmap 
def getCSDDataFrames(dataFile, absolute=None, oscEventInfo=None, verbose=0):
	## This function will return data frames of peak and average CSD amplitudes, for picking cell pops
	### dataFile: str 				--> .pkl file to load, with data from the whole recording
	### absolute: bool 				--> determines whether to output absolute CSD values in dataframe 
	### oscEvenfInfo: dict 			--> dict w/ chan, left, right, minT, maxT, alignoffset
	### verbose: bool 

	if absolute is None:
		absolute=1

	sim.load(dataFile, instantiate=False)
	# Get args for getCSDdata
	dt = sim.cfg.recordStep 	# sim.cfg.recordStep/1000.0
	sampr = 1.0/(dt/1000.0) 	# divide by 1000.0 to turn denominator from units of ms to s  # sampr = 1.0 / dt 
	spacing_um = 100 

	# Get all cell pops (cortical)
	thalPops = ['TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']
	allPops = list(sim.net.allPops.keys())
	pops = [pop for pop in allPops if pop not in thalPops] 			## exclude thal pops 

	## Get all electrodes
	evalElecs = []
	evalElecs.extend(list(range(int(sim.net.recXElectrode.nsites))))
	### add 'avg' to electrode list 
	evalElecs.append('avg')


	## get CSD data
	csdPopData = {}

	for pop in pops:
		csdPopData[pop] = {}

		popCSDdataFULL_origShape_dict = getCSDdata(dataFile=dataFile, dt=dt, sampr=sampr, pop=pop, spacing_um=spacing_um, outputType=[])
		popCSDdataFULL_origShape = popCSDdataFULL_origShape_dict['csdData']
		popCSDdataFULL = np.transpose(popCSDdataFULL_origShape)	### TRANSPOSE THIS so (20,230000) --> (230000, 20)

		if oscEventInfo is not None:
			minT = oscEventInfo['minT']
			maxT = oscEventInfo['maxT']
			timeRange = [minT, maxT]
			popCSDdata = popCSDdataFULL[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]
		else:
			popCSDdata = popCSDdataFULL.copy()


		for i, elec in enumerate(evalElecs): ### HOW IS THIS GOING TO WORK WITH CSD VS LFP?? HM -- VAKNIN GOOD ENOUGH TO SOLVE THIS PROBLEM? 
			if elec == 'avg':
				csdPopData[pop]['avg'] = {}

				avgPopData = np.mean(popCSDdata, axis=1)	## csd data (from 1 pop) at each timepoint, averaged over all electrodes

				avgAvgCSD = np.average(avgPopData)
				csdPopData[pop]['avg']['avg'] = avgAvgCSD	## time-average of CSD data (from 1 pop) that has been averaged in space (over all electrodes)

				peakAvgCSD = np.amax(avgPopData)
				csdPopData[pop]['avg']['peak'] = peakAvgCSD	## highest datapoint of all the CSD data (from 1 pop) that has been averaged in space (over all electrodes)

			elif isinstance(elec, Number):
				elecKey = 'elec' + str(elec)
				csdPopData[pop][elecKey] = {}
				csdPopData[pop][elecKey]['avg'] = np.average(popCSDdata[:, elec])	## CSD data from 1 pop, averaged in time, over 1 electrode 
				csdPopData[pop][elecKey]['peak'] = np.amax(popCSDdata[:, elec])		## Maximum CSD value from 1 pop, over time, recorded at 1 electrode 



	#### PEAK CSD AMPLITUDES, DATA FRAME ####
	peakValues = {}
	peakValues['pops'] = []
	peakValues['peakCSD'] = [[] for i in range(len(pops))]  # should be 36? 
	p=0
	for pop in pops:
		peakValues['pops'].append(pop)
		for i, elec in enumerate(evalElecs): 
			if isinstance(elec, Number):
				elecKey = 'elec' + str(elec)
			elif elec == 'avg': 
				elecKey = 'avg'
			# peakValues['peakCSD'][p].append(csdPopData[pop][elecKey]['peak'])
			if absolute:
				peakValues['peakCSD'][p].append(np.absolute(csdPopData[pop][elecKey]['peak']))
			else:
				peakValues['peakCSD'][p].append(csdPopData[pop][elecKey]['peak'])
		p+=1
	dfPeak = pd.DataFrame(peakValues['peakCSD'], index=pops)


	#### AVERAGE LFP AMPLITUDES, DATA FRAME ####
	avgValues = {}
	avgValues['pops'] = []
	avgValues['avgCSD'] = [[] for i in range(len(pops))]
	q=0
	for pop in pops:
		avgValues['pops'].append(pop)
		for i, elec in enumerate(evalElecs):
			if isinstance(elec, Number):
				elecKey = 'elec' + str(elec)
			elif elec == 'avg':
				elecKey = 'avg'
			# avgValues['avgCSD'][q].append(csdPopData[pop][elecKey]['avg'])
			if absolute:
				avgValues['avgCSD'][q].append(np.absolute(csdPopData[pop][elecKey]['avg']))
			else:
				avgValues['avgCSD'][q].append(csdPopData[pop][elecKey]['avg'])
		q+=1
	dfAvg = pd.DataFrame(avgValues['avgCSD'], index=pops)


	# return csdPopData
	if verbose:
		return dfPeak, dfAvg, peakValues, avgValues, csdPopData 
	else:
		return dfPeak, dfAvg
def plotDataFrames(dataFrame, absolute=None, electrodes=None, pops=None, title=None, cbarLabel=None, figSize=None, savePath=None, saveFig=True):
	#### --> This function will plot a heatmap of the peak or average LFP amplitudes across electrodes & cell populations
	### dataFrame: pandas dataFrame  --> These can be obtained from getDataFrames function above)
	### absolute=1 			--> determines if absolute CSD values will be plotted or not 
	### electrodes: list 	--> DEFAULT: use all electrodes + 'avg'
	### pops: list 			--> DEFAULT: all cortical pops 
	### title: str  		--> Optional; title of the entire figure
	### cbarLabel: str 		--> DEFAULT: 'LFP amplitudes (mV)'  -->  (label on the color scale bar)
	### figSize: tuple 		--> DEFAULT: (12,6)
	### savePath: str 		--> path to directory where figures should be saved 
	### saveFig: bool 		-->  DEFAULT: True 

	### NOTE: I SHOULD HAVE THIS FUNCTION CALL CSD DATA FRAMES FX!!! 

	if absolute is None:
		absolute=1

	## Set label for color scalebar 
	if cbarLabel is None:
		cbarLabel = 'LFP amplitude (mV)'
	elif cbarLabel is 'CSD':
		if absolute:
			cbarLabel = 'Absolute value of CSD amplitude (' + r'$\frac{mV}{mm^2}$' + ')'
		else:
			cbarLabel = 'CSD amplitude (' + r'$\frac{mV}{mm^2}$' + ')'

	## Create lists of electrode (columns and labels)
	if electrodes is None:
		electrodeColumns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
		electrodeLabels = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 'avg']
	else:
		if 'avg' in electrodes:
			electrodeColumns = electrodes.copy()
			avgIndex = electrodeColumns.index('avg')
			electrodeColumns[avgIndex] = 20  				## <-- Ensures that the correct index is used to access the 'avg' data
		else:
			electrodeColumns = electrodes.copy() 
		electrodeLabels = electrodes.copy() 


	## dataFrame subset, according to the electrodes specified in the 'electrodes' argument! 
	dataFrame = dataFrame[electrodeColumns]

	## Create list of cell populations 
	if pops is None:  ### POTENTIAL ISSUE WITH ORDERING OR NO???
		pops = ['NGF1', 'IT2', 'SOM2', 'PV2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 
		'ITP4', 'ITS4', 'SOM4', 'PV4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'SOM5A', 'PV5A', 'VIP5A', 'NGF5A', 
		'IT5B', 'CT5B', 'PT5B', 'SOM5B', 'PV5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'SOM6', 'PV6', 'VIP6', 'NGF6']


	## dataFrame subset according to cell populations specified in argument! 
	dataFramePops = dataFrame[dataFrame.index.isin(pops)]

	## TRANSPOSE DATA FRAME 
	pivotedDataFrame = dataFramePops.T

	## Set Font Sizes for the Heatmap Plot 
	titleFontSize = 30	#25	#20
	labelFontSize = 25	#20	#15
	tickFontSize = 20	#15 #10


	## Create lists of x and y axis labels 
	x_axis_labels = pops.copy() 
	y_axis_labels = electrodeLabels.copy() 


	## Set size of figure 
	if figSize is None:
		figSize = (12,6) 
	plt.figure(figsize = figSize)


	## Set title of figure 
	if title is not None:
		plt.title(title, fontsize=titleFontSize)


	## Create heatmap! 
	if absolute:
		ax = sns.heatmap(pivotedDataFrame, xticklabels=x_axis_labels, yticklabels=y_axis_labels, 
						linewidth=0.4, cmap='jet', cbar_kws={'label': cbarLabel})  # center=7, 
	else:
		ax = sns.heatmap(pivotedDataFrame, xticklabels=x_axis_labels, yticklabels=y_axis_labels, 
						linewidth=0.4, cmap='jet', cbar_kws={'label': cbarLabel}) 	# center=4, 

	## Change fontsize of heatmap label
	ax.figure.axes[-1].yaxis.label.set_size(18)
	ax.figure.axes[-1].tick_params(labelsize=15)

	## Set labels on x and y axes 
	plt.xlabel('Cell populations', fontsize=labelFontSize)
	plt.xticks(rotation=45, fontsize=tickFontSize)
	plt.ylabel('Channel', fontsize=labelFontSize)
	plt.yticks(rotation=0, fontsize=tickFontSize) 

	plt.tight_layout()

	plt.show()

	return ax
def plotHeatmap(based, absolute=1):
	thetaOscEventInfo = {'chan': 8, 'minT': 2785.22321038684, 
					'maxT': 3347.9278996316607, 'alignoffset':-3086.95, 'left': 55704, 'right':66958,
					'w2': 3376} 

	timeRange, dataFile, channel = getWaveletInfo('theta', based=based, verbose=1)

	ECortPops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B', 'PT5B', 'IT6', 'CT6']

	dfCSDPeak, dfCSDAvg = getCSDDataFrames(dataFile=dataFile, oscEventInfo=thetaOscEventInfo, absolute=absolute)
	avgCSDPlot = plotDataFrames(dfCSDAvg, absolute=absolute, electrodes=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], pops=ECortPops, title='Avg CSD Values', cbarLabel='CSD', figSize=(10,8), savePath=None, saveFig=False) # figSize=(10,7)
					# electrodes=None

## Spike panels 
def getRateSpectrogramData(include=['allCells', 'eachPop'], oscEventInfo=None, binSize=5, minFreq=1, maxFreq=100, stepFreq=1, NFFT=256, noverlap=128, smooth=0, transformMethod = 'morlet', norm=False):
	""" MODIFIED FROM NETPYNE 
	include : list
		<Short description of include>
		**Default:** ``['allCells', 'eachPop']``
		**Options:** ``<option>`` <description of option>

	oscEventInfo : dict 
		Dict with information about the oscillation event 
		--> chan, left, right, minT, maxT, alignoffset, w2

	binSize : int
		<Short description of binSize>
		**Default:** ``5``
		**Options:** ``<option>`` <description of option>

	minFreq : int
		<Short description of minFreq>
		**Default:** ``1``
		**Options:** ``<option>`` <description of option>

	maxFreq : int
		<Short description of maxFreq>
		**Default:** ``100``
		**Options:** ``<option>`` <description of option>

	stepFreq : int
		<Short description of stepFreq>
		**Default:** ``1``
		**Options:** ``<option>`` <description of option>

	NFFT : int
		<Short description of NFFT>
		**Default:** ``256``
		**Options:** ``<option>`` <description of option>

	noverlap : int
		<Short description of noverlap>
		**Default:** ``128``
		**Options:** ``<option>`` <description of option>

	smooth : int
		<Short description of smooth>
		**Default:** ``0``
		**Options:** ``<option>`` <description of option>

	transformMethod : str
		<Short description of transformMethod>
		**Default:** ``'morlet'``
		**Options:** ``<option>`` <description of option>

	norm : bool
		<Short description of norm>
		**Default:** ``False``
		**Options:** ``<option>`` <description of option> 
		"""
	print('Getting firing rate spectrogram data ...')

	## Replace 'eachPop' with list of pops
	if 'eachPop' in include:
		include.remove('eachPop')
		for pop in sim.net.allPops: include.append(pop)

	## Set up outputData dict
	outputData = {}

	## Extract oscillation event info 
	if oscEventInfo is not None:
		# Extract left, right, minT, maxT, alignoffset, w2   # RECALL: HAVE TO CORRECT FOR 6_11 !!! THIS WORKS FOR 0_6 AS IS!!! 
		left = oscEventInfo['left']
		right = oscEventInfo['right']
		minT = oscEventInfo['minT']
		maxT = oscEventInfo['maxT']
		alignoffset = oscEventInfo['alignoffset']
		w2 = oscEventInfo['w2']
		# Calculate idx0 and idx1 for before, and beforeT
		idx0_before = max(0,left - w2)
		idx1_before = left 
		beforeT = (maxT-minT) * (idx1_before - idx0_before) / (right - left + 1)
		# Calculate idx0 and idx1 for after, and afterT
		idx0_after = int(right)
		idx1_after = min(idx0_after + w2, 230000)	# max(csdData.shape[0],csdData.shape[1]))  # <-- believe this would be the number of time points (~230000?)
		afterT = (maxT-minT) * (idx1_after - idx0_after) / (right - left + 1)		#		
		# Update outputData with any values necessary for the plotting function
		outputData.update({'alignoffset': alignoffset})
	# else:
	# 	print('No oevent dict')

	## Calculate time ranges for data gathering 		-->  ### NOTE: the time ranges WITHOUT ALIGNOFFSET ADDED TO EITHER ELEMENT OF THE LISTS ensures that spike counts remain accurate!! 
	# during osc event
	T_during_withoutAlignOffset = [minT, maxT]  # + alignoffset to both elements 
	T_during = [minT+alignoffset, maxT+alignoffset]  
	# during osc event + before & after 
	T_full_withoutAlignOffset = [(minT-beforeT), (maxT+afterT)]	 # + alignoffset to both elements 
	T_full = [(minT-beforeT) + alignoffset, (maxT+afterT) + alignoffset]

	# update outputData dict 
	outputData.update({'T_during': T_during, 'T_full': T_full})


	# histData = []
	histDataDuring = []		# hist data during the oscillation event
	histDataFull = []		# hist data during osc event + buffer periods before and after 


	# allSignal, allFreqs = [], []
	allSignalDuring, allFreqsDuring = [], []
	allSignalFull, allFreqsFull = [], []


	# Plot separate line for each entry in include
	for iplot,subset in enumerate(include):
		from netpyne.analysis.utils import getCellsInclude
		cells, cellGids, netStimLabels = getCellsInclude([subset])
		numNetStims = 0

		# Select cells to include
		if len(cellGids) > 0:
			try:
				spkinds,spkts = list(zip(*[(spkgid,spkt) for spkgid,spkt in zip(sim.allSimData['spkid'],sim.allSimData['spkt']) if spkgid in cellGids]))
			except:
				spkinds,spkts = [],[]
		else:
			spkinds,spkts = [],[]


		# Add NetStim spikes
		spkts, spkinds = list(spkts), list(spkinds)
		numNetStims = 0
		if 'stims' in sim.allSimData:
			for netStimLabel in netStimLabels:
				netStimSpks = [spk for cell,stims in sim.allSimData['stims'].items() \
					for stimLabel,stimSpks in stims.items() for spk in stimSpks if stimLabel == netStimLabel]
				if len(netStimSpks) > 0:
					lastInd = max(spkinds) if len(spkinds)>0 else 0
					spktsNew = netStimSpks
					spkindsNew = [lastInd+1+i for i in range(len(netStimSpks))]
					spkts.extend(spktsNew)
					spkinds.extend(spkindsNew)
					numNetStims += 1

	# Generate list of spike times shifted by alignoffset
	spkts_withoutAlignOffset = spkts.copy()
	outputData.update({'spkts_withoutAlignOffset': spkts_withoutAlignOffset})	#, 'spkinds': spkinds})
	spkts_offset = [spkt+alignoffset for spkt in spkts_withoutAlignOffset]
	outputData.update({'spkts_offset': spkts_offset})


	# histo = np.histogram(spkts, bins = np.arange(timeRange[0], timeRange[1], binSize))
	histoDuring = np.histogram(spkts_offset, bins = np.arange(T_during[0], T_during[1], binSize))
	histoFull = np.histogram(spkts_offset, bins = np.arange(T_full[0], T_full[1], binSize))
	# histoT = histo[1][:-1]+binSize/2
	histoTDuring = histoDuring[1][:-1]+binSize/2
	histoTFull = histoFull[1][:-1]+binSize/2
	# histoCount = histo[0]
	histoCountDuring = histoDuring[0]
	histoCountFull = histoFull[0]
	# histoCount = histoCount * (1000.0 / binSize) / (len(cellGids)+numNetStims) # convert to rates
	histoCountDuring = histoCountDuring * (1000.0 / binSize) / (len(cellGids)+numNetStims) # convert to firing rate
	histoCountFull = histoCountFull * (1000.0 / binSize) / (len(cellGids)+numNetStims) # convert to firing rate
	# histData.append(histoCount)
	histDataDuring.append(histoCountDuring)
	histDataFull.append(histoCountFull)



	# Morlet wavelet transform method
	if transformMethod == 'morlet':
		from morlet import MorletSpec, index2ms

		Fs = 1000.0 / binSize

		# morletSpec = MorletSpec(histoCount, Fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq)
		morletSpecDuring = MorletSpec(histoCountDuring, Fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq)
		morletSpecFull = MorletSpec(histoCountFull, Fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq)
		# freqs = morletSpec.f
		freqsDuring = morletSpecDuring.f
		freqsFull = morletSpecFull.f
		# spec = morletSpec.TFR
		specDuring = morletSpecDuring.TFR
		specFull = morletSpecFull.TFR

		# allSignal.append(spec)
		allSignalDuring.append(specDuring)
		allSignalFull.append(specFull)
		# allFreqs.append(freqs)
		allFreqsDuring.append(freqsDuring)
		allFreqsFull.append(freqsFull)


		## vmin, vmax
		vminDuring = np.array(specDuring).min()
		vmaxDuring = np.array(specDuring).max()
		vcDuring = [vminDuring, vmaxDuring]
		outputData.update({'vcDuring': vcDuring})

		vminFull = np.array(specFull).min()
		vmaxFull = np.array(specFull).max()
		vcFull = [vminFull, vmaxFull]
		outputData.update({'vcFull': vcFull})

	# save figure data
	outputData.update({'allSignalDuring': allSignalDuring, 'allFreqsDuring': allFreqsDuring, 
					'allSignalFull': allSignalFull, 'allFreqsFull': allFreqsFull})


	return outputData
def getSpikeHistData(include=['eachPop', 'allCells'], oscEventInfo=None, binSize=5, graphType='bar', measure='rate', norm=False, smooth=None, filtFreq=None, filtOrder=3, axis=True, **kwargs):
	""" MODIFIED FROM NETPYNE 
	include : list
		Populations and cells to include in the plot.
		**Default:**
		``['eachPop', 'allCells']`` plots histogram for each population and overall average
		**Options:**
		``['all']`` plots all cells and stimulations,
		``['allNetStims']`` plots just stimulations,
		``['popName1']`` plots a single population,
		``['popName1', 'popName2']`` plots multiple populations,
		``[120]`` plots a single cell,
		``[120, 130]`` plots multiple cells,
		``[('popName1', 56)]`` plots a cell from a specific population,
		``[('popName1', [0, 1]), ('popName2', [4, 5, 6])]``, plots cells from multiple populations

	oscEventInfo : dict 
		Dict with information about the oscillation event 
		--> chan, left, right, minT, maxT, alignoffset, w2

	binSize : int
		Size of bin in ms to use for spike histogram.
		**Default:** ``5``
		**Options:** ``<option>`` <description of option>

	graphType : str
		Show histograms as line graphs or bar plots.
		**Default:** ``'bar'``
		**Options:** ``'line'``

	measure : str
		Whether to plot spike freguency (rate) or spike count.
		**Default:** ``'rate'``
		**Options:** ``'count'``

	norm : bool
		Whether to normalize the data or not.
		**Default:** ``False`` does not normalize the data
		**Options:** ``<option>`` <description of option>

	smooth : int
		Window width for smoothing.
		**Default:** ``None`` does not smooth the data
		**Options:** ``<option>`` <description of option>

	filtFreq : int or list
		Frequency for low-pass filter (int) or frequencies for bandpass filter in a list: [low, high]
		**Default:** ``None`` does not filter the data
		**Options:** ``<option>`` <description of option>

	filtOrder : int
		Order of the filter defined by `filtFreq`.
		**Default:** ``3``
		**Options:** ``<option>`` <description of option>

	axis : bool
		Whether to include a labeled axis on the figure.
		**Default:** ``True`` includes a labeled axis
		**Options:** ``False`` includes a scale bar

	kwargs : <type>
		<Short description of kwargs>
		**Default:** *required*
	"""

	# from .. import sim  ### <--- SHOULD ALREADY HAVE THIS

	print('Getting spike histogram data...')

	# Replace 'eachPop' with list of pops
	if 'eachPop' in include:
		include.remove('eachPop')
		for pop in sim.net.allPops: include.append(pop)


	# Set up outputData dict
	outputData = {}


	## Extract oscillation event info 
	if oscEventInfo is not None:
		# Extract left, right, minT, maxT, alignoffset, w2   # RECALL: HAVE TO CORRECT FOR 6_11 !!! THIS WORKS FOR 0_6 AS IS!!! 
		left = oscEventInfo['left']
		right = oscEventInfo['right']
		minT = oscEventInfo['minT']
		maxT = oscEventInfo['maxT']
		alignoffset = oscEventInfo['alignoffset']
		w2 = oscEventInfo['w2']
		# Calculate idx0 and idx1 for before, and beforeT
		idx0_before = max(0,left - w2)
		idx1_before = left 
		beforeT = (maxT-minT) * (idx1_before - idx0_before) / (right - left + 1)
		# Calculate idx0 and idx1 for after, and afterT
		idx0_after = int(right)
		idx1_after = min(idx0_after + w2, 230000)	# max(csdData.shape[0],csdData.shape[1]))  # <-- believe this would be the number of time points (~230000?)
		afterT = (maxT-minT) * (idx1_after - idx0_after) / (right - left + 1)		#		
		# Update outputData with any values necessary for the plotting function
		outputData.update({'alignoffset': alignoffset})
	# else:
	# 	print('No oscillation event data detected!')


	# Calculate time range to gather data for (during oscillation event only, or during time of oscillation event + buffer before and after)
	# before osc event
	T_before_withoutAlignOffset = [(minT-beforeT), minT] # + alignoffset to both elements 
	T_before = [(minT-beforeT) + alignoffset, minT + alignoffset]
	# during osc event
	T_during_withoutAlignOffset = [minT, maxT]  # + alignoffset to both elements 
	T_during = [minT+alignoffset, maxT+alignoffset]  
	# after osc event 
	T_after_withoutAlignOffset = [maxT, (maxT+afterT)] # + alignoffset to both elements
	T_after = [maxT+alignoffset, (maxT+afterT)+alignoffset]
	# update outputData dict 
	outputData.update({'T_before': T_before, 'T_during': T_during, 'T_after': T_after})

	# Histogram data 
	histoDataBefore = []
	histoDataDuring = []
	histoDataAfter = []

	# Plot separate line for each entry in include
	for iplot,subset in enumerate(include):
		from netpyne.analysis.utils import getCellsInclude
		if isinstance(subset, list):
			cells, cellGids, netStimLabels = getCellsInclude(subset)
		else:
			cells, cellGids, netStimLabels = getCellsInclude([subset])
		numNetStims = 0

		# Select cells to include
		if len(cellGids) > 0:
			try:
				spkinds,spkts = list(zip(*[(spkgid,spkt) for spkgid,spkt in zip(sim.allSimData['spkid'],sim.allSimData['spkt']) if spkgid in cellGids]))
			except:
				spkinds,spkts = [],[]
		else:
			spkinds,spkts = [],[]

		# Add NetStim spikes
		spkts, spkinds = list(spkts), list(spkinds)
		numNetStims = 0
		if 'stims' in sim.allSimData:
			for netStimLabel in netStimLabels:
				netStimSpks = [spk for cell,stims in sim.allSimData['stims'].items() \
				for stimLabel,stimSpks in stims.items() for spk in stimSpks if stimLabel == netStimLabel]
				if len(netStimSpks) > 0:
					lastInd = max(spkinds) if len(spkinds)>0 else 0
					spktsNew = netStimSpks
					spkindsNew = [lastInd+1+i for i in range(len(netStimSpks))]
					spkts.extend(spktsNew)
					spkinds.extend(spkindsNew)
					numNetStims += 1

	# Generate list of spike times shifted by alignoffset
	spkts_withoutAlignOffset = spkts.copy()
	outputData.update({'spkts_withoutAlignOffset': spkts_withoutAlignOffset})	#, 'spkinds': spkinds})
	spkts_offset = [spkt+alignoffset for spkt in spkts_withoutAlignOffset]
	outputData.update({'spkts_offset': spkts_offset})


	# histo = np.histogram(spkts, bins = np.arange(timeRange[0], timeRange[1], binSize))
	histoBefore = np.histogram(spkts_offset, bins = np.arange(T_before[0], T_before[1], binSize))	# first arg used to be spkts
	histoDuring = np.histogram(spkts_offset, bins = np.arange(T_during[0], T_during[1], binSize))	# first arg used to be spkts
	histoAfter = np.histogram(spkts_offset, bins = np.arange(T_after[0], T_after[1], binSize))		# first arg used to be spkts
	# histoT = histo[1][:-1]+binSize/2
	histoTBefore = histoBefore[1][:-1]+binSize/2
	histoTDuring = histoDuring[1][:-1]+binSize/2
	histoTAfter = histoAfter[1][:-1]+binSize/2
	# histoCount = histo[0]
	histoCountBefore = histoBefore[0]
	histoCountDuring = histoDuring[0]
	histoCountAfter = histoAfter[0]


	if measure == 'rate':
		# histoCount = histoCount * (1000.0 / binSize) / (len(cellGids)+numNetStims) # convert to firing rate
		histoCountBefore = histoCountBefore * (1000.0 / binSize) / (len(cellGids)+numNetStims) # convert to firing rate
		histoCountDuring = histoCountDuring * (1000.0 / binSize) / (len(cellGids)+numNetStims) # convert to firing rate
		histoCountAfter = histoCountAfter * (1000.0 / binSize) / (len(cellGids)+numNetStims) # convert to firing rate


	if filtFreq:
		from scipy import signal
		fs = 1000.0/binSize
		nyquist = fs/2.0
		if isinstance(filtFreq, list): # bandpass
			Wn = [filtFreq[0]/nyquist, filtFreq[1]/nyquist]
			b, a = signal.butter(filtOrder, Wn, btype='bandpass')
		elif isinstance(filtFreq, Number): # lowpass
			Wn = filtFreq/nyquist
			b, a = signal.butter(filtOrder, Wn)
		# histoCount = signal.filtfilt(b, a, histoCount)
		histoCountBefore = signal.filtfilt(b, a, histoCountBefore)
		histoCountDuring = signal.filtfilt(b, a, histoCountDuring)
		histoCountAfter = signal.filtfilt(b, a, histoCountAfter)

	if norm:
		# histoCount /= max(histoCount)
		histoCountBefore /= max(histoCountBefore)
		histoCountDuring /= max(histoCountDuring)
		histoCountAfter /= max(histoCountAfter)

	if smooth:
		# histoCount = _smooth1d(histoCount, smooth)[:len(histoT)]  ## get smooth1d from netpyne.analysis.utils if necessary
		histoCountBefore = _smooth1d(histoCountBefore, smooth)[:len(histoTBefore)]  ## get smooth1d from netpyne.analysis.utils if necessary
		histoCountDuring = _smooth1d(histoCountDuring, smooth)[:len(histoTDuring)]  ## get smooth1d from netpyne.analysis.utils if necessary
		histoCountAfter = _smooth1d(histoCountAfter, smooth)[:len(histoTAfter)]  ## get smooth1d from netpyne.analysis.utils if necessary

	# histoData.append(histoCount)  
	histoDataBefore.append(histoCountBefore)
	histoDataDuring.append(histoCountDuring)
	histoDataAfter.append(histoCountAfter)

	# save figure data
	outputData.update({'histoDataBefore': histoDataBefore, 'histoDataDuring': histoDataDuring, 'histoDataAfter': histoDataAfter, 
				'histoTBefore': histoTBefore, 'histoTDuring': histoTDuring, 'histoTAfter': histoTAfter, 
				'include': include, 'binSize': binSize})

	return outputData 
def getSpikeData(dataFile, pop, graphType, oscEventInfo): 
	### dataFile: path to .pkl data file to load 
	### pop: list or str --> which pop to include 
	### graphType: str --> either 'hist' or 'spect'
	### oscEventInfo: dict 

	# Load data file
	sim.load(dataFile, instantiate=False)

	# Pops
	if type(pop) is str:
		popList = [pop]
	elif type(pop) is list:
		popList = pop

	# Set up which kind of data -- i.e. spectrogram or histogram 
	if graphType is 'spect':
		spikeDict = getRateSpectrogramData(include=popList, oscEventInfo=oscEventInfo)
	elif graphType is 'hist':
		spikeDict = getSpikeHistData(include=popList, oscEventInfo=oscEventInfo, binSize=5, graphType='bar', measure='rate') ## sim.analysis.getSpikeHistData

	return spikeDict 
def plotCombinedSpike(pop, colorDict, plotTypes=['spectrogram', 'histogram'], hasBefore=1, hasAfter=1, spectDict=None, histDict=None, figSize=(10,7), colorMap='jet', minFreq=None, maxFreq=None, vmaxContrast=None, savePath=None, saveFig=True):
	### pop: str or list of length 1 	--> population to include 
	### colorDict: dict 				--> dict that corresponds pops to colors 
	### plotTypes: list 				--> ['spectrogram', 'histogram']
	### hasBefore: bool 				--> plot buffer before the osc event
	### hasAfter: bool 					--> plot buffer after the osc event
	### spectDict: dict 				--> can be gotten with getSpikeData(graphType='spect')
	### histDict: dict  				--> can be gotten with getSpikeData(graphType='hist')
	### figSize: tuple 					--> DEFAULT: (10,7)
	### colorMap: str 					--> DEFAULT: 'jet' 	--> cmap for ax.imshow lines --> Options are currently 'jet' or 'viridis' 
	### minFreq: int 					--> whole number that determines the minimum frequency plotted on the spectrogram 
	### maxFreq: int 					--> whole number that determines the maximum frequency plotted on the spectrogram 
	### vmaxContrast: float or int 		--> Denominator This will help with color contrast if desired!!!, e.g. 1.5 or 3
	### savePath: str   				--> Path to directory where fig should be saved; DEFAULT: '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/popContribFigs/'
	### saveFig: bool 					--> DEFAULT: True

 
	# Get relevant pop
	if type(pop) is str:
		popToPlot = pop
	elif type(pop) is list:
		popToPlot = pop[0]

	# Create figure 
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figSize)

	# Set font sizes
	labelFontSize = 20	# 12
	titleFontSize = 25	# 20

	#### SPECTROGRAM CALCULATIONS ####-------------------------------------------------------
	if 'spectrogram' in plotTypes:

		if hasBefore and hasAfter:
			allSignal = spectDict['allSignalFull']
			allFreqs = spectDict['allFreqsFull']
			timeRange = spectDict['T_full']
			vc = spectDict['vcFull']
		else:
			allSignal = spectDict['allSignalDuring']
			allFreqs = spectDict['allFreqsDuring']
			timeRange = spectDict['T_during']
			vc = spectDict['vcDuring']



		# Set color contrast parameters (vmin / vmax)
		orig_vmin = vc[0]
		orig_vmax = vc[1]

		if vmaxContrast is None:
			vmin = orig_vmin 	# np.amin(imshowSignal)		# None
			vmax = orig_vmax	# np.amax(imshowSignal)		# None
		else:
			vmin = orig_vmin 
			vmax = orig_vmax / vmaxContrast 
			vc = [vmin, vmax]


		# Set frequencies to be plotted (minFreq / maxFreq)
		if minFreq is None:
			minFreq = np.amin(allFreqs[0])			# DEFAULT: 1
		if maxFreq is None:
			maxFreq = np.amax(allFreqs[0])			# DEFAULT: 100
		elif maxFreq is not None:					# maxFreq must be an integer!! 
			if type(maxFreq) is not int:
				maxFreq = round(maxFreq)

		# Set up imshowSignal
		imshowSignal = allSignal[0][minFreq:maxFreq+1]


	#### HISTOGRAM CALCULATIONS ####-------------------------------------------------------
	if 'histogram' in plotTypes:
		# histoT = histDict['histoT']
		histoTBefore = histDict['histoTBefore']
		histoTDuring = histDict['histoTDuring']
		histoTAfter = histDict['histoTAfter']

		# histoCount = histDict['histoData']
		histoCountBefore = histDict['histoDataBefore']
		histoCountDuring = histDict['histoDataDuring']
		histoCountAfter = histDict['histoDataAfter']

		# get alignoffset
		alignoffset = histDict['alignoffset']



	#### PLOTTING ####-------------------------------------------------------
	# For filename naming, if/when saving figure --> 
	if hasBefore and hasAfter:
		figSaveStr='_full_'
	else:
		figSaveStr='_duringOsc_'

	# Plot both spectrogram and histogram -----------------------------------
	if 'spectrogram' in plotTypes and 'histogram' in plotTypes:
		# Plot Spectrogram 
		ax1 = plt.subplot(211)
		img = ax1.imshow(imshowSignal, extent=(np.amin(timeRange), np.amax(timeRange), minFreq, maxFreq), origin='lower', 
				interpolation='None', aspect='auto', cmap=plt.get_cmap(colorMap), vmin=vc[0], vmax=vc[1])
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad = 0.2)
		cax1.tick_params(labelsize=15) #12)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar to be in scientific notation
		fmt.set_powerlimits((0,0))
		cbar=plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		cbar.ax.set_ylabel('Power', size=15)
		ax1.set_title('Spike Rate Spectrogram for ' + popToPlot, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		## Set the spectrogram tick params (font size)
		ax1.tick_params(axis='both', labelsize=15) #12)

		# Plot Histogram 
		ax2 = plt.subplot(212)
		if hasBefore and hasAfter:
			ax2.bar(histoTBefore, histoCountBefore[0], width = 5, color='k', fill=False)
			ax2.bar(histoTDuring, histoCountDuring[0], width = 5, color=colorDict[popToPlot], fill=True)
			ax2.bar(histoTAfter, histoCountAfter[0], width = 5, color='k', fill=False)
		else:
			ax2.bar(histoTDuring, histoCountDuring[0], width = 5, color=colorDict[popToPlot], fill=True)
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes('right', size='3%', pad = 0.2)
		cax2.axis('off')
		ax2.set_title('Spike Rate Histogram for ' + popToPlot, fontsize=titleFontSize)
		ax2.set_xlabel('Time (ms)', fontsize=labelFontSize)
		ax2.set_ylabel('Rate (Hz)', fontsize=labelFontSize) # CLARIFY Y AXIS
		# Set the histogram tick params (font size)
		ax2.tick_params(axis='both', labelsize=15) #12)
		if hasBefore and hasAfter:
			T_before = histDict['T_before']
			T_after = histDict['T_after']
			ax2.set_xlim(left=T_before[0], right=T_after[1])
		else:
			T_during = histDict['T_during']
			ax2.set_xlim(left=T_during[0], right=T_during[1])

		# For figure saving
		figFilename = popToPlot + figSaveStr + '_combinedSpike.png'


	# Plot only spectrogram -----------------------------------------
	elif 'spectrogram' in plotTypes and 'histogram' not in plotTypes:
		# Plot Spectrogram 
		ax1 = plt.subplot(111)
		img = ax1.imshow(imshowSignal, extent=(np.amin(timeRange), np.amax(timeRange), minFreq, maxFreq), origin='lower', 
				interpolation='None', aspect='auto', cmap=plt.get_cmap(colorMap), vmin=vc[0], vmax=vc[0])
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad = 0.2)
		cax1.tick_params(labelsize=15) #12)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar to be in scientific notation
		fmt.set_powerlimits((0,0))
		cbar=plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		cbar.ax.set_ylabel('Power', size=15)
		ax1.set_title('Spike Rate Spectrogram for ' + popToPlot, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq)
		## Set the spectrogram tick params (font size)
		ax1.tick_params(axis='both', labelsize=15) #12)

		# For figure saving
		figFilename = popToPlot + figSaveStr + '_spike_spectrogram.png'


	# Plot only histogram -------------------------------------------
	elif 'spectrogram' not in plotTypes and 'histogram' in plotTypes:
		# Plot Histogram 
		ax2 = plt.subplot(111)
		if hasBefore and hasAfter:
			ax2.bar(histoTBefore, histoCountBefore[0], width = 5, color='k', fill=False)
			ax2.bar(histoTDuring, histoCountDuring[0], width = 5, color=colorDict[popToPlot], fill=True)
			ax2.bar(histoTAfter, histoCountAfter[0], width = 5, color='k', fill=False)
		else:
			ax2.bar(histoTDuring, histoCountDuring[0], width = 5, color=colorDict[popToPlot], fill=True)
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes('right', size='3%', pad = 0.2)
		cax2.axis('off')
		ax2.set_title('Spike Rate Histogram for ' + popToPlot, fontsize=titleFontSize)
		ax2.set_xlabel('Time (ms)', fontsize=labelFontSize)
		ax2.set_ylabel('Rate (Hz)', fontsize=labelFontSize) # CLARIFY Y AXIS
		# ax2.set_xlim(left=timeRange[0], right=timeRange[1])
		# Set the histogram tick params (font size)
		ax2.tick_params(axis='both', labelsize=15) #12)
		if hasBefore and hasAfter:
			T_before = histDict['T_before']
			T_after = histDict['T_after']
			ax2.set_xlim(left=T_before[0], right=T_after[1])
		else:
			T_during = histDict['T_during']
			ax2.set_xlim(left=T_during[0], right=T_during[1])

		# For potential figure saving
		figFilename = popToPlot + figSaveStr + '_spike_histogram.png'

	plt.tight_layout()
def plotSpikePanels(based):
	colorDictCustom = {}
	colorDictCustom['ITP4'] = 'magenta'
	colorDictCustom['ITS4'] = 'dodgerblue'
	colorDictCustom['IT5A'] = 'lime'

	timeRange, dataFile, channel = getWaveletInfo('theta', based=based, verbose=1)

	### OSC EVENT INFO DICTS !!
	thetaOscEventInfo = {'chan': 8, 'minT': 2785.22321038684, 
				'maxT': 3347.9278996316607, 'alignoffset':-3086.95, 'left': 55704, 'right':66958,
				'w2': 3376} 

	includePops=['ITS4', 'ITP4', 'IT5A']

	for pop in includePops:
		## Get dictionaries with spiking data for spectrogram and histogram plotting 
		spikeSpectDict = getSpikeData(dataFile, graphType='spect', pop=pop, oscEventInfo=thetaOscEventInfo)
		histDict = getSpikeData(dataFile, graphType='hist', pop=pop, oscEventInfo=thetaOscEventInfo)

		## Then call plotting function 
		if pop=='ITS4':
			popFigSize=(10,7)
		elif pop=='ITP4':
			popFigSize=(9.8,7)
		elif pop=='IT5A':
			popFigSize=(9.5,7)
		
		plotCombinedSpike(spectDict=spikeSpectDict, histDict=histDict, colorDict=colorDictCustom, plotTypes=['spectrogram', 'histogram'],
		hasBefore=1, hasAfter=1, pop=pop, figSize=popFigSize, colorMap='jet', vmaxContrast=2, maxFreq=None, saveFig=0) # figSize=(10,7)


	plt.show()



####### DATA DIRECTORY ####### 
based = '../data/v34_batch57/'  # Change this if necessary with path to data dir 



# --------------------------
# Main
# --------------------------
if __name__ == '__main__':
	# # Fig 8 -- heatmap
	# dfCSDPeak, dfCSDAvg, peakValues, avgValues, csdPopData= 
	plotHeatmap(based, absolute=0)
	# # Fig 8 -- CSD panels 
	# plotCSDPanels(based)
	# # Fig 8 -- Spike data panels 
	# plotSpikePanels(based)




