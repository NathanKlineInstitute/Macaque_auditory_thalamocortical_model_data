"""
simDataAnalysis.py

Functions to produce wavelet and CSD/LFP plotting for oscillation events in A1 model, using NetPyNE 

Contributors: ericaygriffith@gmail.com 
"""

### IMPORTS ### 
from netpyne import sim
import os
import utils
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker  ## for colorbar 
import matplotlib.image as mpimg
from PIL import Image
import numpy as np
import netpyne
#from netpyne.analysis import csd
from csd import *
from numbers import Number
import seaborn as sns 
import pandas as pd 
import pickle
import morlet
from morlet import MorletSpec, index2ms
from mpl_toolkits.axes_grid1 import make_axes_locatable
from loadSelect import * 
#from loadSelect2 import * 


## for COLORS in DATA PLOTTING!! 
colorList = [[0.42,0.67,0.84], [0.90,0.76,0.00], [0.42,0.83,0.59], [0.90,0.32,0.00],
			[0.34,0.67,0.67], [0.90,0.59,0.00], [0.42,0.82,0.83], [1.00,0.85,0.00],
			[0.33,0.67,0.47], [1.00,0.38,0.60], [0.57,0.67,0.33], [0.5,0.2,0.0],
			[0.71,0.82,0.41], [0.0,0.2,0.5], [0.70,0.32,0.10]]*3


######################################################################
##### FUNCTIONS IN PROGRESS ##### 
def evalPopsRelative():
	includePopsRel = []
	return includePopsRel
def evalPopsAbsolute():
	## CONDITIONS: 
	## (1) lfp signal avg'ed over all electrodes for each pop
	## (2) lfp signal for particular electrodes for each pop 
	## ADDRESS (2) LATER!!

	includePopsAbs = []
	return includePopsAbs
def evalWaveletsByBand(based, dfPklFile):
	## NOT IN USE RIGHT NOW --> ## freqBand: str 			--> e.g. 'alpha', 'beta' ,'theta', 'delta', 'gamma'
	## NOT IN USE RIGHT NOW --> dlmsPklFile: .pkl file 	--> from dlms.pkl file, saved from load.py 
	## based: str 				--> Beginning of path to the .pkl files 
	## dfPklFile: .pkl file 	--> from df.pkl file, saved from load.py 

	print('Evaluating all oscillation events in a given frequency band')

	if based is None:
		based = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/wavelets/sim/spont/' ### COULD make this an arg!! # '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/wavelets/' 
	subjectDir = dfPklFile.split('_df.pkl')[0]

	# Load df file 
	dfFullPath = based + subjectDir + '/' + dfPklFile
	df = pd.read_pickle(dfFullPath)

	# # Load dlms file
	# dlmsFullPath = based + subjectDir + '/' + dlmsPklFile
	# dlmsFile = open(dlmsFullPath, 'rb')
	# dlms = pickle.load(dlmsFile)
	# dlmsFile.close()

	return df
## plotting functions
def plotTimeSeries(timeSeriesDict, dataType=['LFP', 'CSD']):
	### dataType: list of str 		--> Indicates which dataType will be input / plotted
	### timeSeriesDict: dict 		--> Output of ____ (which exact funtions)
	print('Plot time series for LFP or CSD data')
def plotSpectrogram(spectDict, dataType=['LFP', 'CSD', 'Spike']):
	### dataType: list of str 		--> Indicates which dataType will be input / plotted
	### spectDict: dict 			--> output of ____ which functions?
	print('Plot spectrogram for LFP, CSD, or Spiking data')
def plotHistogram(histDict):
	### histDict: dict 			--> output of ___which function exactly? 
	print('Plot histogram of spiking data')
######################################################################




#################################################################################################################
######### NetPyNE Functions that have been modified for use here!! ##############################################
#################################################################################################################

### LFP ### 
def getLFPData(pop=None, timeRange=None, electrodes=['avg', 'all'], plots=['timeSeries', 'spectrogram'], inputLFP=None, NFFT=256, noverlap=128, nperseg=256, 
	minFreq=1, maxFreq=100, stepFreq=1, smooth=0, separation=1.0, logx=False, logy=False, normSignal=False, normSpec=False, filtFreq=False, filtOrder=3, detrend=False, 
	transformMethod='morlet'):
	"""
	Function for plotting local field potentials (LFP)

	Parameters

	pop: str (NOTE: for now) 
		Population to plot lfp data for (sim.allSimData['LFPPops'][pop] <-- requires LFP pop saving!)
		``None`` plots overall LFP data 

	timeRange : list [start, stop]
		Time range to plot.
		**Default:**
		``None`` plots entire time range

	electrodes : list
		List of electrodes to include; ``'avg'`` is the average of all electrodes; ``'all'`` is each electrode separately.
		**Default:** ``['avg', 'all']``

	plots : list
		List of plot types to show.
		**Default:** ``['timeSeries', 'spectrogram']`` <-- added 'PSD' (EYG, 2/10/2022) <-- 'PSD' taken out on 3/16/22

	NFFT : int (power of 2)
		Number of data points used in each block for the PSD and time-freq FFT.
		**Default:** ``256``

	noverlap : int (<nperseg)
		Number of points of overlap between segments for PSD and time-freq.
		**Default:** ``128``

	nperseg : int
		Length of each segment for time-freq.
		**Default:** ``256``

	minFreq : float
		Minimum frequency shown in plot for PSD and time-freq.
		**Default:** ``1``

	maxFreq : float
		Maximum frequency shown in plot for PSD and time-freq.
		**Default:** ``100``

	stepFreq : float
		Step frequency.
		**Default:** ``1``

	smooth : int
		Window size for smoothing LFP; no smoothing if ``0``
		**Default:** ``0``

	separation : float
		Separation factor between time-resolved LFP plots; multiplied by max LFP value.
		**Default:** ``1.0``

	logx : bool
		Whether to make x-axis logarithmic
		**Default:** ``False``
		**Options:** ``<option>`` <description of option>

	logy : bool
		Whether to make y-axis logarithmic
		**Default:** ``False``

	normSignal : bool
		Whether to normalize the signal.
		**Default:** ``False``

	normSpec : bool
		Needs documentation.
		**Default:** ``False``

	filtFreq : int or list
		Frequency for low-pass filter (int) or frequencies for bandpass filter in a list: [low, high]
		**Default:** ``False`` does not filter the data

	filtOrder : int
		Order of the filter defined by `filtFreq`.
		**Default:** ``3``

	detrend : bool
		Whether to detrend.
		**Default:** ``False``

	transformMethod : str
		Transform method.
		**Default:** ``'morlet'``
		**Options:** ``'fft'``

	Returns
	-------
	(dict)
		A tuple consisting of the Matplotlib figure handles and a dictionary containing the plot data
	"""

	# from .. import sim   ### <-- should already be there!!!! in imports!! 

	print('Getting LFP data ...')


	# time range
	if timeRange is None:
		timeRange = [0,sim.cfg.duration]

	# populations
	#### CLEAN THIS UP.... go through all the possibilities implied by these if statements and make sure all are accounted for! 
	if inputLFP is None: 
		if pop is None:
			if timeRange is None:
				lfp = np.array(sim.allSimData['LFP'])
			elif timeRange is not None: 
				lfp = np.array(sim.allSimData['LFP'])[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]

		elif pop is not None:
			if type(pop) is str:
				popToPlot = pop
			elif type(pop) is list and len(pop)==1:
				popToPlot = pop[0]
			# elif type(pop) is list and len(pop) > 1:  #### USE THIS AS JUMPING OFF POINT TO EXPAND FOR LIST OF MULTIPLE POPS?? 
			lfp = np.array(sim.allSimData['LFPPops'][popToPlot])[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]

	#### MAKE SURE THAT THIS ADDITION DOESN'T MAKE ANYTHING ELSE BREAK!! 
	elif inputLFP is not None:  ### DOING THIS FOR PSD FOR SUMMED LFP SIGNAL !!! 
		lfp = inputLFP 

		# if timeRange is None:
		#     lfp = inputLFP 
		# elif timeRange is not None:
		#     lfp = inputLFP[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:] ### hmm. 


	if filtFreq:
		from scipy import signal
		fs = 1000.0/sim.cfg.recordStep
		nyquist = fs/2.0
		if isinstance(filtFreq, list): # bandpass
			Wn = [filtFreq[0]/nyquist, filtFreq[1]/nyquist]
			b, a = signal.butter(filtOrder, Wn, btype='bandpass')
		elif isinstance(filtFreq, Number): # lowpass
			Wn = filtFreq/nyquist
			b, a = signal.butter(filtOrder, Wn)
		for i in range(lfp.shape[1]):
			lfp[:,i] = signal.filtfilt(b, a, lfp[:,i])

	if detrend:
		from scipy import signal
		for i in range(lfp.shape[1]):
			lfp[:,i] = signal.detrend(lfp[:,i])

	if normSignal:
		for i in range(lfp.shape[1]):
			offset = min(lfp[:,i])
			if offset <= 0:
				lfp[:,i] += abs(offset)
			lfp[:,i] /= max(lfp[:,i])

	# electrode selection
	print('electrodes: ' + str(electrodes))
	if electrodes is None:
		print('electrodes is None -- improve this')
	elif type(electrodes) is list:
		if 'all' in electrodes:
			electrodes.remove('all')
			electrodes.extend(list(range(int(sim.net.recXElectrode.nsites))))

	data = {'lfp': lfp}  # returned data


	# time series -----------------------------------------
	if 'timeSeries' in plots:
		ydisp = np.absolute(lfp).max() * separation
		offset = 1.0*ydisp
		t = np.arange(timeRange[0], timeRange[1], sim.cfg.recordStep)


		for i,elec in enumerate(electrodes):
			if elec == 'avg':
				lfpPlot = np.mean(lfp, axis=1)
				color = 'k'
				lw=1.0
			elif isinstance(elec, Number) and (inputLFP is not None or elec <= sim.net.recXElectrode.nsites):
				lfpPlot = lfp[:, elec]
				color = 'k' #colors[i%len(colors)]
				lw = 1.0

			if len(t) < len(lfpPlot):
				lfpPlot = lfpPlot[:len(t)]


		data['lfpPlot'] = lfpPlot
		data['ydisp'] =  ydisp
		data['t'] = t

	# Spectrogram ------------------------------
	if 'spectrogram' in plots:
		import matplotlib.cm as cm
		numCols = 1 #np.round(len(electrodes) / maxPlots) + 1

		# Morlet wavelet transform method
		if transformMethod == 'morlet':
			spec = []
			freqList = None
			if logy:
				freqList = np.logspace(np.log10(minFreq), np.log10(maxFreq), int((maxFreq-minFreq)/stepFreq))

			for i,elec in enumerate(electrodes):
				if elec == 'avg':
					lfpPlot = np.mean(lfp, axis=1)
				elif isinstance(elec, Number) and (inputLFP is not None or elec <= sim.net.recXElectrode.nsites):
					lfpPlot = lfp[:, elec]
				fs = int(1000.0 / sim.cfg.recordStep)
				t_spec = np.linspace(0, morlet.index2ms(len(lfpPlot), fs), len(lfpPlot))
				spec.append(MorletSpec(lfpPlot, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq, lfreq=freqList))

			f = freqList if freqList is not None else np.array(range(minFreq, maxFreq+1, stepFreq))   # only used as output for user

			vmin = np.array([s.TFR for s in spec]).min()
			vmax = np.array([s.TFR for s in spec]).max()

			for i,elec in enumerate(electrodes):
				T = timeRange
				F = spec[i].f
				if normSpec:
					spec[i].TFR = spec[i].TFR / vmax
					S = spec[i].TFR
					vc = [0, 1]
				else:
					S = spec[i].TFR
					vc = [vmin, vmax]


		# FFT transform method
		elif transformMethod == 'fft':

			from scipy import signal as spsig
			spec = []

			for i,elec in enumerate(electrodes):
				if elec == 'avg':
					lfpPlot = np.mean(lfp, axis=1)
				elif isinstance(elec, Number) and elec <= sim.net.recXElectrode.nsites:
					lfpPlot = lfp[:, elec]
				# creates spectrogram over a range of data
				# from: http://joelyancey.com/lfp-python-practice/
				fs = int(1000.0/sim.cfg.recordStep)
				f, t_spec, x_spec = spsig.spectrogram(lfpPlot, fs=fs, window='hanning',
				detrend=mlab.detrend_none, nperseg=nperseg, noverlap=noverlap, nfft=NFFT,  mode='psd')
				x_mesh, y_mesh = np.meshgrid(t_spec*1000.0, f[f<maxFreq])
				spec.append(10*np.log10(x_spec[f<maxFreq]))

			vmin = np.array(spec).min()
			vmax = np.array(spec).max()


	outputData = {'LFP': lfp, 'lfpPlot': lfpPlot, 'electrodes': electrodes, 'timeRange': timeRange}
	### Added lfpPlot to this, because that usually has the post-processed electrode-based lfp data

	if 'timeSeries' in plots:
		outputData.update({'t': t})

	if 'spectrogram' in plots:
		outputData.update({'spec': spec, 't': t_spec*1000.0, 'freqs': f[f<=maxFreq]})

	# if 'PSD' in plots:
	# 	outputData.update({'allFreqs': allFreqs, 'allSignal': allSignal})


	return outputData
def plotLFP(pop=None, timeRange=None, electrodes=['avg', 'all'], plots=['timeSeries', 'PSD', 'spectrogram', 'locations'], inputLFP=None, NFFT=256, noverlap=128, nperseg=256, minFreq=1, maxFreq=100, stepFreq=1, smooth=0, separation=1.0, includeAxon=True, logx=False, logy=False, normSignal=False, normPSD=False, normSpec=False, filtFreq=False, filtOrder=3, detrend=False, transformMethod='morlet', maxPlots=8, overlay=False, colors=None, figSize=(8, 8), fontSize=14, lineWidth=1.5, dpi=200, saveData=None, saveFig=None, showFig=True):
	"""
	Parameters

	pop: str (NOTE: for now) 
		Population to plot lfp data for (sim.allSimData['LFPPops'][pop] <-- requires LFP pop saving!)
		``None`` plots overall LFP data 

	timeRange : list [start, stop]
		Time range to plot.
		**Default:**
		``None`` plots entire time range

	electrodes : list
		List of electrodes to include; ``'avg'`` is the average of all electrodes; ``'all'`` is each electrode separately.
		**Default:** ``['avg', 'all']``

	plots : list
		List of plot types to show.
		**Default:** ``['timeSeries', 'PSD', 'spectrogram', 'locations']``

	NFFT : int (power of 2)
		Number of data points used in each block for the PSD and time-freq FFT.
		**Default:** ``256``

	noverlap : int (<nperseg)
		Number of points of overlap between segments for PSD and time-freq.
		**Default:** ``128``

	nperseg : int
		Length of each segment for time-freq.
		**Default:** ``256``

	minFreq : float
		Minimum frequency shown in plot for PSD and time-freq.
		**Default:** ``1``

	maxFreq : float
		Maximum frequency shown in plot for PSD and time-freq.
		**Default:** ``100``

	stepFreq : float
		Step frequency.
		**Default:** ``1``

	smooth : int
		Window size for smoothing LFP; no smoothing if ``0``
		**Default:** ``0``

	separation : float
		Separation factor between time-resolved LFP plots; multiplied by max LFP value.
		**Default:** ``1.0``

	includeAxon : bool
		Whether to show the axon in the location plot.
		**Default:** ``True``

	logx : bool
		Whether to make x-axis logarithmic
		**Default:** ``False``
		**Options:** ``<option>`` <description of option>

	logy : bool
		Whether to make y-axis logarithmic
		**Default:** ``False``

	normSignal : bool
		Whether to normalize the signal.
		**Default:** ``False``

	normPSD : bool
		Whether to normalize the power spectral density.
		**Default:** ``False``

	normSpec : bool
		Needs documentation.
		**Default:** ``False``

	filtFreq : int or list
		Frequency for low-pass filter (int) or frequencies for bandpass filter in a list: [low, high]
		**Default:** ``False`` does not filter the data

	filtOrder : int
		Order of the filter defined by `filtFreq`.
		**Default:** ``3``

	detrend : bool
		Whether to detrend.
		**Default:** ``False``

	transformMethod : str
		Transform method.
		**Default:** ``'morlet'``
		**Options:** ``'fft'``

	maxPlots : int
		Maximum number of subplots. Currently unused.
		**Default:** ``8``

	overlay : bool
		Whether to overlay plots or use subplots.
		**Default:** ``False`` overlays plots.

	colors : list
		List of normalized RGB colors to use.
		**Default:** ``None`` uses standard colors

	figSize : list [width, height]
		Size of figure in inches.
		**Default:** ``(10, 8)``

	fontSize : int
		Font size on figure.
		**Default:** ``14``

	lineWidth : int
		Line width.
		**Default:** ``1.5``

	dpi : int
		Resolution of figure in dots per inch.
		**Default:** ``100``

	saveData : bool or str
		Whether and where to save the data used to generate the plot.
		**Default:** ``False``
		**Options:** ``True`` autosaves the data,
		``'/path/filename.ext'`` saves to a custom path and filename, valid file extensions are ``'.pkl'`` and ``'.json'``

	saveFig : bool or str
		Whether and where to save the figure.
		**Default:** ``False``
		**Options:** ``True`` autosaves the figure,
		``'/path/filename.ext'`` saves to a custom path and filename, valid file extensions are ``'.png'``, ``'.jpg'``, ``'.eps'``, and ``'.tiff'``

	showFig : bool
		Shows the figure if ``True``.
		**Default:** ``True``

	Returns
	-------
	(figs, dict)
		A tuple consisting of the Matplotlib figure handles and a dictionary containing the plot data

	"""

	# from .. import sim  ### <-- should already be here 
	# from ..support.scalebar import add_scalebar
	import testScalebar

	print('Plotting LFP ...')

	if not colors: colors = colorList

	# set font size
	plt.rcParams.update({'font.size': fontSize})

	# time range
	if timeRange is None:
		timeRange = [0,sim.cfg.duration]

	# populations
	if pop is None:
		if inputLFP is not None:
			lfp = inputLFP #[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]
		else:
			lfp = np.array(sim.allSimData['LFP'])[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]
	elif pop is not None:
		lfp = np.array(sim.allSimData['LFPPops'][pop])[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]

	print('lfp shape:' + str(lfp.shape))


	if filtFreq:
		from scipy import signal
		fs = 1000.0/sim.cfg.recordStep
		nyquist = fs/2.0
		if isinstance(filtFreq, list): # bandpass
			Wn = [filtFreq[0]/nyquist, filtFreq[1]/nyquist]
			b, a = signal.butter(filtOrder, Wn, btype='bandpass')
		elif isinstance(filtFreq, Number): # lowpass
			Wn = filtFreq/nyquist
			b, a = signal.butter(filtOrder, Wn)
		for i in range(lfp.shape[1]):
			lfp[:,i] = signal.filtfilt(b, a, lfp[:,i])

	if detrend:
		from scipy import signal
		for i in range(lfp.shape[1]):
			lfp[:,i] = signal.detrend(lfp[:,i])

	if normSignal:
		for i in range(lfp.shape[1]):
			offset = min(lfp[:,i])
			if offset <= 0:
				lfp[:,i] += abs(offset)
			lfp[:,i] /= max(lfp[:,i])

	# electrode selection
	if electrodes is None:
		print('electrodes is None: improve this -- ') ### FOR PSD PLOTTING / INFO FOR SUMMED LFP SIGNAL!!
	elif type(electrodes) is list:
		if 'all' in electrodes:
			electrodes.remove('all')
			electrodes.extend(list(range(int(sim.net.recXElectrode.nsites))))

	# plotting
	figs = []

	data = {'lfp': lfp}  # returned data


	# time series -----------------------------------------
	if 'timeSeries' in plots:
		ydisp = np.absolute(lfp).max() * separation
		offset = 1.0*ydisp
		t = np.arange(timeRange[0], timeRange[1], sim.cfg.recordStep)

		if figSize:
			figs.append(plt.figure(figsize=figSize))

		for i,elec in enumerate(electrodes):
			if elec == 'avg':
				lfpPlot = np.mean(lfp, axis=1)
				color = 'k'
				lw=1.0
			elif isinstance(elec, Number) and (inputLFP is not None or elec <= sim.net.recXElectrode.nsites):
				lfpPlot = lfp[:, elec]
				color = colors[i%len(colors)]
				lw = 1.0

			if len(t) < len(lfpPlot):
				lfpPlot = lfpPlot[:len(t)]

			plt.plot(t[0:len(lfpPlot)], -lfpPlot+(i*ydisp), color=color, linewidth=lw)
			if len(electrodes) > 1:
				plt.text(timeRange[0]-0.07*(timeRange[1]-timeRange[0]), (i*ydisp), elec, color=color, ha='center', va='top', fontsize=fontSize, fontweight='bold')

		ax = plt.gca()

		data['lfpPlot'] = lfpPlot
		data['ydisp'] =  ydisp
		data['t'] = t

		# format plot
		if len(electrodes) > 1:
			plt.text(timeRange[0]-0.14*(timeRange[1]-timeRange[0]), (len(electrodes)*ydisp)/2.0, 'LFP electrode', color='k', ha='left', va='bottom', fontSize=fontSize, rotation=90)
			plt.ylim(-offset, (len(electrodes))*ydisp)
		else:
			if pop is None:
				plt.suptitle('LFP Signal', fontSize=fontSize, fontweight='bold')
			elif pop is not None:
				timeSeriesTitle = 'LFP Signal of ' + pop + ' population'
				plt.suptitle(timeSeriesTitle, fontSize=fontSize, fontweight='bold')
		ax.invert_yaxis()
		plt.xlabel('time (ms)', fontsize=fontSize)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['left'].set_visible(False)
		plt.subplots_adjust(bottom=0.1, top=1.0, right=1.0)

		# calculate scalebar size and add scalebar
		round_to_n = lambda x, n, m: int(np.ceil(round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1)) / m)) * m
		scaley = 1000.0  # values in mV but want to convert to uV
		m = 10.0
		sizey = 100/scaley
		while sizey > 0.25*ydisp:
			try:
				sizey = round_to_n(0.2*ydisp*scaley, 1, m) / scaley
			except:
				sizey /= 10.0
			m /= 10.0
		labely = '%.3g $\mu$V'%(sizey*scaley)#)[1:]
		if len(electrodes) > 1:
			add_scalebar(ax,hidey=True, matchy=False, hidex=False, matchx=False, sizex=0, sizey=-sizey, labely=labely, unitsy='$\mu$V', scaley=scaley, loc=3, pad=0.5, borderpad=0.5, sep=3, prop=None, barcolor="black", barwidth=2)
		else:
			add_scalebar(ax, hidey=True, matchy=False, hidex=True, matchx=True, sizex=None, sizey=-sizey, labely=labely, unitsy='$\mu$V', scaley=scaley, unitsx='ms', loc=3, pad=0.5, borderpad=0.5, sep=3, prop=None, barcolor="black", barwidth=2)
		# save figure
		if saveFig:
			if isinstance(saveFig, basestring):
				filename = saveFig
			else:
				filename = sim.cfg.filename + '_LFP_timeseries.png'
			plt.savefig(filename, dpi=dpi)

	# Spectrogram ------------------------------
	if 'spectrogram' in plots:
		import matplotlib.cm as cm
		numCols = 1 #np.round(len(electrodes) / maxPlots) + 1
		figs.append(plt.figure(figsize=(figSize[0]*numCols, figSize[1])))

		# Morlet wavelet transform method
		if transformMethod == 'morlet':
			spec = []
			freqList = None
			if logy:
				freqList = np.logspace(np.log10(minFreq), np.log10(maxFreq), int((maxFreq-minFreq)/stepFreq))

			for i,elec in enumerate(electrodes):
				if elec == 'avg':
					lfpPlot = np.mean(lfp, axis=1)
				elif isinstance(elec, Number) and (inputLFP is not None or elec <= sim.net.recXElectrode.nsites):
					lfpPlot = lfp[:, elec]
				fs = int(1000.0 / sim.cfg.recordStep)
				t_spec = np.linspace(0, morlet.index2ms(len(lfpPlot), fs), len(lfpPlot))
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

				plt.imshow(S, extent=(np.amin(T), np.amax(T), np.amin(F), np.amax(F)), origin='lower', interpolation='None', aspect='auto', vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap('viridis'))
				if normSpec:
					plt.colorbar(label='Normalized power')
				else:
					plt.colorbar(label='Power')
				plt.ylabel('Hz')
				plt.tight_layout()
				if len(electrodes) > 1:
					plt.title('Electrode %s' % (str(elec)), fontsize=fontSize - 2)

		# FFT transform method
		elif transformMethod == 'fft':

			from scipy import signal as spsig
			spec = []

			for i,elec in enumerate(electrodes):
				if elec == 'avg':
					lfpPlot = np.mean(lfp, axis=1)
				elif isinstance(elec, Number) and elec <= sim.net.recXElectrode.nsites:
					lfpPlot = lfp[:, elec]
				# creates spectrogram over a range of data
				# from: http://joelyancey.com/lfp-python-practice/
				fs = int(1000.0/sim.cfg.recordStep)
				f, t_spec, x_spec = spsig.spectrogram(lfpPlot, fs=fs, window='hanning',
				detrend=mlab.detrend_none, nperseg=nperseg, noverlap=noverlap, nfft=NFFT,  mode='psd')
				x_mesh, y_mesh = np.meshgrid(t_spec*1000.0, f[f<maxFreq])
				spec.append(10*np.log10(x_spec[f<maxFreq]))

			vmin = np.array(spec).min()
			vmax = np.array(spec).max()

			for i,elec in enumerate(electrodes):
				plt.subplot(np.ceil(len(electrodes)/numCols), numCols, i+1)
				plt.pcolormesh(x_mesh, y_mesh, spec[i], cmap=cm.viridis, vmin=vmin, vmax=vmax)
				plt.colorbar(label='dB/Hz', ticks=[np.ceil(vmin), np.floor(vmax)])
				if logy:
					plt.yscale('log')
					plt.ylabel('Log-frequency (Hz)')
					if isinstance(logy, list):
						yticks = tuple(logy)
						plt.yticks(yticks, yticks)
				else:
					plt.ylabel('(Hz)')
				if len(electrodes) > 1:
					plt.title('Electrode %s'%(str(elec)), fontsize=fontSize-2)

		plt.xlabel('time (ms)', fontsize=fontSize)
		plt.tight_layout()
		if pop is None:
			plt.suptitle('LFP spectrogram', size=fontSize, fontweight='bold')
		elif pop is not None:
			spectTitle = 'LFP spectrogram of ' + pop + ' population'
			plt.suptitle(spectTitle, size=fontSize, fontweight='bold')

		plt.subplots_adjust(bottom=0.08, top=0.90)

		# save figure
		if saveFig:
			if isinstance(saveFig, basestring):
				filename = saveFig
			else:
				filename = sim.cfg.filename + '_LFP_timefreq.png'
			plt.savefig(filename, dpi=dpi)

	# locations ------------------------------
	if 'locations' in plots:
		try:
			cvals = [] # used to store total transfer resistance

			for cell in sim.net.compartCells:
				trSegs = list(np.sum(sim.net.recXElectrode.getTransferResistance(cell.gid)*1e3, axis=0)) # convert from Mohm to kilohm
				if not includeAxon:
					i = 0
					for secName, sec in cell.secs.items():
						nseg = sec['hObj'].nseg #.geom.nseg
						if 'axon' in secName:
							for j in range(i,i+nseg): del trSegs[j]
						i+=nseg
				cvals.extend(trSegs)

			includePost = [c.gid for c in sim.net.compartCells]
			fig = sim.analysis.plotShape(includePost=includePost, showElectrodes=electrodes, cvals=cvals, includeAxon=includeAxon, dpi=dpi,
			fontSize=fontSize, saveFig=saveFig, showFig=showFig, figSize=figSize)[0]
			figs.append(fig)
		except:
			print('  Failed to plot LFP locations...')

	# PSD ----------------------------------
	if 'PSD' in plots:
		if overlay:
			figs.append(plt.figure(figsize=figSize))
		else:
			numCols = 1 # np.round(len(electrodes) / maxPlots) + 1
			figs.append(plt.figure(figsize=(figSize[0]*numCols, figSize[1])))
			#import seaborn as sb

		allFreqs = []
		allSignal = []
		data['allFreqs'] = allFreqs
		data['allSignal'] = allSignal

		for i,elec in enumerate(electrodes):
			print('elec: ' + str(elec))
			if elec == 'avg':
				lfpPlot = np.mean(lfp, axis=1)
			elif isinstance(elec, Number) and (inputLFP is not None or elec <= sim.net.recXElectrode.nsites):
				lfpPlot = lfp[:, elec]

			# Morlet wavelet transform method
			if transformMethod == 'morlet':
				# from morlet import MorletSpec, index2ms

				Fs = int(1000.0/sim.cfg.recordStep)

				#t_spec = np.linspace(0, index2ms(len(lfpPlot), Fs), len(lfpPlot))
				morletSpec = MorletSpec(lfpPlot, Fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq)
				freqs = F = morletSpec.f
				spec = morletSpec.TFR
				signal = np.mean(spec, 1)
				ylabel = 'Power'

			# FFT transform method
			elif transformMethod == 'fft':
				Fs = int(1000.0/sim.cfg.recordStep)
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
				plt.subplot(int(np.ceil(len(electrodes)/numCols)), numCols,i+1)
			if elec == 'avg':
				color = 'k'
			elif isinstance(elec, Number) and (inputLFP is not None or elec <= sim.net.recXElectrode.nsites):
				color = colors[i % len(colors)]
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


	outputData = {'LFP': lfp, 'electrodes': electrodes, 'timeRange': timeRange, 'saveData': saveData, 'saveFig': saveFig, 'showFig': showFig}

	if 'timeSeries' in plots:
		outputData.update({'t': t})

	if 'PSD' in plots:
		outputData.update({'allFreqs': allFreqs, 'allSignal': allSignal})

	if 'spectrogram' in plots:
		outputData.update({'spec': spec, 't': t_spec*1000.0, 'freqs': f[f<=maxFreq]})

	#save figure data
	if saveData:
		figData = outputData
		_saveFigData(figData, saveData, 'lfp')

	# show fig
	if showFig: plt.show() #_showFigure()


	return figs, outputData
#################################################################################################################



###################
#### FUNCTIONS ####
###################
def getWaveletInfo(freqBand, based, verbose=0): 
	## freqBand: str  --> e.g. 'delta', 'alpha', 'theta'
	## based: str --> path to directory with the .pkl data files 
	## verbose: bool --> if 0, default to only putting out timeRange and dataFile, if 1 --> include channel as well 

	waveletInfo = {
	'delta': {'dataFile': 'A1_v34_batch65_v34_batch65_0_0_data.pkl', 'timeRange': [1480, 2520], 'channel': 14},
	'beta': {'dataFile': 'A1_v34_batch67_v34_batch67_0_0_data.pkl',	'timeRange': [456, 572], 'channel': 14}, 
	'alpha': {'dataFile': 'A1_v34_batch67_v34_batch67_0_0_data.pkl', 'timeRange': [3111, 3325], 'channel': 9}, 
	'theta': {'dataFile': 'A1_v34_batch67_v34_batch67_1_1_data.pkl', 'timeRange': [2785, 3350], 'channel': 8}}

	timeRange = waveletInfo[freqBand]['timeRange']
	dataFileNoPath = waveletInfo[freqBand]['dataFile']
	dataFile = based + dataFileNoPath
	channel = waveletInfo[freqBand]['channel']

	if verbose:
		return timeRange, dataFile, channel
	else:
		return timeRange, dataFile

# getSimLayers and geteventprop are dependencies for getOscEventInfo
def getSimLayers():
	layers = {'supra':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 'gran':[10, 11], 'infra':[12, 13, 14, 15, 16, 17, 18, 19]} ## NETPYNE SIMS
	return layers
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
# ONLY FOR SIM SUBJECTS FOR NOW!! 
### TO DO: (1) expand to NHP (2) add in capability for 'all' regions 
### NOTE: can / should this even be expanded to NHP? 
def getOscEventInfo(frequencyBands, waveletPath, subjects=None):
	# subjects: list 
	# frequencyBands: list 
	# waveletPath: str 

	layers = getSimLayers()   # layers = {'supra':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 'gran':[10, 11], 'infra':[12, 13, 14, 15, 16, 17, 18, 19]} ## NETPYNE SIMS


	if subjects is None:
		os.chdir(waveletPath)
		dirContents = os.listdir()
		subjs = []
		for dirContent in dirContents:
			if os.path.isdir(dirContent) and 'batch' in dirContent:
				subjs.append(dirContent)
		subjects = subjs

	print('Analyzing subjects: ' + str(subjects))


	## organize a dict with all the wavelets from the above subjects:
	oscEventInfo = {}
	for band in frequencyBands:
		oscEventInfo[band] = {}

		for layerRegion in layers:
			oscEventInfo[band][layerRegion] = {}

			for subj in subjects:
				oscEventInfo[band][layerRegion][subj] = {}

				for chan in layers[layerRegion]:
					chanDir = 'chan_' + str(chan) + '/'

					dfsDir = waveletPath + subj  + '/' + chanDir
					# print('dfsDir aka chanDir: ' + str(dfsDir))  # TESTING LINE!! 
					os.chdir(dfsDir)  
					chanDirFiles = os.listdir()
					dfsFiles = []
					for chanDirFile in chanDirFiles:
						if '.pkl' in chanDirFile:
							dfsFiles.append(chanDirFile)   
					## so now here we have a list w/ the dfs.pkl files from a particular chan in a particular subject 

					for dfsFile in dfsFiles:
						# read this file for event idx info etc 
						dfsPkl = pd.read_pickle(dfsDir + dfsFile)

						for idx in dfsPkl.index:
							if band in dfsFile:
								oscEventInfo[band][layerRegion][subj][idx] = {}

								eventIdx = idx
								chan = dfsPkl.at[idx,'chan']
								minT = dfsPkl.at[idx,'minT']
								maxT = dfsPkl.at[idx,'maxT']
								WaveletPeakT = dfsPkl.at[idx,'WaveletPeakT']
								alignoffset = -WaveletPeakT
								left = dfsPkl.at[idx,'left']
								right = dfsPkl.at[idx,'right']
								## w2 calculations
								w2 = int((right-left+1)/2.)
								# Resize w2 to match the load.py calculation for the osc event plotting (in def draw() in class eventviewer)
								w2 = int(w2*0.6)
								## print('w2: ' + str(w2))

								# add in the above info into dict 
								oscEventInfo[band][layerRegion][subj][idx]['eventIdx'] = eventIdx
								oscEventInfo[band][layerRegion][subj][idx]['chan'] = chan
								oscEventInfo[band][layerRegion][subj][idx]['minT'] = minT
								oscEventInfo[band][layerRegion][subj][idx]['maxT'] = maxT
								oscEventInfo[band][layerRegion][subj][idx]['alignoffset'] = alignoffset
								oscEventInfo[band][layerRegion][subj][idx]['left'] = left
								oscEventInfo[band][layerRegion][subj][idx]['right'] = right
								oscEventInfo[band][layerRegion][subj][idx]['w2'] = w2 

	return oscEventInfo


## Heatmaps for LFP data ## 
def getDataFrames(dataFile, timeRange, verbose=0):  ### Make this work with arbitrary input data, not just LFP so can look at CSD as well!!!! 
	## This function will return data frames of peak and average LFP amplitudes, for picking cell pops
	### dataFile: str 		--> .pkl file to load, with data from the whole recording
	### timeRange: list 	--> e.g. [start, stop]
	### verbose: bool 		--> if 0, return only the data frames; if 1 - return all lists and dataframes 

	## Load data file
	sim.load(dataFile, instantiate=False)

	## Get all cell pops (cortical)
	thalPops = ['TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']
	allPops = list(sim.net.allPops.keys())
	pops = [pop for pop in allPops if pop not in thalPops] 			## exclude thal pops 

	## Get all electrodes 
	evalElecs = []
	evalElecs.extend(list(range(int(sim.net.recXElectrode.nsites))))
	### add 'avg' to electrode list 
	evalElecs.append('avg')


	## Create dict with lfp population data, and calculate average and peak amplitudes for each pop at each electrode!
	lfpPopData = {}

	for pop in pops:
		lfpPopData[pop] = {}

		popLFPdata = np.array(sim.allSimData['LFPPops'][pop])[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]

		for i, elec in enumerate(evalElecs):
			if elec == 'avg':
				lfpPopData[pop]['avg'] = {}

				avgPopData = np.mean(popLFPdata, axis=1) 		## lfp data (from 1 pop) at each timepoint, averaged over all electrodes

				avgAvgLFP = np.average(avgPopData)
				lfpPopData[pop]['avg']['avg'] = avgAvgLFP 		## time-average of lfp data (from 1 pop) that has been averaged in space (over all electrodes)

				peakAvgLFP = np.amax(avgPopData)
				lfpPopData[pop]['avg']['peak'] = peakAvgLFP 	## highest datapoint of all the lfp data (from 1 pop) that has been averaged in space (over all electrodes)

			elif isinstance(elec, Number):
				elecKey = 'elec' + str(elec)
				lfpPopData[pop][elecKey] = {}
				lfpPopData[pop][elecKey]['avg'] = np.average(popLFPdata[:, elec]) 	## LFP data (from 1 pop) averaged in time, over 1 electrode 
				lfpPopData[pop][elecKey]['peak'] = np.amax(popLFPdata[:, elec]) 	## Maximum LFP value from 1 pop, over time, recorded at 1 electrode 


	#### PEAK LFP AMPLITUDES, DATA FRAME ####
	peakValues = {}
	peakValues['pops'] = []
	peakValues['peakLFP'] = [[] for i in range(len(pops))]  # should be 36! 
	p=0
	for pop in pops:
		peakValues['pops'].append(pop)
		for i, elec in enumerate(evalElecs): 
			if isinstance(elec, Number):
				elecKey = 'elec' + str(elec)
			elif elec == 'avg': 
				elecKey = 'avg'
			peakValues['peakLFP'][p].append(lfpPopData[pop][elecKey]['peak'])
		p+=1
	dfPeak = pd.DataFrame(peakValues['peakLFP'], index=pops)


	#### AVERAGE LFP AMPLITUDES, DATA FRAME ####
	avgValues = {}
	avgValues['pops'] = []
	avgValues['avgLFP'] = [[] for i in range(len(pops))]
	q=0
	for pop in pops:
		avgValues['pops'].append(pop)
		for i, elec in enumerate(evalElecs):
			if isinstance(elec, Number):
				elecKey = 'elec' + str(elec)
			elif elec == 'avg':
				elecKey = 'avg'
			avgValues['avgLFP'][q].append(lfpPopData[pop][elecKey]['avg'])
		q+=1
	dfAvg = pd.DataFrame(avgValues['avgLFP'], index=pops)


	if verbose:
		return dfPeak, dfAvg, peakValues, avgValues, lfpPopData 
	else:
		return dfPeak, dfAvg
def plotDataFrames(dataFrame, electrodes=None, pops=None, title=None, cbarLabel=None, figSize=None, savePath=None, saveFig=True):
	#### --> This function will plot a heatmap of the peak or average LFP amplitudes across electrodes & cell populations
	### dataFrame: pandas dataFrame  --> These can be obtained from getDataFrames function above)
	### electrodes: list 	--> DEFAULT: use all electrodes + 'avg'
	### pops: list 			--> DEFAULT: all cortical pops 
	### title: str  		--> Optional; title of the entire figure
	### cbarLabel: str 		--> DEFAULT: 'LFP amplitudes (mV)'  -->  (label on the color scale bar)
	### figSize: tuple 		--> DEFAULT: (12,6)
	### savePath: str, path to directory where figures should be saved  --> DEFAULT: '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/'
	### saveFig: bool 		-->  DEFAULT: True 


	## Set label for color scalebar 
	if cbarLabel is None:
		cbarLabel = 'LFP amplitude (mV)'
	elif cbarLabel == 'CSD':
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
		# print('electrodeLabels' + str(electrodeLabels)) 	## <-- TESTING LINE (to make sure electrode labels for the plot are coming out as intended) 


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
	titleFontSize = 20
	labelFontSize = 15
	tickFontSize = 10


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
	ax = sns.heatmap(pivotedDataFrame, xticklabels=x_axis_labels, yticklabels=y_axis_labels, linewidth=0.4, cbar_kws={'label': cbarLabel})

	## Set labels on x and y axes 
	plt.xlabel('Cell populations', fontsize=labelFontSize)
	plt.xticks(rotation=45, fontsize=tickFontSize)
	plt.ylabel('Channel', fontsize=labelFontSize)
	plt.yticks(rotation=0, fontsize=tickFontSize) 

	if saveFig:
		if savePath is None:
			prePath = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/LFP_heatmaps/'
		else:
			prePath = savePath
		fileName = 'heatmap.png'
		pathToFile = prePath + fileName
		plt.savefig(pathToFile, dpi=300)

	plt.show()

	return ax

## Pops with the highest contributions to LFP / CSD signal at a particular electrode ## 
def evalPops(dataFrame, electrode, verbose=0):   #  TODO --> elecOnly=1
	##### -->  This function outputs a dict with the pops & values from the main & adjacent electrodes specified in the 'electrode' arg 
	## dataFrame: pandas dataFrame --> with peak or avg LFP or CSD data of each pop at each electrode --> output of getDataFrames
	## electrode: int 	--> electrode where the oscillation event of interest occurred - focus on data from this electrode plus the ones immediately above and below it
	## verbose: bool 	--> Default: 0 (False), only returns maxPopsValues. True returns dict + main electrode data frames. 

	#### SUBSET THE DATAFRAME ACCORDING TO SPECIFIED ELECTRODE & ADJACENT ELECTRODE(S) ####  
	if electrode == 0:									# if electrode is 0, then this is topmost, so only include 0, 1  	# dataFrameSubset = dataFrame[[elec, bottomElec]]
		# Determine electrodes to use for subsetting 
		elec = electrode
		bottomElec = electrode + 1
		topElec = None
		# Subset the dataFrame & create another with absolute values of the subset 
		dataFrameSubsetBottomElec = dataFrame[[bottomElec]]
		dataFrameSubsetBottomElecAbs = dataFrameSubsetBottomElec.abs()

	elif electrode == 19:	# if electrode is 19, this is bottom-most, so only include 18, 19 	# dataFrameSubset = dataFrame[[topElec, elec]]
		# Determine electrodes to use for subsetting 
		elec = electrode
		bottomElec = None
		topElec = electrode - 1
		# Subset the dataFrame & create another with absolute values of the subset 
		dataFrameSubsetTopElec = dataFrame[[topElec]]
		dataFrameSubsetTopElecAbs = dataFrameSubsetTopElec.abs()

	elif electrode > 0 and electrode < 19:														# dataFrameSubset = dataFrame[[topElec, elec, bottomElec]]
		# Determine electrodes to use for subsetting 
		elec = electrode
		bottomElec = electrode + 1
		topElec = electrode - 1
		# Subset the dataFrame & create another with absolute values of the subset 
		dataFrameSubsetTopElec = dataFrame[[topElec]]
		dataFrameSubsetTopElecAbs = dataFrameSubsetTopElec.abs()
		dataFrameSubsetBottomElec = dataFrame[[bottomElec]]
		dataFrameSubsetBottomElecAbs = dataFrameSubsetBottomElec.abs()

	dataFrameSubsetElec = dataFrame[[elec]]
	dataFrameSubsetElecAbs = dataFrameSubsetElec.abs()

	# So now I have: 
	## dataframes with electrode subset and adjacent electrode(s) subset (e.g. dataFrameSubsetElec, dataFrameSubsetBottomElec, dataFrameSubsetTopElec)
	## dataframes with absolute values of electrode subset and adjacent electrode(s) subset (e.g. dataFrameSubsetElecAbs, dataFrameSubsetBottomElecAbs)


	#### MAIN ELECTRODE #### 
	## Make dict for storage 
	maxPopsValues = {}
	maxPopsValues['elec'] = {}
	maxPopsValues['elec']['electrodeNumber'] = electrode
	## Values / pops associated with main electrode 
	dfElecAbs = dataFrameSubsetElecAbs.sort_values(by=[elec], ascending=False) # elec OR maxPopsValues['elec']['electrodeNumber']
		#### SUBSET VIA HARDCODED THRESHOLD --> 
	# dfElecSub = dfElecAbs[dfElecAbs[elec]>0.1]
		#### OR SUBSET BY TAKING ANY VALUE > (MAX VALUE / 3) ---> 
	# maxValueElec = list(dfElecAbs.max())[0]
	# dfElecSub = dfElecAbs[dfElecAbs[elec]>(maxValueElec/3)]
		#### OR SUBSET BY TAKING THE TOP 5 (OR TOP 3) CONTRIBUTORS ----> 
	dfElecSub = dfElecAbs.head(3) # dfElecAbs.head(5)
	elecPops = list(dfElecSub.index)
	## Find the non-absolute values of the above, and store this info along with the associated populations, into the storage dict:
	for pop in elecPops:
		maxPopsValues['elec'][pop] = list(dataFrameSubsetElec.loc[pop])[0]

	#### BOTTOM ADJACENT ELECTRODE #### 
	if bottomElec is not None:
		## Dict for storage 
		maxPopsValues['bottomElec'] = {}
		maxPopsValues['bottomElec']['electrodeNumber'] = bottomElec
		## Values / pops associated with bottom adjacent electrode 
		dfBottomElecAbs = dataFrameSubsetBottomElecAbs.sort_values(by=[bottomElec], ascending=False) # bottomElec OR maxPopsValues['bottomElec']['electrodeNumber']
		dfBottomElecSub = dfBottomElecAbs.head(3) # dfBottomElecAbs.head(5)		# dfBottomElecSub = dfBottomElecAbs[dfBottomElecAbs[bottomElec]>0.1]
		bottomElecPops = list(dfBottomElecSub.index)
		## Find the non-absolute values of the above, and store this info along with the associated populations, into the storage dict:
		for pop in bottomElecPops:
			maxPopsValues['bottomElec'][pop] = list(dataFrameSubsetBottomElec.loc[pop])[0]

	#### TOP ADJACENT ELECTRODE #### 
	if topElec is not None:
		## Dict for storage 
		maxPopsValues['topElec'] = {}
		maxPopsValues['topElec']['electrodeNumber'] = topElec
		## Values / pops associated with top adjacent electrode 
		dfTopElecAbs = dataFrameSubsetTopElecAbs.sort_values(by=[topElec], ascending=False) # topElec OR maxPopsValues['topElec']['electrodeNumber']
		dfTopElecSub = dfTopElecAbs.head(3) # dfTopElecAbs.head(5) 		# dfTopElecSub = dfTopElecAbs[dfTopElecAbs[topElec]>0.1]
		topElecPops = list(dfTopElecSub.index)
		## Find the non-absolute values of the above, and store this info along with the associated populations, into the storage dict:
		for pop in topElecPops:
			maxPopsValues['topElec'][pop] = list(dataFrameSubsetTopElec.loc[pop])[0]

	if verbose:
		return maxPopsValues, dfElecSub, dataFrameSubsetElec
	else:
		return maxPopsValues

## Spike Activity: data and plotting ## 
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
		# Print out the values for view in terminal 
		print('left: ' + str(left))
		print('right: ' + str(right))
		print('minT: ' + str(minT))
		print('maxT: ' + str(maxT))
		print('alignoffset: ' + str(alignoffset))
		print('w2: ' + str(w2))
		# Calculate idx0 and idx1 for before, and beforeT
		idx0_before = max(0,left - w2)
		idx1_before = left 
		beforeT = (maxT-minT) * (idx1_before - idx0_before) / (right - left + 1)
		print('beforeT: ' + str(beforeT))
		# Calculate idx0 and idx1 for after, and afterT
		idx0_after = int(right)
		idx1_after = min(idx0_after + w2, 230000)	# max(csdData.shape[0],csdData.shape[1]))  # <-- believe this would be the number of time points (~230000?)
		afterT = (maxT-minT) * (idx1_after - idx0_after) / (right - left + 1)		#		
		print('afterT: ' + str(afterT))
		# Update outputData with any values necessary for the plotting function
		outputData.update({'alignoffset': alignoffset})
	else:
		print('No oscillation event data detected!')

	## Calculate time ranges for data gathering 		-->  ### NOTE: the time ranges WITHOUT ALIGNOFFSET ADDED TO EITHER ELEMENT OF THE LISTS ensures that spike counts remain accurate!! 
	# during osc event
	T_during_withoutAlignOffset = [minT, maxT]  # + alignoffset to both elements 
	T_during = [minT+alignoffset, maxT+alignoffset]  
	# during osc event + before & after 
	T_full_withoutAlignOffset = [(minT-beforeT), (maxT+afterT)]	 # + alignoffset to both elements 
	T_full = [(minT-beforeT) + alignoffset, (maxT+afterT) + alignoffset]
	# print out time ranges to check 
	print('T_during: ' + str(T_during))
	print('T_full: ' + str(T_full))
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
		print('vminDuring: ' + str(vminDuring))
		vmaxDuring = np.array(specDuring).max()
		print('vmaxDuring: ' + str(vmaxDuring))
		vcDuring = [vminDuring, vmaxDuring]
		outputData.update({'vcDuring': vcDuring})

		vminFull = np.array(specFull).min()
		print('vminFull: ' + str(vminFull))
		vmaxFull = np.array(specFull).max()
		print('vmaxFull: ' + str(vmaxFull))
		vcFull = [vminFull, vmaxFull]
		outputData.update({'vcFull': vcFull})

	# save figure data
	# figData = {'histData': histData, 'histT': histoT, 'include': include, 'timeRange': timeRange, 'binSize': binSize}
	outputData.update({'allSignalDuring': allSignalDuring, 'allFreqsDuring': allFreqsDuring, 
					'allSignalFull': allSignalFull, 'allFreqsFull': allFreqsFull})


	return outputData # {'allSignal': allSignal, 'allFreqs':allFreqs}
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
		# Print out the values for view in terminal 
		print('left: ' + str(left))
		print('right: ' + str(right))
		print('minT: ' + str(minT))
		print('maxT: ' + str(maxT))
		print('alignoffset: ' + str(alignoffset))
		print('w2: ' + str(w2))
		# Calculate idx0 and idx1 for before, and beforeT
		idx0_before = max(0,left - w2)
		idx1_before = left 
		beforeT = (maxT-minT) * (idx1_before - idx0_before) / (right - left + 1)
		print('beforeT: ' + str(beforeT))
		# Calculate idx0 and idx1 for after, and afterT
		idx0_after = int(right)
		idx1_after = min(idx0_after + w2, 230000)	# max(csdData.shape[0],csdData.shape[1]))  # <-- believe this would be the number of time points (~230000?)
		afterT = (maxT-minT) * (idx1_after - idx0_after) / (right - left + 1)		#		
		print('afterT: ' + str(afterT))
		# Update outputData with any values necessary for the plotting function
		outputData.update({'alignoffset': alignoffset})
	else:
		print('No oscillation event data detected!')


	# Calculate time range to gather data for (during oscillation event only, or during time of oscillation event + buffer before and after)
	### NOTE: the time ranges WITHOUT ALIGNOFFSET ADDED TO EITHER ELEMENT OF THE LISTS ensures that spike counts remain accurate!! 
	# before osc event
	T_before_withoutAlignOffset = [(minT-beforeT), minT] # + alignoffset to both elements 
	T_before = [(minT-beforeT) + alignoffset, minT + alignoffset]
	# during osc event
	T_during_withoutAlignOffset = [minT, maxT]  # + alignoffset to both elements 
	T_during = [minT+alignoffset, maxT+alignoffset]  
	# after osc event 
	T_after_withoutAlignOffset = [maxT, (maxT+afterT)] # + alignoffset to both elements
	T_after = [maxT+alignoffset, (maxT+afterT)+alignoffset]
	# print out time ranges to check 
	print('T_before: ' + str(T_before))
	print('T_during: ' + str(T_during))
	print('T_after: ' + str(T_after))
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
	# figData = {'histoData': histoData, 'histoT': histoT, 'include': include, 'binSize': binSize}	# 'timeRange': timeRange,  ## RENAME figData --> outputData
	outputData.update({'histoDataBefore': histoDataBefore, 'histoDataDuring': histoDataDuring, 'histoDataAfter': histoDataAfter, 
				'histoTBefore': histoTBefore, 'histoTDuring': histoTDuring, 'histoTAfter': histoTAfter, 
				'include': include, 'binSize': binSize})

	return outputData 				# {'histoData': histoData, 'histoT': histoT, 'include': include}	# 'timeRange': timeRange, 
#  def getSpikeData outputs spike data dicts for use in plotCombinedSpike
def getSpikeData(dataFile, pop, graphType, oscEventInfo): 
	### dataFile: path to .pkl data file to load 
	### pop: list or str --> which pop to include 
	### graphType: str --> either 'hist' or 'spect'
	### oscEventInfo: dict 
				## DEPRECATED --> ### timeRange: list --> e.g. [start, stop]

	# Load data file
	sim.load(dataFile, instantiate=False)

	# Pops
	if type(pop) is str:
		popList = [pop]
	elif type(pop) is list:
		popList = pop

	# Set up which kind of data -- i.e. spectrogram or histogram 
	if graphType == 'spect':
		spikeDict = getRateSpectrogramData(include=popList, oscEventInfo=oscEventInfo)
	elif graphType == 'hist':
		spikeDict = getSpikeHistData(include=popList, oscEventInfo=oscEventInfo, binSize=5, graphType='bar', measure='rate') ## sim.analysis.getSpikeHistData

	return spikeDict 
def plotCombinedSpike(pop, colorDict, plotTypes=['spectrogram', 'histogram'], hasBefore=1, hasAfter=1, spectDict=None, histDict=None, figSize=(10,7), colorMap='jet', minFreq=None, maxFreq=None, vmaxContrast=None, savePath=None, saveFig=True):
	# DEPRECATED --> ### timeRange: list 				--> e.g. [start, stop]
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
	labelFontSize = 12
	titleFontSize = 20

	#### SPECTROGRAM CALCULATIONS ####-------------------------------------------------------
	if 'spectrogram' in plotTypes:
		# allSignal = spectDict['allSignal']
		# allFreqs = spectDict['allFreqs']

		if hasBefore and hasAfter:
			allSignal = spectDict['allSignalFull'] ### --> [0] ?? --> okay yeah that happens later! 
			allFreqs = spectDict['allFreqsFull']
			timeRange = spectDict['T_full']
			print('timeRange: ' + str(timeRange))
			vc = spectDict['vcFull']
		else:
			allSignal = spectDict['allSignalDuring']
			allFreqs = spectDict['allFreqsDuring']
			timeRange = spectDict['T_during']
			print('timeRange: ' + str(timeRange))
			vc = spectDict['vcDuring']



		# Set color contrast parameters (vmin / vmax)
		orig_vmin = vc[0]
		orig_vmax = vc[1]

		if vmaxContrast is None:
			vmin = orig_vmin 	# np.amin(imshowSignal)		# None
			vmax = orig_vmax	# np.amax(imshowSignal)		# None
		else:
			vmin = orig_vmin 	# np.amin(imshowSignal)
			vmax = orig_vmax / vmaxContrast # np.amax(imshowSignal) / vmaxContrast
			vc = [vmin, vmax]
		# ## TRIAL VMIN VMAX LINES
		# vmin = 1.0e-07
		# vmax = 0.2
		# vc = [vmin, vmax]

		# Set frequencies to be plotted (minFreq / maxFreq)
		if minFreq is None:
			minFreq = np.amin(allFreqs[0])			# DEFAULT: 1
		if maxFreq is None:
			maxFreq = np.amax(allFreqs[0])			# DEFAULT: 100
		elif maxFreq is not None:					# maxFreq must be an integer!! 
			if type(maxFreq) is not int:
				maxFreq = round(maxFreq)

		# Set up imshowSignal
		imshowSignal = allSignal[0][minFreq:maxFreq+1]  ## Can't tell if the +1 is necessary or not here - leave it in for now!! 


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
				interpolation='None', aspect='auto', cmap=plt.get_cmap(colorMap), vmin=vc[0], vmax=vc[1])	# vmin=vmin, vmax=vmax)
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad = 0.2)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar to be in scientific notation
		fmt.set_powerlimits((0,0))
		plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		ax1.set_title('Spike Rate Spectrogram for ' + popToPlot, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq) ## Redundant if ax1.imshow extent arg uses minFreq and maxFreq (and if imshowSignal is sliced as well)

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
		if hasBefore and hasAfter:
			T_before = histDict['T_before']
			T_after = histDict['T_after']
			ax2.set_xlim(left=T_before[0], right=T_after[1]) 		# # ax2.set_xlim(left=timeRange[0], right=timeRange[1])
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
				interpolation='None', aspect='auto', cmap=plt.get_cmap(colorMap), vmin=vc[0], vmax=vc[0]) # vmin=vmin, vmax=vmax) ## CHANGE maxFreq // np.amax(allFreqs[0])
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad = 0.2)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar to be in scientific notation
		fmt.set_powerlimits((0,0))
		plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		ax1.set_title('Spike Rate Spectrogram for ' + popToPlot, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq)

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

	## Save figure
	if saveFig:
		if savePath is None:
			prePath = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/popContribFigs/' 	# popContribFigs_cmapJet/'
		else:
			prePath = savePath
		# fileName = pop + '_combinedSpike.png'
		pathToFile = prePath + figFilename	# fileName
		plt.savefig(pathToFile, dpi=600)

	## Show figure
	plt.show()

## LFP: data and plotting ## 
def getLFPDataDict(dataFile, pop, plotType, timeRange, electrode):
	### dataFile: str 			--> path to .pkl data file to load for analysis 
	### pop: str or list  		--> cell population to get the LFP data for 
	### plotType: str or list 	--> 'spectrogram' or 'timeSeries'
	### timeRange: list  		-->  e.g. [start, stop]
	### electrode: list or int or str designating the electrode of interest --> e.g. 10, [10], 'avg'
		### NOT IN USE CURRENTLY: filtFreq: list --> DEFAULT: None; frequencies to bandpass filter lfp data  <-- supposed to come from def getBandpassRange() 

	# Load data file 
	sim.load(dataFile, instantiate=False)

	# Pops
	if type(pop) is str:
		popList = [pop]
	elif type(pop) is list:
		popList = pop

	# Electrodes
	if type(electrode) is not list: 
		electrodeList = [electrode]
	else:
		electrodeList = electrode


	# Set up which kind of data -- i.e. timeSeries or spectrogram 
	if type(plotType) is str:
		plots = [plotType]
	elif type(plotType) is list:
		plots = plotType

	lfpOutput = getLFPData(pop=popList, timeRange=timeRange, electrodes=electrodeList, plots=plots) # sim.analysis.getLFPData # filtFreq=filtFreq (see above; in args)

	return lfpOutput
def plotCombinedLFP(timeRange, pop, colorDict, plotTypes=['spectrogram', 'timeSeries'], spectDict=None, timeSeriesDict=None, figSize=(10,7), colorMap='jet', minFreq=None, maxFreq=None, vmaxContrast=None, electrode=None, savePath=None, saveFig=True): # electrode='avg',
	### timeRange: list 						--> [start, stop]
	### pop: list or str 						--> relevant population to plot data for 
	### colorDict: dict 						--> corresponds pop to color 
	### plotTypes: list 						--> DEFAULT: ['spectrogram', 'timeSeries'] 
	### spectDict: dict 						--> Contains spectrogram data; output of getLFPDataDict
	### timeSeriesDict: dict 					--> Contains timeSeries data; output of getLFPDataDict
	### figSize: tuple 							--> DEFAULT: (10,7)
	### colorMap: str 							--> DEFAULT: 'jet' 	--> cmap for ax.imshow lines --> Options are currently 'jet' or 'viridis' 
	### maxFreq: int 							--> whole number that determines the maximum frequency plotted on the spectrogram 
	### minFreq: int 							--> whole number that determines the minimum frequency plotted on the spectrogram
	### vmaxContrast: float or int 				--> Denominator This will help with color contrast if desired!!!, e.g. 1.5 or 3
	### electrode: int 					--> electrode at which to plot the CSD data 
		# ^^^ REPLACEMENT ### titleElectrode: str or (1-element) list	-->  FOR USE IN PLOT TITLES !! --> This is for the electrode that will appear in the title 
	### savePath: str 	  						--> Path to directory where fig should be saved; DEFAULT: '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/popContribFigs/'
	### saveFig: bool 							--> DEFAULT: True 
		### NOTE: RE-WORK HOW THE ELECTRODE GETS FACTORED INTO THE TITLE AND SAVING!! 

	# Get relevant pop
	if type(pop) is str:
		popToPlot = pop
	elif type(pop) is list:
		popToPlot = pop[0]

	## Create figure 
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figSize)

	# Set electrode variable, for plot title(s)
	if type(electrode) is list:
		electrode = electrode[0]


	## Set font size 
	labelFontSize = 15 		# NOTE: in plotCombinedSpike, labelFontSize = 12
	titleFontSize = 20


	#### SPECTROGRAM CALCULATIONS ####-------------------------------------------------------
	if 'spectrogram' in plotTypes:
		spec = spectDict['spec']

		S = spec[0].TFR
		F = spec[0].f
		T = timeRange

		## Set up vmin / vmax color contrasts 
		orig_vmin = np.array([s.TFR for s in spec]).min()
		orig_vmax = np.array([s.TFR for s in spec]).max()

		if vmaxContrast is None:
			vmin = orig_vmin
			vmax = orig_vmax
		else:
			vmin = orig_vmin
			vmax = orig_vmax / vmaxContrast 
			print('original vmax: ' + str(orig_vmax))
			print('new vmax: ' + str(vmax))

		vc = [vmin, vmax]


		## Set up minFreq and maxFreq for spectrogram
		if minFreq is None:
			minFreq = np.amin(F) 	# 1
		if maxFreq is None:
			maxFreq = np.amax(F)	# 100
		print('minFreq: ' + str(minFreq) + ' Hz')
		print('maxFreq: ' + str(maxFreq) + ' Hz')

		## Set up imshowSignal
		imshowSignal = S[minFreq:maxFreq+1] ## NOTE: Is the +1 necessary here or not? Same question for spiking data. Leave it in for now. 



	#### TIME SERIES CALCULATIONS ####-------------------------------------------------------
	if 'timeSeries' in plotTypes:
		t = timeSeriesDict['t']
		lfpPlot = timeSeriesDict['lfpPlot']



	#### PLOTTING ####-------------------------------------------------------
	# Plot both 
	if 'spectrogram' in plotTypes and 'timeSeries' in plotTypes:

		### PLOT SPECTROGRAM ###
		# spectrogram title #
		spectTitle = 'LFP Spectrogram for ' + popToPlot + ', electrode ' + str(electrode)   ## NOTE: make condition for avg elec? 
		# plot and format # 
		ax1 = plt.subplot(2, 1, 1)
		img = ax1.imshow(imshowSignal, extent=(np.amin(T), np.amax(T), minFreq, maxFreq), origin='lower', interpolation='None', aspect='auto', 
			vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap(colorMap))	 # imshowSignal --> S 	 	# minFreq --> np.amin(F)  	# maxFreq --> np.amax(F) 
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad=0.2)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar scientific notation
		fmt.set_powerlimits((0,0))
		plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		ax1.set_title(spectTitle, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq)		## Redundant if ax1.imshow extent reflects minFreq and maxFreq


		### PLOT TIME SERIES ###
		# timeSeries title # 
		timeSeriesTitle = 'LFP Signal for ' + popToPlot + ', electrode ' + str(electrode)
		# plot and format #
		lw = 1.0
		ax2 = plt.subplot(2, 1, 2)
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes('right', size='3%', pad=0.2)
		cax2.axis('off')
		ax2.plot(t[0:len(lfpPlot)], lfpPlot, color=colorDict[popToPlot], linewidth=lw)
		ax2.set_title(timeSeriesTitle, fontsize=titleFontSize)
		ax2.set_xlabel('Time (ms)', fontsize=labelFontSize)
		ax2.set_xlim(left=timeRange[0], right=timeRange[1])
		ax2.set_ylabel('LFP Amplitudes (mV)', fontsize=labelFontSize)

		# For potential figure saving
		figFilename = popToPlot + '_combinedLFP_elec_' + str(electrode) + '.png'

	# Plot only spectrogram 
	elif 'spectrogram' in plotTypes and 'timeSeries' not in plotTypes:
		### PLOT SPECTROGRAM ###
		# spectrogram title # 
		spectTitle = 'LFP Spectrogram for ' + popToPlot + ', electrode ' + str(electrode)
		# plot and format # 
		ax1 = plt.subplot(1, 1, 1)
		img = ax1.imshow(imshowSignal, extent=(np.amin(T), np.amax(T), minFreq, maxFreq), origin='lower', interpolation='None', aspect='auto', 
			vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap(colorMap))  # imshowSignal --> S 	 	# minFreq --> np.amin(F)  	# maxFreq --> np.amax(F) 
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad=0.2)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar scientific notation
		fmt.set_powerlimits((0,0))
		plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		ax1.set_title(spectTitle, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq)		## Redundant if ax1.imshow extent reflects minFreq and maxFreq

		# For potential figure saving
		figFilename = popToPlot + '_LFP_spectrogram_elec_' + str(electrode) + '.png'

	# Plot only timeSeries 
	elif 'spectrogram' not in plotTypes and 'timeSeries' in plotTypes:
		### PLOT TIME SERIES ###
		# timeSeries title # 
		timeSeriesTitle = 'LFP Signal for ' + popToPlot + ', electrode ' + str(electrode)
		# plot and format # 
		lw = 1.0
		ax2 = plt.subplot(1, 1, 1)
		divider2 = make_axes_locatable(ax2)
		cax2 = divider2.append_axes('right', size='3%', pad=0.2)
		cax2.axis('off')
		ax2.plot(t[0:len(lfpPlot)], lfpPlot, color=colorDict[popToPlot], linewidth=lw)
		ax2.set_title(timeSeriesTitle, fontsize=titleFontSize)
		ax2.set_xlabel('Time (ms)', fontsize=labelFontSize)
		ax2.set_xlim(left=timeRange[0], right=timeRange[1])
		ax2.set_ylabel('LFP Amplitudes (mV)', fontsize=labelFontSize)

		# For potential figure saving
		figFilename = popToPlot + '_LFP_timeSeries_elec_' + str(electrode) + '.png'

	plt.tight_layout()

	## Save figure 
	if saveFig:
		if savePath is None:
			prePath = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/popContribFigs/' 		# popContribFigs_cmapJet/'
		else:
			prePath = savePath
		pathToFile = prePath + figFilename
		plt.savefig(pathToFile, dpi=300)

	## Show figure 
	plt.show()

def getSumLFP(dataFile, pops, elecs=False, timeRange=None, showFig=False):
	# THIS FUNCTON ALLOWS YOU TO ADD TOGETHER LFP CONTRIBUTIONS FROM ARBITRARY POPS AT SPECIFIED ELECTRODES
	### dataFile: str --> .pkl file to load w/ simulation data 
	### pops: list OR dict --> list if just a list of pops, dict if including electrodes (e.g. {'IT3': 1, 'IT5A': 10, 'PT5B': 11})
			### popElecDict: dict --> e.g. {'IT3': 1, 'IT5A': 10, 'PT5B': 11}
	### timeRange: list --> e.g. [start, stop]
	### showFig: bool --> Determines whether or not plt.show() will be called 
	### elecs: bool ---> False by default; means that electrodes will NOT be taken into account and output will be the total LFP.

	print('Calculating combined LFP signal')

	sim.load(dataFile, instantiate=False)  ### Figure out efficiency here!!  

	if timeRange is None:   	# if timeRange is None, use the whole simulation! 
		timeRange = [0, sim.cfg.duration]


	lfpData = {} 	# dictionary to store lfp data 

	## Calculating sum of LFP signals at specified electrodes! 
	if elecs:
		if type(pops) is list:  	# We need pops w/ corresponding electrodes if we want to calculate summed LFP at specfic electrodes! 
			print('pops argument needs to be a dict with pop as key & corresponding electrode as value, e.g. {\'IT3\': 1, \'IT5A\': 10, \'PT5B\': 11}')
			lfpData = None
		elif type(pops) is dict:
			for pop in pops:
				lfpData[pop] = {}
				elec = pops[pop]		# electrode that corresponds with the current pop being evaluated! 
				popLFPdata = np.array(sim.allSimData['LFPPops'][pop])[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]
				lfpData[pop]['elec'] = popLFPdata[:, elec]

			popList = list(pops.keys()) 	## Make a list from the 'pops' dict keys -- this will be a list of the relevant populations 
			lfpData['sum'] = np.zeros(lfpData[popList[0]]['elec'].shape) 	## have to establish an array of zeros first so the below 'for' loop functions properly
			for i in range(len(pops)):
				lfpData['sum'] += lfpData[popList[i]]['elec']
	## Calculating sum of LFP signals over ALL electrodes! 
	else:
		print('Calculating summed LFP signal from ' + str(pops) + ' over all electrodes')
		for pop in pops:
			lfpData[pop] = np.array(sim.allSimData['LFPPops'][pop])[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]
		lfpData['sum'] = np.zeros(lfpData[pops[0]].shape)  ## have to establish an array of zeros first so the below 'for' loop functions properly
		for i in range(len(pops)):
			lfpData['sum'] += lfpData[pops[i]]

	### PLOTTING ### 
	if showFig:
		t = np.arange(timeRange[0], timeRange[1], sim.cfg.recordStep)
		plt.figure(figsize = (12,7))
		plt.xlabel('Time (ms)')
		### Create title of plot ### 
		titlePreamble = 'Summed LFP signal from'
		popsInTitle = ''
		if elecs and type(pops) is dict:
			for pop in pops:
				popsInTitle += ' ' + pop + '(elec ' + str(pops[pop]) + ')'
			plt.ylabel('LFP Amplitude (mV)')	### y axis label 
			plt.plot(t, lfpData['sum'])   #### PLOT HERE 
			lfpSumTitle = titlePreamble + popsInTitle 
			plt.title(lfpSumTitle)
			plt.legend()
			plt.show()
		elif not elecs:
			print('Use plotLFP')
			# for pop in pops:
			# 	popsInTitle += ' ' + pop + ' '
			# for i in range(lfpData['sum'].shape[1]):  ## number of electrodes 
			# 	lfpPlot = lfpData['sum'][:, i]
			# 	lw = 1.0 ## line-width
			# 	ydisp = np.absolute(lfpPlot).max() 
			# 	plt.plot(t, -lfpPlot+(i*ydisp), linewidth=lw, label='electrode ' + str(i))

	return lfpData 

## CSD: data and plotting ## 
def getCSDdata(dataFile=None, outputType=['timeSeries', 'spectrogram'], oscEventInfo=None, dt=None, sampr=None, pop=None, norm=None, spacing_um=100, minFreq=1, maxFreq=110, stepFreq=0.25):
	#### Outputs a dict with CSD and other relevant data for plotting! 
	## dataFile: str     					--> .pkl file with recorded simulation 
	## outputType: list of strings 			--> options are 'timeSeries' +/- 'spectrogram' --> OR could be empty, if want csdData from all electrodes!! 
	## oscEventInfo: dict 					---> dict w/ chan, left, right, minT, maxT, alignoffset
		## RIGHT NOW THIS IS *NECESSARY* TO RUN THIS!!! 
		# DEPRECATED ## timeRange: list 						--> e.g. [start, stop]
		# DEPRECATED ## channel: list or None 				--> e.g. [8], None
	## dt: time step of the simulation 		--> (usually --> sim.cfg.recordStep)
	## sampr: sampling rate (Hz) 			--> (usually --> 1/(dt/1000))
	## pop: str or list 					--> e.g. 'ITS4' or ['ITS4']
	## norm: bool 							--> e.g. norm=True, norm=False; DEFAULT: True 
	## spacing_um: int 						--> 100 by DEFAULT (spacing between electrodes in MICRONS)
	## minFreq: float / int 				--> DEFAULT: 1 Hz  
	## maxFreq: float / int 				--> DEFAULT: 110 Hz 
	## stepFreq: float / int 				--> DEFAULT: 0.25 Hz 
			## TO DO: 
			###  --> Should I also have an lfp_input option so that I can get the CSD data of summed LFP data...?
			###  --> Should I make it such that output can include both timeSeries and spectrogram so don't have to run this twice? test this!! 

	## load .pkl simulation file 
	if dataFile:
		sim.load(dataFile, instantiate=False)
	else:
		print('No dataFile; will use data from dataFile already loaded elsewhere!')



	## Determine timestep, sampling rate, and electrode spacing 
	dt = sim.cfg.recordStep  	# or should I divide by 1000.0 up here, and then just do 1.0/dt below for sampr? 
	sampr = 1.0/(dt/1000.0) 	# divide by 1000.0 to turn denominator from units of ms to s
	#dt = sim.cfg.recordStep/1000.0	# This is in Seconds 
	#sampr = 1.0 / dt 				# Hz 
	spacing_um = spacing_um			# 100 um by default # 


	## Get LFP data   # ----> NOTE: SHOULD I MAKE AN LFP INPUT OPTION?????? FOR SUMMED LFP DATA??? (also noted above in arg descriptions)
	if pop is None:
		lfpData = sim.allSimData['LFP']
	else:
		if type(pop) is list:
			pop = pop[0]
		lfpData = sim.allSimData['LFPPops'][pop]

	## Set up dictionary for output data 
	outputData = {}


	#### ALL CSD DATA -- ALL CHANS, ALL TIMEPOINTS!!! 
	if norm is None:
		norm=True
	#csdData = csd.getCSD(LFP_input_data=lfpData, dt=dt, sampr=sampr, spacing_um=spacing_um, norm=norm, vaknin=True)
	csdData = getCSDa1dat(lfps=lfpData,sampr=sampr,spacing_um=spacing_um,minf=0.05,maxf=300,norm=norm,vaknin=True)
	tt = np.linspace(0, sim.cfg.duration, len(csdData[1])) 
	# Input full CSD data (and time array) into outputData dict 
	outputData.update({'csdData': csdData})
	outputData.update({'tt': tt}) ## Keep in mind this is in milliseconds! 


	## Extract oscillation event info 
	if oscEventInfo is not None: ### RIGHT NOW THIS ARG IS NECESSARY FOR RUNNING THIS FUNCTION!!! FIX THIS!!!!
		# Extract chan, left, right, minT, maxT, alignoffset, w2   # RECALL: HAVE TO CORRECT FOR 6_11 !!! THIS WORKS FOR 0_6 AS IS!!! 
		chan = oscEventInfo['chan']
		left = oscEventInfo['left']
		right = oscEventInfo['right']
		minT = oscEventInfo['minT']
		maxT = oscEventInfo['maxT']
		alignoffset = oscEventInfo['alignoffset']
		w2 = oscEventInfo['w2']
		# # Print out the values for view in terminal 
		# print('channel: ' + str(chan))
		# print('left: ' + str(left))
		# print('right: ' + str(right))
		# print('minT: ' + str(minT))
		# print('maxT: ' + str(maxT))
		# print('alignoffset: ' + str(alignoffset))
		# print('w2: ' + str(w2))
		# # Calculate idx0 and idx1 for before, and beforeT
		idx0_before = max(0,left - w2)
		idx1_before = left 
		beforeT = (maxT-minT) * (idx1_before - idx0_before) / (right - left + 1)
		outputData.update({'beforeT': beforeT})
		# Calculate idx0 and idx1 for after, and afterT
		idx0_after = int(right)
		idx1_after = min(idx0_after + w2,max(csdData.shape[0],csdData.shape[1]))
		afterT = (maxT-minT) * (idx1_after - idx0_after) / (right - left + 1)		#		print('afterT: ' + str(afterT))
		outputData.update({'afterT': afterT})
	else:
		chan=None
		print('No oscillation event data detected!')




	#### ALL CHANS, OSC EVENT TIMEPOINTS!! 
	csdDataAllChans = csdData[:,left:right]	 	# NOTE: THIS ONLY WORKS IF left and right EXIST!!! AKA IF oscEventInfo is not None!!!				# idx0_before:idx1_after]  
	# Input full CSD data into outputData dict 
	outputData.update({'csdDataAllChans': csdDataAllChans})  

	#### ALL CHANS, OSC EVENT + BUFFER!!
	csdDataAllChans_plusTimeBuffer = csdData[:,idx0_before:idx1_after]
	# Input this data into outputData dict 
	outputData.update({'csdDataAllChans_plusTimeBuffer': csdDataAllChans_plusTimeBuffer})  # REMEMBER: ALL CHANS, OSC EVENT TIMEPOINTS!! 


	# timeSeries --------------------------------------------
	if 'timeSeries' in outputType:  ### make case for when it IS None  # and oscEventInfo is not None
		print('Returning timeSeries data')

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
		tt_during = np.linspace(minT,maxT,len(csdDuring)) + alignoffset				# 		duringT = np.linspace(minT,maxT,len(csdDuring)) + alignoffset 
		# (3) Input time and CSD data for during osc event into outputData dict 
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



		### NOTE: ^^ Replaced "full" with OscChan_plusTimeBuffer <--- CHECK FOR DEPENDENCIES AFFECTED ELSEWHERE!!! 

		# Update outputData with the channel info as well
		outputData.update({'chan': chan})

		# Get xlim for plotting
		xl = (minT-beforeT + alignoffset, maxT+afterT + alignoffset)
		outputData.update({'xl': xl})


	# spectrogram -------------------------------------------
	if 'spectrogram' in outputType:
		print('Returning spectrogram data')

		spec = []
		freqList = None


		## Spectrogram Data Calculations ## 
		fs = sampr # int(1000.0 / sim.cfg.recordStep)



		##############################
		#### DURING THE OSC EVENT #### 
		csdDuring = csdData[chan,left:right]
		specDuring = []
		specDuring.append(MorletSpec(csdDuring, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq, getphase=True,width=7.0)) # lfreq=freqList, # # Seems this is only used for the fft circumstance...? 
		## vmin, vmax 
		vminDuring = np.array([s.TFR for s in specDuring]).min()
		vmaxDuring = np.array([s.TFR for s in specDuring]).max()
		vcDuring = [vminDuring, vmaxDuring]
		## T 
		T_during = [minT + alignoffset, maxT + alignoffset]					# tt_during = np.linspace(minT,maxT,len(csdDuring)) + alignoffset
		## F, S 															# commented out the normspec lines 
		F_during = specDuring[0].f
		S_during = specDuring[0].TFR
		## outputData update 
		outputData.update({'T_during': T_during, 'F_during': F_during, 'S_during': S_during, 'vcDuring': vcDuring})
		outputData.update({'specDuring': specDuring})


		###############################################
		#### DURING THE OSC EVENT + BEFORE & AFTER ####   ### CHANGED THE TERMINOLOGY BUT MAYBE SHOULDN'T HAVE - SEE ABOVE!!! 
		csdFull = csdData[chan, idx0_before:idx1_after] 
		specFull = []
		specFull.append(MorletSpec(csdFull, fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq, getphase=True,width=7.0))  # lfreq=freqList, 
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
		f = freqList if freqList is not None else np.arange(minFreq, maxFreq+stepFreq, stepFreq)	# np.array(range(minFreq, maxFreq+1, stepFreq))   # only used as output for user # with arange stepFreq can be a non-integer value!! 
		outputData.update({'freqs': f[f<=maxFreq]})

		# Update outputData with the channel info as well
		outputData.update({'chan': chan})


	return outputData
def plotCombinedCSD(pop, electrode, colorDict, customPopTitle=None, timeSeriesDict=None, spectDict=None, alignoffset=None, vmaxContrast=None, colorMap='jet', figSize=(10,7), minFreq=None, maxFreq=None, plotTypes=['timeSeries', 'spectrogram'], hasBefore=1, hasAfter=1, savePath=None, saveFig=True):
	### pop: list or str 				--> relevant population to plot data for 
	### electrode: int 					--> electrode at which to plot the CSD data 
	### colorDict: dict 				--> corresponds pop to color 
	### customPopTitle: str 
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
		### TO DO: CHANGE ALL INSTANCES OF "ELECTRODE" TO "CHANNEL"

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
	labelFontSize = 20 # 15  ## NOTE: spike has this as 12, lfp plotting has this as 15 
	titleFontSize = 30
	tickFontSize = 17 #15

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


		# # These lines will work as long as the getCSDdata function that retrieves the spectDict had only 1 electrode in electrode arg!!
		# S = spectDict['S']
		# F = spectDict['F']
		# T = spectDict['T']				# timeRange

		# # vmin and vmax  -- adjust for color contrast purposes, if desired 
		# vc = spectDict['vc']


		orig_vmin = vc[0]
		orig_vmax = vc[1]

		if vmaxContrast is None:
			vmin = orig_vmin
			vmax = orig_vmax
		else:
			print('original vmax: ' + str(orig_vmax))
			vmin = orig_vmin
			vmax = orig_vmax / vmaxContrast
			print('new vmax: ' + str(vmax))
			vc = [vmin, vmax]


		## Set up minFreq and maxFreq for spectrogram
		if minFreq is None:
			minFreq = np.amin(F) 	# 1
		if maxFreq is None:
			maxFreq = np.amax(F)	# 100

		## Set up imshowSignal
		imshowSignal = S 			#[minFreq:maxFreq+1] ## NOTE: Is the +1 necessary here or not? Same question for spiking data. Leave it in for now. 


	#### TIME SERIES CALCULATIONS ####-------------------------------------------------------
	if 'timeSeries' in plotTypes:  ### NOTE!!!! --> should I do "hasbefore" and 'hasafter' as args?!??! 
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
			if customPopTitle:
				spectTitle = customPopTitle	#'CSD Spectrogram for ' + popToPlot + ', channel ' + str(spectDict['chan'])
		# plot and format 
		ax1 = plt.subplot(2, 1, 1)
		img = ax1.imshow(imshowSignal, extent=(np.amin(T), np.amax(T), minFreq, maxFreq), origin='lower', interpolation='None', aspect='auto', 
			vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap(colorMap))  # minFreq --> np.amin(F)  	# maxFreq --> np.amax(F)	# imshowSignal --> S
		divider1 = make_axes_locatable(ax1)
		cax1 = divider1.append_axes('right', size='3%', pad=0.2)
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar scientific notation
		fmt.set_powerlimits((0,0))
		cbar = plt.colorbar(img, cax = cax1, orientation='vertical', format=fmt)  # label='Power', 
		cbar.set_label(label='Power', size=17)
		cbar.ax.tick_params(labelsize=17)#12)
		ax1.set_title(spectTitle, fontsize=titleFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=T[0], right=T[1]) 			# ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq)				# Uncomment this if using the commented-out ax1.imshow (with S, and with np.amin(F) etc.)
		ax1.tick_params(axis='both', labelsize=tickFontSize)		# Increase tick size 


		### PLOT TIMESERIES ###
		# timeSeries title
		if popToPlot is None:
			timeSeriesTitle = 'CSD Signal, channel ' + str(timeSeriesDict['chan'])
		else:
			timeSeriesTitle = 'Model ' + popToPlot + ', CSD time series' #'CSD Signal for ' + popToPlot + ', channel ' + str(timeSeriesDict['chan'])

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
		ax2.tick_params(axis='both', labelsize=tickFontSize)


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
		fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)		## fmt lines are for colorbar scientific notation
		fmt.set_powerlimits((0,0))
		plt.colorbar(img, cax = cax1, orientation='vertical', label='Power', format=fmt)
		ax1.set_title(spectTitle, fontsize=titleFontSize)
		ax1.set_xlabel('Time (ms)', fontsize=labelFontSize)
		ax1.set_ylabel('Frequency (Hz)', fontsize=labelFontSize)
		ax1.set_xlim(left=T[0], right=T[1]) 			# ax1.set_xlim(left=timeRange[0], right=timeRange[1])
		# ax1.set_ylim(minFreq, maxFreq) 				# Uncomment this if using the commented-out ax1.imshow (with S, and with np.amin(F) etc.)

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

		# For potential saving 
		if popToPlot is None:
			figFilename = 'CSD_timeSeries_chan_' + str(timeSeriesDict['chan']) + '.png'
		else:
			figFilename = popToPlot + '_CSD_timeSeries_chan_' + str(timeSeriesDict['chan']) + '.png'
 
	plt.tight_layout()

	## Save figure 
	if saveFig:
		if savePath is None:
			prePath = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/figs/popContribFigs/' 		# popContribFigs_cmapJet/'
		else:
			prePath = savePath
		pathToFile = prePath + figFilename
		plt.savefig(pathToFile, dpi=600)

	## Show figure
	plt.show()
# Heatmaps for CSD data ##    NOTE: SHOULD COMBINE THIS WITH LFP DATA HEATMAP FUNCTIONS IN THE FUTURE!!
def getCSDDataFrames(dataFile, timeRange=None, oscEventInfo=None, norm=None, pops=None, verbose=0):
	#### NOTE: ADD POPS AS ARG?? 
	## This function will return data frames of peak and average CSD amplitudes, for picking cell pops
	### dataFile: str 		--> .pkl file to load, with data from the whole recording
	### oscEvenfInfo: dict 			--> dict w/ chan, left, right, minT, maxT, alignoffset
		## RIGHT NOW THIS IS *NECESSARY* TO RUN THIS SINCE NECESSARY IN getCSDdata !!!! 
	### norm: bool ---> Determines if CSD data should be normalized or not!! 
	### pops: list --> list of pops to get data frames for 
	### verbose: bool 
		### TO DO: Make evalElecs an argument, so don't have to do all of them, if that's not what we want! 

	# Load .pkl data file...? Is this necessary? Yes if I end up commenting this out for the getCSDdata function! 
	sim.load(dataFile, instantiate=False)
	# Get args for getCSDdata
	dt = sim.cfg.recordStep
	sampr = 1.0/(dt/1000.0) 	# divide by 1000.0 to turn denominator from units of ms to s
	spacing_um = 100 

	# Get all cell pops (cortical)
	thalPops = ['TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']
	ECortPops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B', 'PT5B', 'IT6', 'CT6']
	ICortPops = ['NGF1', 
			'PV2', 'SOM2', 'VIP2', 'NGF2', 
			'PV3', 'SOM3', 'VIP3', 'NGF3',
			'PV4', 'SOM4', 'VIP4', 'NGF4',
			'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',
			'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',
			'PV6', 'SOM6', 'VIP6', 'NGF6']
	allPops = list(sim.net.allPops.keys())
	if pops is None:
		pops = [pop for pop in allPops if pop not in thalPops]
	#pops = [pop for pop in allPops if pop not in thalPops] 			## exclude thal pops  ## COMMENTED OUT 11/02 TO TEST ONLY EXCITATORY CORT POPS!! 
	#pops = ICortPops  ## Changed to ICortPops from ECortPops on 11/04/22 ## ADDED IN HERE 11/02 FOR TESTING!!! 


	## Get all electrodes 	## SEE TO-DO ABOVE!! (make an if statement)
	evalElecs = []
	evalElecs.extend(list(range(int(sim.net.recXElectrode.nsites))))
	### add 'avg' to electrode list 
	evalElecs.append('avg')


	## get CSD data
	csdPopData = {}

	for pop in pops:  ## do this for ALL THE CELL POPULATIONS -- the pop selection will occur in plotting 
		print('Calculating CSD for pop ' + pop)

		csdPopData[pop] = {}
		### MAYBE CHANGE dataFile arg here?!?!?! 
		## norm:
		if norm is None:
			norm=True
		popCSDdataFULL_origShape_dict = getCSDdata(oscEventInfo=oscEventInfo, dt=dt, sampr=sampr, dataFile=dataFile, pop=pop, norm=norm, spacing_um=spacing_um, outputType=[]) ## HAVE TO ADD oscEventInfo arg since right now this arg is * NECESSARY * TO RUN getCSDdata!!!! 	# popCSDdataFULL_origShape = getCSDdata(dataFile, pop=pop) 
		popCSDdataFULL_origShape = popCSDdataFULL_origShape_dict['csdData']#['csd']
		popCSDdataFULL = np.transpose(popCSDdataFULL_origShape)	### TRANSPOSE THIS so (20,230000) --> (230000, 20)

		if timeRange is None:
			popCSDdata = popCSDdataFULL.copy()
		else:
			popCSDdata = popCSDdataFULL[int(timeRange[0]/sim.cfg.recordStep):int(timeRange[1]/sim.cfg.recordStep),:]


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
			avgValues['avgCSD'][q].append(csdPopData[pop][elecKey]['avg'])
		q+=1
	dfAvg = pd.DataFrame(avgValues['avgCSD'], index=pops)


	# return csdPopData
	if verbose:
		return dfPeak, dfAvg, peakValues, avgValues, csdPopData 
	else:
		return dfPeak, dfAvg

## PSD: data and plotting ## 
def getPSDdata(dataFile, inputData, inputDataType='spectrogram', duringOsc=1, minFreq=1, maxFreq=110, stepFreq=0.25):
	## Look at the power spectral density of a given data set (e.g. CSD, LFP, summed LFP, etc.)
	### dataFile:str 				--> .pkl file with simulation recording 
	### inputData 					--> data to be analyzed 
	### inputDataType: str			--> 'spectrogram' or 'timeSeries'
	### duringOsc: bool 			--> evaluating during an oscillation event or no? 
	### minFreq
	### maxFreq
	### stepFreq

	## TO DO: 
		## Go thru the comments for smaller to-do tasks
		## Make this generalizable for CSD and LFP -- figure out more efficient way to do this
		## Be able to extract minFreq and maxFreq and stepFreq from extant spectrogram data input 

	# load simulation .pkl file 
	sim.load(dataFile, instantiate=False)  ## Loading this just to get the sim.cfg.recordStep !! 


	# Set up output data dictionary
	psdData = {}

	#########################################
	#### Morlet wavelet transform method ####
	#########################################

	allFreqs = []
	allSignal = []

	Fs = int(1000.0/sim.cfg.recordStep)

	## Calculations for raw LFP or CSD timeSeries
	if inputDataType == 'timeSeries':
		# freqList = np.arange(minFreq, maxFreq+stepFreq, stepFreq)
		morletSpec = MorletSpec(inputData, Fs, freqmin=minFreq, freqmax=maxFreq, freqstep=stepFreq)#, lfreq=freqList) # specDuring[0]
		freqs = F = morletSpec.f 		# F_during = specDuring[0].f
		spec = morletSpec.TFR 			# S_during = specDuring[0].TFR
		psdData.update({'spec': spec})
		signal = np.mean(spec, axis=1)  
		# signal = mednorm(spec) # 		lmsnorm = [mednorm(ms.TFR) for ms in lms]




	## Calculations for spectrogram input 
	elif inputDataType == 'spectrogram':
		## This assumes a dict w/ an output structure that matches getCSDdata output when outputType is spectrogram
		if duringOsc:
			freqs = inputData['F_during']
			spec = inputData['S_during']
			signal = np.mean(spec, axis=1)
		elif not duringOsc:
			freqs = inputData['F_full']
			spec = inputData['S_full']
			signal = np.mean(spec, axis=1)


	allFreqs.append(freqs)
	allSignal.append(signal)

	# put allFreqs and allSignal into output data dictionary 
	psdData.update({'allFreqs': allFreqs, 'allSignal': allSignal})


	### To print out frequency w/ max power:  ## TO DO: BETTER NAMING 
	maxSignalIndex_unformatted = np.where(signal==np.amax(signal))
	x = list(maxSignalIndex_unformatted[0])
	maxSignalIndex = x[0]
	maxPowerFrequency = freqs[maxSignalIndex]
	psdData.update({'maxSignalIndex': maxSignalIndex, 'maxPowerFrequency': maxPowerFrequency})

	print('max power frequency in signal: ' + str(maxPowerFrequency))

	return psdData 
def plotPSD(psdData, minFreq=1, maxFreq=110, freqStep=0.25, lineWidth=1.0, fontSize=12, color='k', figSize=(10,7)):
	### 	----> NOTE: MAKE OVERLAY OPTION POSSIBLE? 
	### This function should plot the PSD data 
	### psdData -->  output of def getPSD()
	### maxFreq --> 
	### lineWidth --> 
	### fontSize --> 
	### color --> Default: 'k' (black)
	### figSize --> Default: (10,7)

	# Get signal & frequency data
	signalList = psdData['allSignal']
	signal = signalList[0]

	freqsList = psdData['allFreqs']
	freqs = freqsList[0]
	if type(freqs) is not list:
		freqs = list(freqs)


	freqsToPlot = [freq for freq in freqs if freq >= minFreq and freq <= maxFreq]
	signalToPlot = signal[freqs.index(freqsToPlot[0]):freqs.index(freqsToPlot[-1])+1]		# signalToPlot = signal[freqs.index(minFreq):freqs.index(maxFreq)+1]
	if len(freqsToPlot) == len(signalToPlot):
		plt.figure(figsize=figSize)
		plt.plot(freqsToPlot, signalToPlot, linewidth=lineWidth, color=color)				# plt.plot(freqs[freqs<maxFreq], signal[freqs<maxFreq], linewidth=lineWidth, color=color)


	# format plot
	plt.xlim([minFreq, maxFreq])
	plt.xticks(np.arange(minFreq, maxFreq, step=freqStep))
	plt.xlabel('Frequency (Hz)', fontsize=fontSize)
	ylabel = 'Power'
	plt.ylabel(ylabel, fontsize=fontSize)
	plt.tight_layout()
	plt.suptitle('Power Spectral Density', fontsize=fontSize, fontweight='bold')
	plt.show()

## peakF calculations ## 
def getPeakF(dataFile, inputData, csdAllChans, timeData=None, chan=8, freqmin=1.0, freqmax=110, freqStep=0.25, plotTest=True, plotNorm=0): #  lchan=None, 
	### This function will calculate the peakF of a given dataset (oscEvent) using OEvent methodology (see load.py)
	## dataFile: .pkl file 
	## inputData: list / array with data to be analyzed
	## csdAllChans: list / array of csd data with all chans, over same timepoints as inputData
	## timeData
	## chan: int 
	## freqmin: float			--> default: 0.25 (see load.py)
	## freqmax: float			--> default: 110  (see load.py)
	## freqStep: float 			--> default: 0.25 (see load.py)
	## plotTest: bool 			--> test if correct image is being evaluated by plotting 


	print('Getting peakF using OEvent methods')

	# Establish output data dictionary 
	outputData = {}


	# Load dataFile to determine sampling rate
	if dataFile:
		sim.load(dataFile, instantiate=False)		## load .pkl simulation file 
	else:
		print('No dataFile; will use data from dataFile already loaded elsewhere!')

	# Calculate sampling rate from simulation timestep 
	dt = sim.cfg.recordStep  	# or should I divide by 1000.0 up here, and then just do 1.0/dt below for sampr? 
	sampr = 1.0/(dt/1000.0) 	# divide by 1000.0 to turn denominator from units of ms to s

	# Argument values for getmorletwin() below 
	winsz = 10 #565 #10 				# window size
	getPhase = True
	noiseampCSD = 200.0 / 10.0 # amplitude cutoff for CSD noise; was 200 before units fix
	noiseamp=noiseampCSD 
	useloglfreq=False
	mspecwidth=7.0

	# This line is from getIEIstatsbyBand in load.py 
	## "get morlet specgrams on windows of dat time series (window size in samples = winsz)"
	lms,lnoise,lsidx,leidx = getmorletwin(inputData,int(winsz*sampr),sampr,freqmin=freqmin,
		freqmax=freqmax,freqstep=freqStep,getphase=getPhase,useloglfreq=useloglfreq,mspecwidth=mspecwidth,
		noiseamp=noiseamp) # inputData <-- dat[chan, :] OR dat[:, chan]  # dat[:,chan]

	# msn <-- lmsnorm # This is from getIEIstatsbyBand, where normop == mednorm (per the arguments) 	## lmsnorm = [normop(ms.TFR) for ms in lms]
	lmsnorm = [mednorm(ms.TFR) for ms in lms] ##unitnorm or mednorm 


	# Update output data dictionary with output from getmorletwin & lmsnorm 
	outputData.update({'lms': lms, 'lmsnorm': lmsnorm, 'lnoise': lnoise, 'lsidx': lsidx, 'leidx': leidx})

	## Trying getspecevents line 
	evthresh = medthresh = 4.0
	MUA=None
	chan = chan
	overlapth = 0.5
	getphase=True  ## TRYING FALSE FOR NOW BC of all the issues with getspecevents -- don't feel like debugging now 
	endfctr=0.5
	# sig --> full CSD data (not selected for chan)

	## for normalized spectrogram
	llevent_norm = getspecevents_norm(lms,lmsnorm,lnoise,evthresh,lsidx,leidx,inputData,MUA,chan,sampr,overlapth=overlapth,getphase=getphase,endfctr=endfctr) # get the spectral events
	outputData.update({'llevent_norm': llevent_norm})

	## for non-normalized spectrogram 
	llevent = getspecevents_nonNorm(lms,lmsnorm,lnoise,evthresh,lsidx,leidx,inputData,MUA,chan,sampr,overlapth=overlapth,getphase=getphase,endfctr=endfctr)
	outputData.update({'llevent': llevent})

	############################################


	# Argument values for getblobsfrompeaks()
	medthresh = 4.0 # 20 #4.0
	endfctr = 0.5

	# These lines are from getspecevents (which is called in getIEIstatsbyBand)
	for windowidx,offidx,ms,msn,noise in zip(np.arange(len(lms)),lsidx,lms,lmsnorm,lnoise): 
		imgpk = detectpeaks(msn) # detect the 2D local maxima ### NOTE: This should probably give me *ONE* peak!! 
		imgpk_nonNorm = detectpeaks(ms.TFR)
		print('imgpk detected')
		lblob_norm = getblobsfrompeaks(msn,imgpk,ms.TFR,medthresh,endfctr=endfctr,T=ms.t,F=ms.f) # cut out the blobs/events
		lblob_nonNorm = getblobsfrompeaks(ms.TFR,imgpk_nonNorm,ms.TFR,medthresh,endfctr=endfctr,T=ms.t,F=ms.f) # cut out the blobs/events
												
		print('lblob gotten')

	outputData.update({'imgpk_nonNorm': imgpk_nonNorm, 'lblob_nonNorm': lblob_nonNorm})
	outputData.update({'imgpk': imgpk, 'lblob_norm': lblob_norm})


	## TEST PLOTTING
	if plotTest:

		fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,7)) # fig = plt.figure() # figsize=figSize
		ax1 = plt.subplot(1, 1, 1)

		# vc = [1,80]
		# vmin = np.array([s.TFR for s in lms]).min()
		# vmax = np.array([s.TFR for s in lms]).max()
		# vc = [vmin, vmax]
		# # print('vmin: ' + str(vmin))
		# # print('vmax: ' + str(vmax))

		
		if plotNorm:
			imshowSignalNorm = lmsnorm[0] # .TFR  
			imshowSignal=imshowSignalNorm
			vmin = np.array([s.TFR for s in lms][1]).min()
			vmax = np.array([s.TFR for s in lms][1]).max()
			vc = [vmin, vmax]
		else:
			S = lms[0].TFR
			imshowSignal = S # imshowSignal 	# S #### [freqmin:freqmax+1]		# S[minFreq:maxFreq+1]
			vmin = np.array([s.TFR for s in lms][0]).min()
			vmax = np.array([s.TFR for s in lms][0]).max()
			vc = [vmin, vmax]


		T = lms[0].t # timeData # ms.t # lms[0].t ### CORRECT USING MINT, ALIGNOFFSET, ETC. 
		F = lms[0].f

		img = ax1.imshow(imshowSignal, extent=(np.amin(T), np.amax(T), np.amin(F), np.amax(F)), origin='lower', interpolation='None', aspect='auto', 
				vmin=vc[0], vmax=vc[1], cmap=plt.get_cmap('jet')) # np.amin(f) <--freqmin <-- minFreq # np.amax(f) <-- freqmax <-- maxFreq # 'jet' <-- colorMap


		## Overlay peak points
		if plotNorm:
			peaks = np.where(imgpk==True)
			plt.title('Normalized Spectrogram w/ peak points overlaid')
		else:
			peaks = np.where(imgpk_nonNorm==True)
			plt.title('Non-normalized Spectrogram w/ peak points overlaid')

		# time x-axis point
		time_inds = peaks[1]
		# spectrogram y-axis point
		spec_inds = peaks[0]

		x = []
		y = []
		for i in range(len(peaks[0])):
			if T[time_inds[i]] > np.amin(T) and T[time_inds[i]] < np.amax(T) and F[spec_inds[i]] > np.amin(F) and F[spec_inds[i]] < np.amax(F):
				x.append(T[time_inds[i]])
				y.append(F[spec_inds[i]])

		plt.scatter(x,y, c='yellow', s=5)


		plt.show()

	return outputData
















