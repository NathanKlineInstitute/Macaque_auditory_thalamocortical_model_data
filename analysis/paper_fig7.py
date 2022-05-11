import os
import pandas as pd
import matplotlib.pyplot as plt


def sortSubjectFiles(waveletPath, sim):
	### waveletPath: str --> path to directory with wavelet files
	### sim: bool --> to decide which prefix to use for sorting 

	os.chdir(waveletPath)
	subjectNames = os.listdir()

	subjFiles = []

	if sim:
		filePrefix = 'v34_batch57_'
		splitStr = '_data'
	else:
		filePrefix = '2-'
		splitStr = '_'

	for name in subjectNames:
		if filePrefix in name:
			subjFiles.append(name)

	subjFiles.sort()
	subjSubjects = {}
	for subjFile in subjFiles:
		subjName = subjFile.split(splitStr)[0]
		if subjName not in subjSubjects.keys():
			subjSubjects[subjName] = []
		subjSubjects[subjName].append(subjFile)

	return subjSubjects
def getSimLayers():
	layers = {'supra':[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 'gran':[10, 11], 'infra':[12, 13, 14, 15, 16, 17, 18, 19]} ## NETPYNE SIMS
	return layers
def getNHPLayers(waveletPath, subjects):
	### waveletPath: str
	### subjects: list or dict (probably dict) -- but timeRange == 1 if dict, 0 if not 

	layers = {} # {'supra': [], 'gran': [], 'infra': []}		# have to reset this for each NHP subject 
	regions = ['supra', 'gran', 'infra']

	if type(subjects) is dict:
		for subject in subjects:
			layers[subject] = {}

			for region in regions:
				layers[subject][region] = []

			channels = []
			allSubFiles = os.listdir(waveletPath + subjects[subject][0])
			for n in allSubFiles:
				if 'chan_' in n:
					channels.append(int(n.split('_')[-1]))
			channels.sort()

			for layerKey in layers[subject]:
				if layerKey == 'supra':
					layers[subject][layerKey].append(channels[0])
					layers[subject][layerKey].append(channels[1])
				elif layerKey == 'gran':
					layers[subject][layerKey].append(channels[2])
				elif layerKey == 'infra':
					layers[subject][layerKey].append(channels[-2])
					layers[subject][layerKey].append(channels[-1])

	# elif type(subjects) is list:  # FILL THIS IN

	return layers 
def getStats(df, evidx,align='bywaveletpeak',verbose=False):
	dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore = geteventprop(df, evidx, align)   #= self.getallprop(evidx,align) 

	if verbose:
		return dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore
	else:
		return dur,peakF,ncycle,band
def getPropDicts(subjects, frequencyBands, waveletPath, propType='all', sim=1):
	### subjects: list (if no timeRange) or dict (if broken up into timeRanges)
	### layers: dict
	### frequencyBands: list
	### waveletPath: str
	### propType: str ('dur', 'peakF', 'ncycles', 'all')
	### timeRange: bool 

	durDict = {}   		# establish a dict for wavelet duration values 
	peakFDict = {}		# establish a dict for wavelet peakF values 
	ncycleDict = {}		# establish a dict for wavelet ncycles values 

	if sim:  # NetPyNE sim data
		layers = getSimLayers()

		if type(subjects) is list:
			for subject in subjects:
				durDict[subject] = {}
				peakFDict[subject] = {}
				ncycleDict[subject] = {}

				for layerKey in layers:
					durDict[subject][layerKey] = {}
					peakFDict[subject][layerKey] = {}
					ncycleDict[subject][layerKey] = {}

					for band in frequencyBands:
						durDict[subject][layerKey][band] = []
						peakFDict[subject][layerKey][band] = []
						ncycleDict[subject][layerKey][band] = []

					for chan in layers[layerKey]:
						allChanFiles = os.listdir(waveletPath + subject + '/chan_' + str(chan))
						pklFiles = []
						for file in allChanFiles:
							if '.pkl' in file:
								pklFiles.append(file)   ## now pklFiles is a list of all the .pkl files in a particular chan subdir 
						#print(pklFiles)

						for pklFile in pklFiles:
							#print(pklFile)  # ---> print line for testing
							pklBand = pklFile.split('_')[-1][:-4]
							for band in frequencyBands:
								if band == pklBand:  	# this prevents gamma / hgamma double-counting (excluding hgamma for now!)
									dfsPkl = pd.read_pickle(waveletPath + subject + '/chan_' + str(chan) + '/' + pklFile)
									#print('band = ' + str(band)) # ---> print line for testing
									#print('chan = ' + str(chan)) # ---> print line for testing
									for idx in dfsPkl.index:
										dur, peakF, ncycle, band = getStats(df = dfsPkl, evidx = idx, align='bywaveletpeak') 
										durDict[subject][layerKey][band].append(dur)
										peakFDict[subject][layerKey][band].append(peakF)
										ncycleDict[subject][layerKey][band].append(ncycle)

		elif type(subjects) is dict:
			for subject in subjects:   # here subjects is a dict, e.g. nhpSubjects = {'2-rb045046026': ['2-rb045046026_timeRange_0_40', '2-rb045046026_timeRange_120_160', '2-rb045046026_timeRange_160_200']}
				durDict[subject] = {}
				peakFDict[subject] = {}
				ncycleDict[subject] = {}

				for layerKey in layers:
					durDict[subject][layerKey] = {}
					peakFDict[subject][layerKey] = {}
					ncycleDict[subject][layerKey] = {}

					for band in frequencyBands:
						durDict[subject][layerKey][band] = []
						peakFDict[subject][layerKey][band] = []
						ncycleDict[subject][layerKey][band] = []

					## Now get pklFiles from each channel 
					for chan in layers[layerKey]:
						for tRange in subjects[subject]:
							# print(timeRange)  # <-- print line for testing 
							allChanFiles = os.listdir(waveletPath + tRange + '/chan_' + str(chan))
							pklFiles = []
							for file in allChanFiles:
								if '.pkl' in file:
									pklFiles.append(file)   ## now pklFiles is a list of all the .pkl files in a particular chan subdir 
							#print(pklFiles)		# <-- print line for testing 

							for pklFile in pklFiles:
								pklBand = pklFile.split('_')[-1][:-4]
								for band in frequencyBands:
									if band == pklBand:  	# this prevents gamma / hgamma double-counting (excluding hgamma for now!)
										dfsPkl = pd.read_pickle(waveletPath + tRange + '/chan_' + str(chan) + '/' + pklFile)
										#print('band = ' + str(band)) # ---> print line for testing
										#print('chan = ' + str(chan)) # ---> print line for testing
										for idx in dfsPkl.index:
											dur, peakF, ncycle, band = getStats(df = dfsPkl, evidx = idx, align='bywaveletpeak') 
											durDict[subject][layerKey][band].append(dur)
											peakFDict[subject][layerKey][band].append(peakF)
											ncycleDict[subject][layerKey][band].append(ncycle)

	else:   # NHP data 
		if type(subjects) is dict:
			layers = getNHPLayers(waveletPath,subjects)

			for subject in subjects:
				durDict[subject] = {}
				peakFDict[subject] = {}
				ncycleDict[subject] = {}

				for layerKey in layers[subject]:
					durDict[subject][layerKey] = {}
					peakFDict[subject][layerKey] = {}
					ncycleDict[subject][layerKey] = {}

					for band in frequencyBands:
						durDict[subject][layerKey][band] = []
						peakFDict[subject][layerKey][band] = []
						ncycleDict[subject][layerKey][band] = []

					## Now get pklFiles from each channel 
					for chan in layers[subject][layerKey]:
						for tRange in subjects[subject]:
							# print(timeRange)  # <-- print line for testing 
							allChanFiles = os.listdir(waveletPath + tRange + '/chan_' + str(chan))
							pklFiles = []
							for file in allChanFiles:
								if '.pkl' in file:
									pklFiles.append(file) 

							for pklFile in pklFiles:
								pklBand = pklFile.split('_')[-1][:-4]
								for band in frequencyBands:
									if band == pklBand:  	# this prevents gamma / hgamma double-counting (excluding hgamma for now!)
										dfsPkl = pd.read_pickle(waveletPath + tRange + '/chan_' + str(chan) + '/' + pklFile)
										#print('band = ' + str(band)) # ---> print line for testing
										#print('chan = ' + str(chan)) # ---> print line for testing
										for idx in dfsPkl.index:
											dur, peakF, ncycle, band = getStats(df = dfsPkl, evidx = idx, align='bywaveletpeak') 
											durDict[subject][layerKey][band].append(dur)
											peakFDict[subject][layerKey][band].append(peakF)
											ncycleDict[subject][layerKey][band].append(ncycle)


		elif type(subjects) is list:
			print('THIS CONDITION IS NOT COMPLETED YET!!!')


	if sim:
		durDictSim = durDict
		peakFDictSim = peakFDict
		ncycleDictSim = ncycleDict
		if propType == 'all':
			return durDictSim, peakFDictSim, ncycleDictSim
		elif propType == 'dur':
			return durDictSim
		elif propType == 'peakF':
			return peakFDictSim
		elif propType == 'ncycle':
			return ncycleDictSim 
	else:
		durDictNHP = durDict
		peakFDictNHP = peakFDict
		ncycleDictNHP = ncycleDict
		if propType == 'all':
			return durDictNHP, peakFDictNHP, ncycleDictNHP
		elif propType == 'dur':
			return durDictNHP
		elif propType == 'peakF':
			return peakFDictNHP
		elif propType == 'ncycle':
			return ncycleDictNHP 
def getPropStats(propDict, region, frequencyBand, avg=1):
	### propDict: dict
	### region: str --> 'all', 'supra', 'gran', 'infra'
	### frequencyBand: str --> 'alpha', 'beta', 'delta', 'theta', 'gamma'
	### avg: bool 

	if region is not 'all':
		prop = []
		for subject in propDict:
			prop.append(propDict[subject][region][frequencyBand])

	else:
		prop = []
		for subject in propDict:
			prop.append(propDict[subject]['supra'][frequencyBand])
			prop.append(propDict[subject]['gran'][frequencyBand])
			prop.append(propDict[subject]['infra'][frequencyBand])

	propFlat = [item for subList in prop for item in subList]
	avgProp = sum(propFlat) / len(propFlat)

	if avg:
		return avgProp
	else:
		return propFlat 
def getPropLists(propDict, regions, frequencyBands): 
	### propDict: dict
	### regions: list --> e.g. ['supra', 'gran', 'infra']
	### frequencyBands: list --> e.g. ['alpha', 'beta', 'delta', 'theta', 'gamma']
	### avg: bool 


	propLists = {}

	for band in frequencyBands:
		propLists[band] = {}
		for region in regions:
			propLists[band][region] = getPropStats(propDict, region, band, 0)

	return propLists
def geteventprop (dframe,evidx,align):
	# get all major event properties, used for drawing the event or other...
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
def statsBoxplotALL(frequencyBands, simListsDict, nhpListsDict, dataCategories, figsize=None, colors=None):  
	### frequencyBands: list, e.g. ['alpha', 'beta']
	### simListsDict: dict, e.g. ____
	### nhpListsDict: dict, e.g. ____
	### dataCategories: list --> ['dur', 'peakF', 'nCycle']
	### titles: list of strs, e.g. ['WAVELET DURATION', 'NUM CYCLES'] <-- hmm probably change this 
	### figsize: tuple, e.g. (11,7)
	### colors: list, e.g. ['blue', 'lightgreen']

	if figsize is None:
		figsize = (13,7) #(11,7)
	fig = plt.figure(figsize=figsize)

	if colors is None:
		# colors = ['yellow','lightblue']
		colorsDur = ['yellow','lightblue']	#['purple', 'green']
		colorsPeakF = ['purple', 'green']
		colorsNCycle = ['red', 'blue']

	nrows = len(dataCategories)  ### should make it 
	ncols = len(frequencyBands)

	i=1
	for category in dataCategories: 
		r=1
		for band in frequencyBands:
			## Lists & colors of boxes ##
			if category == 'dur':
				colors = colorsDur
				simLists = simListsDict['dur']
				nhpLists = nhpListsDict['dur']
				cat_yLabel = 'DURATION\n(ms)'
			elif category == 'peakF':
				colors = colorsPeakF
				simLists = simListsDict['peakF']
				nhpLists = nhpListsDict['peakF']
				cat_yLabel = 'PEAK\nFREQUENCY\n(Hz)'
			elif category == 'nCycle':
				colors = colorsNCycle
				simLists = simListsDict['nCycle']
				nhpLists = nhpListsDict['nCycle']
				cat_yLabel = 'NUM CYCLES'

			ax = fig.add_subplot(nrows, ncols, i)
			simDataPlot = simLists[band]['all'] 			# simLists[band][region]
			nhpDataPlot = nhpLists[band]['all']				# nhpLists[band][region]
			bp = ax.boxplot((simDataPlot,nhpDataPlot),patch_artist=True)#,showfliers=False)


			for patch, color in zip(bp['boxes'], colors):
				patch.set_facecolor(color)
			## style of fliers ##
			for flier in bp['fliers']:
				flier.set(marker='o', markersize=3)
			## axes and titles ##
			ax.set_xticklabels(['MODEL', 'NHP'], fontsize=9)	# (['SIM', 'NHP'], fontsize=9)		# ax.xaxis.set_visible(False)
			ax.tick_params(axis='y', labelsize=6)
			if r:
				ax.set_ylabel(cat_yLabel,labelpad=40,rotation=0,fontsize=10,verticalalignment='center',fontweight='bold')
			if category == dataCategories[0]: 
				ax.set_title(band, fontsize=10, fontweight='bold')
			i+=1
			r=0

	plt.subplots_adjust(top=0.85, bottom=0.05, wspace=0.3, hspace=0.3)
	fig.suptitle('COMPARISON OF OSCILLATION EVENT PROPERTIES', fontsize=14, fontweight='bold', horizontalalignment='center', y=0.95)

	plt.show()
def plotStats():
	simSubjects = sortSubjectFiles(waveletPath=waveletPath, sim=1)
	nhpSubjects = sortSubjectFiles(waveletPath=waveletPath, sim=0)

	### Frequency bands & region ###
	frequencyBands = ['delta', 'theta', 'alpha', 'beta', 'gamma']	# , 'hgamma'] ## DON'T DO HGAMMA FOR NOW 
	regions = ['supra', 'gran', 'infra', 'all']  


	###########################
	#### SIM WAVELET STATS ####
	###########################

	durDictSim, peakFDictSim, ncycleDictSim = getPropDicts(simSubjects, frequencyBands, waveletPath, propType='all', sim=1)

	durSimLists = getPropLists(durDictSim, regions, frequencyBands)
	peakFSimLists = getPropLists(peakFDictSim, regions, frequencyBands)
	ncycleSimLists = getPropLists(ncycleDictSim, regions, frequencyBands)

	simListsDict = {'dur': durSimLists, 'peakF': peakFSimLists, 'nCycle': ncycleSimLists}



	##########################
	### NHP WAVELET STATS ####
	##########################

	durDictNHP, peakFDictNHP, ncycleDictNHP = getPropDicts(nhpSubjects, frequencyBands, waveletPath, propType='all', sim=0)

	durNHPLists = getPropLists(durDictNHP, regions, frequencyBands)
	peakFNHPLists = getPropLists(peakFDictNHP, regions, frequencyBands)
	ncycleNHPLists = getPropLists(ncycleDictNHP, regions, frequencyBands)

	nhpListsDict = {'dur': durNHPLists, 'peakF': peakFNHPLists, 'nCycle': ncycleNHPLists}



	#### GENERATE BOXPLOTS ####
	statsBoxplotALL(frequencyBands, simListsDict=simListsDict, nhpListsDict=nhpListsDict, dataCategories=['dur', 'peakF', 'nCycle'], figsize=None, colors=None)  ### <-- regions can take 'all' or no? 


######### PATH TO DATA DIRECTORY #########
waveletPath = '../data/'  ## NOTE: Change this to wherever Macaque_auditory_thalamocortical_model_data/data directory is located 


# --------------------------
# Main
# --------------------------
if __name__ == '__main__':
	# Fig 7 
	plotStats()





