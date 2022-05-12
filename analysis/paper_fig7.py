import os
import pandas as pd
import matplotlib.pyplot as plt


def getSubjectNames(dataPath, sim):
	### dataPath: str --> path to directory with wavelet files
	### sim: bool --> to decide which prefix to use for sorting 

	os.chdir(dataPath)
	dataDirContents = os.listdir()

	subjects = []

	if sim:
		prefix = 'v34_batch57_'
	else:
		prefix = '2-'

	for obj in dataDirContents:
		if prefix in obj and os.path.isdir(obj):
			subjects.append(obj)

	return subjects
def getStats(df, evidx,align='bywaveletpeak',verbose=False):
	dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore = geteventprop(df, evidx, align)   #= self.getallprop(evidx,align) 

	if verbose:
		return dur,chan,hasbefore,hasafter,windowidx,offidx,left,right,minT,maxT,peakT,minF,maxF,peakF,avgpowevent,ncycle,WavePeakT,WaveTroughT,WaveletPeakT,WaveletLeftTroughT,WaveletRightTroughT,w2,left,right,band,alignoffset,filtsigcor,Foct,cycnpeak,ERPscore,OSCscore
	else:
		return dur,peakF,ncycle,band
def getPropDicts(subjects, frequencyBands, dataPath, propType='all'): 
	### subjects: list
	### layers: dict
	### frequencyBands: list
	### dataPath: str
	### propType: str ('dur', 'peakF', 'ncycles', 'all')
	### timeRange: bool 

	durDict = {}   		# establish a dict for wavelet duration values 
	peakFDict = {}		# establish a dict for wavelet peakF values 
	ncycleDict = {}		# establish a dict for wavelet ncycles values 

	regions = ['supra', 'gran', 'infra']

	for subject in subjects:
		durDict[subject] = {}
		peakFDict[subject] = {}
		ncycleDict[subject] = {}

		for region in regions:
			durDict[subject][region] = {}
			peakFDict[subject][region] = {}
			ncycleDict[subject][region] = {}

			for band in frequencyBands:
				durDict[subject][region][band] = []
				peakFDict[subject][region][band] = []
				ncycleDict[subject][region][band] = []

				allRegionFiles = os.listdir(dataPath + subject + '/' + region)
				pklFiles = []
				for file in allRegionFiles:
					if '.pkl' in file:
						pklFiles.append(file)

				for pklFile in pklFiles:
					pklBand = pklFile.split('_')[-1][:-4]
					# for band in frequencyBands:
					if band == pklBand:  	# this prevents gamma / hgamma double-counting
						dfsPkl = pd.read_pickle(dataPath + subject + '/' + region + '/' + pklFile)
						for idx in dfsPkl.index:
							dur, peakF, ncycle, band = getStats(df = dfsPkl, evidx = idx, align='bywaveletpeak') 
							durDict[subject][region][band].append(dur)
							peakFDict[subject][region][band].append(peakF)
							ncycleDict[subject][region][band].append(ncycle)


	if propType == 'all':
		return durDict, peakFDict, ncycleDict
	elif propType == 'dur':
		return durDict
	elif propType == 'peakF':
		return peakFDict
	elif propType == 'ncycle':
		return ncycleDict
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
	simPath = dataPathPrefix + 'v34_batch57/'
	nhpPath = dataPathPrefix + 'NHP_data/spont/'

	simSubjects = getSubjectNames(dataPath=simPath, sim=1)
	nhpSubjects = getSubjectNames(dataPath=nhpPath, sim=0)

	# print('simSubjects: ' + str(simSubjects)) ## TESTING LINES --> WORKS
	# print('nhpSubjects: ' + str(nhpSubjects)) ## TESTING LINES --> WORKS


	### Frequency bands & region ###
	frequencyBands = ['delta', 'theta', 'alpha', 'beta', 'gamma']
	regions = ['supra', 'gran', 'infra', 'all']  



	###########################
	#### SIM WAVELET STATS ####
	###########################

	durDictSim, peakFDictSim, ncycleDictSim = getPropDicts(simSubjects, frequencyBands, simPath, propType='all')

	durSimLists = getPropLists(durDictSim, regions, frequencyBands)
	peakFSimLists = getPropLists(peakFDictSim, regions, frequencyBands)
	ncycleSimLists = getPropLists(ncycleDictSim, regions, frequencyBands)

	simListsDict = {'dur': durSimLists, 'peakF': peakFSimLists, 'nCycle': ncycleSimLists}



	# ##########################
	# ### NHP WAVELET STATS ####
	# ##########################

	durDictNHP, peakFDictNHP, ncycleDictNHP = getPropDicts(nhpSubjects, frequencyBands, nhpPath, propType='all')

	durNHPLists = getPropLists(durDictNHP, regions, frequencyBands)
	peakFNHPLists = getPropLists(peakFDictNHP, regions, frequencyBands)
	ncycleNHPLists = getPropLists(ncycleDictNHP, regions, frequencyBands)

	nhpListsDict = {'dur': durNHPLists, 'peakF': peakFNHPLists, 'nCycle': ncycleNHPLists}


	#### GENERATE BOXPLOTS ####
	statsBoxplotALL(frequencyBands, simListsDict=simListsDict, nhpListsDict=nhpListsDict, dataCategories=['dur', 'peakF', 'nCycle'], figsize=None, colors=None)

######### PATH TO DATA DIRECTORY #########
dataPathPrefix = '../data/'  ## NOTE: Change this to wherever Macaque_auditory_thalamocortical_model_data/data directory is located, e.g. '/Home/Desktop/Macaque_auditory_thalamocortical_model_data/data/'


# --------------------------
# Main
# --------------------------
if __name__ == '__main__':
	# Fig 7 
	plotStats()





