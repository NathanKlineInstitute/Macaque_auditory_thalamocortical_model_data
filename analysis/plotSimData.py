from simDataAnalysis import *
import json
import os
import pickle 

##########################
#### USEFUL VARIABLES ####
##########################
## set layer bounds:
layerBounds= {'L1': 100, 'L2': 160, 'L3': 950, 'L4': 1250, 'L5A': 1334, 'L5B': 1550, 'L6': 2000}

## Cell populations: 
allpops = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4',
'PV4', 'SOM4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B',
'SOM5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']#, 'IC']
L1pops = ['NGF1']
L2pops = ['IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2']
L3pops = ['IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3']
L4pops = ['ITP4', 'ITS4', 'PV4', 'SOM4', 'VIP4', 'NGF4']
L5Apops = ['IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A']
L5Bpops = ['IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B']
L6pops = ['IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6']
thalPops = ['TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']
ECortPops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B', 'PT5B', 'IT6', 'CT6']
ICortPops = ['NGF1', 
			'PV2', 'SOM2', 'VIP2', 'NGF2', 
			'PV3', 'SOM3', 'VIP3', 'NGF3',
			'PV4', 'SOM4', 'VIP4', 'NGF4',
			'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',
			'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',
			'PV6', 'SOM6', 'VIP6', 'NGF6']
AllCortPops = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4',
'PV4', 'SOM4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B',
'SOM5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6']
EThalPops = ['TC', 'TCM', 'HTC']				# TEpops = ['TC', 'TCM', 'HTC']
IThalPops = ['IRE', 'IREM', 'TI', 'TIM']		# TIpops = ['IRE', 'IREM', 'TI', 'TIM']
reticPops = ['IRE', 'IREM']
matrixPops = ['TCM', 'TIM', 'IREM']
corePops = ['TC', 'HTC', 'TI', 'IRE']


## electrodes <--> layer
L1electrodes = [0]
L2electrodes = [1]
L3electrodes = [2,3,4,5,6,7,8,9]
L4electrodes = [9,10,11,12]
L5Aelectrodes = [12,13]
L5Belectrodes = [13,14,15]
L6electrodes = [15,16,18,19]
allElectrodes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
supra = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
gran = [9, 10, 11, 12]
infra = [12, 13, 14, 15, 16, 17, 18, 19]

## for COLORS in DATA PLOTTING!! 
colorList = [[0.42,0.67,0.84], [0.90,0.76,0.00], [0.42,0.83,0.59], [0.90,0.32,0.00],
			[0.34,0.67,0.67], [0.90,0.59,0.00], [0.42,0.82,0.83], [1.00,0.85,0.00],
			[0.33,0.67,0.47], [1.00,0.38,0.60], [0.57,0.67,0.33], [0.5,0.2,0.0],
			[0.71,0.82,0.41], [0.0,0.2,0.5], [0.70,0.32,0.10]]*3

colorDict = {}
for p in range(len(allpops)):
	colorDict[allpops[p]] = colorList[p]




##########################################
##### LOOKING AT INDIVIDUAL WAVELETS ##### 
##########################################

## Set base directory for run data files (not wavelet data files!)
local = 1
if local:
	based = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/' 


##  Set wavelet example frequency band
delta = 0
beta = 	0
alpha = 0
theta = 0
# gamma = 0 

if delta:
	timeRange, dataFile, waveletElectrode = getWaveletInfo('delta', based=based, verbose=1)
	wavelet='delta' ### MOVE THESE EVENTUALLY -- BEING USED FOR peakTitle
	# ylim = [1, 40]
	# maxFreq = ylim[1]  ## maxFreq is for use in plotCombinedLFP, for the spectrogram plot 
elif beta:
	timeRange, dataFile, waveletElectrode = getWaveletInfo('beta', based=based, verbose=1)
	wavelet='beta'
	# maxFreq=None
elif alpha:
	timeRange, dataFile, waveletElectrode = getWaveletInfo('alpha', based=based, verbose=1) ## recall timeRange issue (see nb)
	wavelet='alpha'
	# maxFreq=None
elif theta:
	# timeRange, dataFile = getWaveletInfo('theta', based)
	timeRange, dataFile, waveletElectrode = getWaveletInfo('theta', based=based, verbose=1)
	wavelet='theta'
	# maxFreq=None
# elif gamma:
# 	print('Cannot analyze gamma wavelet at this time')


### OSC EVENT INFO DICTS !!
thetaOscEventInfo = {'chan': 8, 'minT': 2785.22321038684, 
					'maxT': 3347.9278996316607, 'alignoffset':-3086.95, 'left': 55704, 'right':66958,
					'w2': 3376}  # 

# betaOscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': 0, 'left': 0, 'right': 100, 'w2':100}


###############################
####### Evaluating Pops #######
###############################

## Improve these two sections..... 

## Evaluate pops / wavelets by frequency band ## --> IN PROGRESS 
evalWavelets_Band = 0

if evalWavelets_Band:
	basedPkl = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/wavelets/sim/spont/' 			# '/figs/wavelets/'
	dlmsPklFile = 'v34_batch57_3_4_data_timeRange_0_6_dlms.pkl'
	dfPklFile = 'v34_batch57_3_4_data_timeRange_0_6_df.pkl'   ### AUTOMATE / CONDENSE THIS SOMEHOW... 
	# dlmsData, dfData = evalWaveletsByBand(based=basedPkl, dlmsPklFile=dlmsPklFile, dfPklFile=dfPklFile)
	dfData = evalWaveletsByBand(based=basedPkl, dfPklFile=dfPklFile)  # THIS FUNCTION IS IN PROGRESS 


## Automated pop selection based on max (peak or avg) CSD values ## 
evalPopsBool = 0
if evalPopsBool:

	# Put any test info here # 
	beta2_getWaveletInfo = {'dataFile': 'v34_batch67_CINECA_0_0_data.pkl', 'channel': 19, 'timeRange': [2149.6607483037415, 2332.7116635583175]}
	based = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/'
	timeRange = beta2_getWaveletInfo['timeRange']
	dataFile = based + beta2_getWaveletInfo['dataFile']
	waveletElectrode = beta2_getWaveletInfo['channel']

	beta2_oscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': 0, 'left': 0, 'right': 100, 'w2': 100}


	print('timeRange: ' + str(timeRange))
	print('dataFile: ' + str(dataFile))
	print('channel: ' + str(waveletElectrode))

	# Get data frames for LFP and CSD data
	### dfPeak_LFP, dfAvg_LFP = getDataFrames(dataFile=dataFile, timeRange=timeRange)			# dfPeak, dfAvg, peakValues, avgValues, lfpPopData = getDataFrames(dataFile=dataFile, timeRange=timeRange, verbose=1)
	#dfPeak_CSD, dfAvg_CSD = getCSDDataFrames(dataFile=dataFile, timeRange=timeRange, oscEventInfo = thetaOscEventInfo)
	dfPeak_CSD, dfAvg_CSD = getCSDDataFrames(dataFile=dataFile, timeRange=timeRange, oscEventInfo = beta2_oscEventInfo)


	# Get the pops with the max contributions 
	maxPopsValues_peakCSD = evalPops(dataFrame=dfPeak_CSD, electrode=waveletElectrode)
	maxPopsValues_avgCSD = evalPops(dataFrame=dfAvg_CSD, electrode=waveletElectrode)

	# maxPopsValues_avgCSD['elec']


########################################################################
####### LOOKING AT ALL OSC EVENTS BY BAND -- FOR THESIS PROPOSAL #######
########################################################################
popsByBand = 0

if popsByBand:
	#####################
	#### LOCAL PATHS #### 
	#####################
	## based local ---> ##
	# based = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/'
	## waveletPathSim local ---> ##
	# waveletPathSim = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/wavelets/sim/spont/'

	######################
	### NEUROSIM PATHS ###
	######################
	### BASED NEUROSIM: 
	based = '/u/ericag/a1dat/spontData_A1/data_pklFiles/'
	### waveletPathSim NEUROSIM: 
	waveletPathSim = '/u/ericag/a1dat/spontData_A1/oscEvents/'


	### NORM EXPERIMENT ###
	norm=True
	#######################

	#######################
	pops = AllCortPops #ICortPops
	#######################


	layers = getSimLayers()

	frequencyBands = ['beta'] #['delta', 'theta', 'alpha', 'beta', 'gamma', 'hgamma']

	simSubjects = ['v34_batch67_CINECA_0_0_data', 'v34_batch67_CINECA_0_1_data', 'v34_batch67_CINECA_0_2_data', 'v34_batch67_CINECA_0_3_data', 'v34_batch67_CINECA_0_4_data',
				   'v34_batch67_CINECA_1_0_data', 'v34_batch67_CINECA_1_1_data', 'v34_batch67_CINECA_1_2_data', 'v34_batch67_CINECA_1_3_data', 'v34_batch67_CINECA_1_4_data',
				   'v34_batch67_CINECA_2_0_data', 'v34_batch67_CINECA_2_1_data', 'v34_batch67_CINECA_2_2_data', 'v34_batch67_CINECA_2_3_data', 'v34_batch67_CINECA_2_4_data',
				   'v34_batch67_CINECA_3_0_data', 'v34_batch67_CINECA_3_1_data', 'v34_batch67_CINECA_3_2_data', 'v34_batch67_CINECA_3_3_data', 'v34_batch67_CINECA_3_4_data',
				   'v34_batch67_CINECA_4_0_data', 'v34_batch67_CINECA_4_1_data', 'v34_batch67_CINECA_4_2_data', 'v34_batch67_CINECA_4_3_data', 'v34_batch67_CINECA_4_4_data']

	# waveletCounts_w_LFPrecording = getWaveletCounts(regions, frequencyBands, subjects_w_LFPrecording, waveletPathSim, sim=1, cutoffs=0)

	oscEventInfo = getOscEventInfo(subjects=simSubjects,frequencyBands=frequencyBands, waveletPath=waveletPathSim) # subjects = simSubjects

	for band in frequencyBands:
		for region in layers.keys():
			for subject in simSubjects:
				print('Analyzing subject ' + str(subject) + '; band: ' + str(band) + ' ; region: ' + str(region))

				dataFilePath = based + subject + '.pkl'
				# print('path to dataFile:' + str(dataFilePath))

				subject_all_oscEventInfo = oscEventInfo[band][region][subject]
				for eventIdx in subject_all_oscEventInfo.keys():
					individual_oscEventInfo = subject_all_oscEventInfo[eventIdx]

					timeRange = list(np.zeros(2))
					timeRange[0] = individual_oscEventInfo['minT']
					timeRange[1] = individual_oscEventInfo['maxT']

					print('Retrieving CSD data frame for subject: ' + str(subject))
					dfPeak_CSD, dfAvg_CSD = getCSDDataFrames(dataFile=dataFilePath, timeRange=timeRange, 
														oscEventInfo = individual_oscEventInfo, norm=norm, pops=pops)

					waveletElectrode = individual_oscEventInfo['chan']
					## max pops contributing to avg CSD 
					maxPopsValues_avgCSD = evalPops(dataFrame=dfAvg_CSD, electrode=waveletElectrode)
					## PRINTING / TESTING LINES 
					print('max pops for ' + subject + '\n'
							+ 'frequencyBand: ' + band + '\n'
							+ 'region: ' + region + '\n'
							+ 'idx: ' + str(eventIdx) + '\n'
							+ 'MAX POPS: ') 
					print(maxPopsValues_avgCSD)
					oscEventInfo[band][region][subject][eventIdx]['maxPops_avgCSD'] = maxPopsValues_avgCSD

					# max pops contributing to peak CSD
					maxPopsValues_peakCSD = evalPops(dataFrame=dfPeak_CSD, electrode=waveletElectrode)
					oscEventInfo[band][region][subject][eventIdx]['maxPops_peakCSD'] = maxPopsValues_peakCSD


	### Now print out a list with the top 3 (avg) pops for each osc event, organized by band and layer region 
	topPopsAvg = {}
	for band in frequencyBands:
		topPopsAvg[band] = {}
		for region in layers.keys():
			print('frequencyBand:' + band)
			print('layer region: ' + region)

			topPopsAvg[band][region] = {}

			for subject in simSubjects:
				eventIdxList = list(oscEventInfo[band][region][subject].keys())

				for eventIdx in eventIdxList:
					topPopsAvg[band][region][eventIdx] = []

					topPopKeys = oscEventInfo[band][region][subject][eventIdx]['maxPops_avgCSD']['elec'].keys()

					# topPops = []
					for key in topPopKeys:
						if not key == 'electrodeNumber':
							#topPops.append(key)
							topPopsAvg[band][region][eventIdx].append(key)

		print(topPopsAvg)

	### Now print out a list with the top 3 (peak) pops for each osc event, organized by band and layer region 
	topPopsPeak = {}
	for band in frequencyBands:
		topPopsPeak[band] = {}
		for region in layers.keys():
			print('frequencyBand:' + band)
			print('layer region: ' + region)

			topPopsPeak[band][region] = {}

			for subject in simSubjects:
				eventIdxList = list(oscEventInfo[band][region][subject].keys())

				for eventIdx in eventIdxList:
					topPopsPeak[band][region][eventIdx] = []

					topPopKeys = oscEventInfo[band][region][subject][eventIdx]['maxPops_peakCSD']['elec'].keys()

					# topPops = []
					for key in topPopKeys:
						if not key == 'electrodeNumber':
							#topPops.append(key)
							topPopsPeak[band][region][eventIdx].append(key)

		print(topPopsPeak)



		##### SAVING ######
		## SAVE topPopsAvg DICTIONARY!! 
		dictSavePath = waveletPathSim + 'oscEventDicts/' + 'AllCortPops/'
		os.chdir(dictSavePath)
		jsonFileAvg = 'topPopsAvg_' + str(band) + '_AllCortPops.json'
		print('Saving topPopsAvg json file!!')
		with open(jsonFileAvg, 'w') as fpAvg:
			json.dump(topPopsAvg, fpAvg)

		## SAVE topPopsPeak DICTIONARY!! 
		dictSavePath = waveletPathSim + 'oscEventDicts/' + 'AllCortPops/'
		os.chdir(dictSavePath)
		jsonFilePeak = 'topPopsPeak_' + str(band) + '_AllCortPops.json'
		print('Saving topPopsPeak json file!!')
		with open(jsonFilePeak, 'w') as fpPeak:
			json.dump(topPopsPeak, fpPeak)

		## SAVE oscEventInfo DICTIONARY!! 
		oscEventFile = 'oscEventInfo_' + str(band) + '_AllCortPops.pkl'
		print('Saving oscEventInfo pkl file!')
		with open(oscEventFile, 'wb') as f:
			pickle.dump(oscEventInfo, f, protocol=pickle.HIGHEST_PROTOCOL)



###################################
######## PLOTTING HEATMAPS ########
###################################


##### CSD HEATMAPS ######
plotCSDheatmaps = 1

if plotCSDheatmaps:   # define electrodes & figSize? 

	dataFile = '../data/v34_batch67/v34_batch67_CINECA_0_0_data.pkl'
	timeRange = [2149.66, 2332.71]
	betaOscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': -2270.0, 'left': 42993, 'right': 46654, 'w2':1098}


	## Get peak and avg dataframes 
	dfCSDPeak, dfCSDAvg = getCSDDataFrames(dataFile=dataFile, timeRange=timeRange, oscEventInfo=betaOscEventInfo) ## going to need oscEventInfo here 

	## Plot peak CSD heatmap 
	peakCSDPlot = plotDataFrames(dfCSDPeak, electrodes=None, pops=ECortPops, title='Peak CSD Values', 
					cbarLabel='CSD', figSize=(10,7), savePath=None, saveFig=False)

	## Plot avg CSD heatmap 
	avgCSDPlot = plotDataFrames(dfCSDAvg, electrodes=None, pops=ECortPops, title='Avg CSD Values', 
					cbarLabel='CSD', figSize=(10,7), savePath='../../A1/data/v34_batch67/', saveFig=True)
 


##### LFP HEATMAPS ######
plotLFPheatmaps = 1

if plotLFPheatmaps:
	dfLFPPeak, dfLFPAvg = getDataFrames(dataFile=dataFile, timeRange=timeRange)
	peakLFPPlot = plotDataFrames(dfLFPPeak, electrodes=None, pops=ECortPops, title='Peak LFP Values', 
					cbarLabel='LFP', figSize=(10,7), savePath=None, saveFig=False)
	avgLFPPlot = plotDataFrames(dfLFPAvg, electrodes=None, pops=ECortPops, title='Avg LFP Values', 
					cbarLabel='LFP', figSize=(10,7), savePath='../../A1/data/v34_batch67/', saveFig=False)




###################################
###### COMBINED LFP PLOTTING ######
###################################


### NEED TO DEFINE:
#### timeRange --> NOTE that this is contingent on getWaveletInfo lines in the specific-wavelet section


plotLFPCombinedData = 0


if plotLFPCombinedData:

	includePops = ['IT3'] #, 'IT5A', 'PT5B']  			## Automate this somehow?? 
	electrode = [10]									## Automate this somehow?? 


	for pop in includePops:

		print('Plotting LFP spectrogram and timeSeries for ' + pop + ' at electrode ' + str(electrode))

		## Get dictionaries with LFP data for spectrogram and timeSeries plotting  
		LFPSpectOutput = getLFPDataDict(dataFile, pop=pop, timeRange=timeRange, plotType=['spectrogram'], electrode=electrode) 
		LFPtimeSeriesOutput = getLFPDataDict(dataFile, pop=pop, timeRange=timeRange, plotType=['timeSeries'], electrode=electrode)

		plotCombinedLFP(timeRange=timeRange, pop=pop, colorDict=colorDict, plotTypes=['spectrogram','timeSeries'], 
			spectDict=LFPSpectOutput, timeSeriesDict=LFPtimeSeriesOutput, figSize=(10,7), colorMap='jet', 
			minFreq=15, maxFreq=None, vmaxContrast=None, electrode=electrode, savePath=None, saveFig=True)





#######################################
######## COMBINED CSD PLOTTING ########
#######################################


plotCSDCombinedData = 0



dataFile = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/v34_batch67_CINECA_0_0_data.pkl'
betaOscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': -2270.0, 'left': 42993, 'right': 46654, 'w2':1098}
oscEventInfo = betaOscEventInfo


if plotCSDCombinedData:

	# print('Plotting Combined CSD data')
	electrode=[19] #[8]
	includePops=['IT6']#['CT6', 'IT6', 'PT5B']	#['ITS4']	#[None, 'ITS4', 'ITP4', 'IT5A'] # ['IT3', 'ITS4', 'ITP4', 'IT5A', 'PT5B']

	minFreq = 1 			# 0.25 # 1 
	maxFreq = 50 #12 			# 110 # 40 # 25 
	stepFreq = 0.25 		# 1 # 0.25 

	peakFDict = {'CT6': 15.25, 'IT6': 15.5}

	for pop in includePops:

		peakFbyPop = peakFDict[pop]
		customPopTitle = 'Model ' + pop + ' (Beta)' #' (Beta, ' + str(peakFbyPop) + ' Hz)'

		print('Plotting CSD spectrogram and timeSeries for ' + pop + ' at electrode ' + str(electrode))

		## Get dictionaries with CSD data for spectrogram and timeSeries plotting 
		timeSeriesDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=oscEventInfo, pop=pop, minFreq=minFreq, maxFreq=maxFreq, stepFreq=stepFreq)
		spectDict = getCSDdata(dataFile=dataFile, outputType=['spectrogram'], oscEventInfo=oscEventInfo, pop=pop, minFreq=minFreq, maxFreq=maxFreq, stepFreq=stepFreq)


		plotCombinedCSD(timeSeriesDict=timeSeriesDict, spectDict=spectDict, colorDict=colorDict, customPopTitle=customPopTitle, pop=pop, electrode=electrode, 
			minFreq=1, maxFreq=maxFreq, vmaxContrast=None, colorMap='jet', figSize=(10,7), plotTypes=['timeSeries', 'spectrogram'], 
			hasBefore=1, hasAfter=1, saveFig=True) # colorDict=colorDictCustom 


##########################################
###### COMBINED SPIKE DATA PLOTTING ######
##########################################


plotCombinedSpikeData = 0

dataFile = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/v34_batch67_CINECA_0_0_data.pkl'
betaOscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': -2270.0, 'left': 42993, 'right': 46654, 'w2':1098}
oscEventInfo = betaOscEventInfo


if plotCombinedSpikeData:

	includePops=['CT6', 'IT6', 'PT5B', 'CT5A'] #['ITS4', 'ITP4', 'IT5A']		# ['IT3', 'ITS4', 'IT5A']  ## AUTOMATE THIS SOMEHOW? i.e. with eval pops fx's? 

	for pop in includePops:
		print('Plotting spike data for ' + pop)

		## Get dictionaries with spiking data for spectrogram and histogram plotting 
		spikeSpectDict = getSpikeData(dataFile, graphType='spect', pop=pop, oscEventInfo=thetaOscEventInfo)		#	timeRange=timeRange)
		histDict = getSpikeData(dataFile, graphType='hist', pop=pop, oscEventInfo=thetaOscEventInfo)			#	timeRange=timeRange)

		## Then call plotting function
		plotCombinedSpike(spectDict=spikeSpectDict, histDict=histDict, colorDict=colorDict, plotTypes=['spectrogram', 'histogram'],
		hasBefore=1, hasAfter=1, pop=pop, figSize=(10,7), colorMap='jet', vmaxContrast=2, maxFreq=None, saveFig=1) 	# colorDictCustom	# timeRange=timeRange, 

# TO DO: Smooth or mess with bin size to smooth out spectrogram for spiking data




################
##### PSD ######
################


### LFP ### ---> Evidently this is code for PSD of a summed LFP signal ; adapt for other contexts? 
summedLFP = 0 	 		# COMBINE TOP 3 LFP SIGNALS  

if summedLFP: 
	includePops = ['IT3', 'IT5A', 'PT5B']
	popElecDict = {'IT3': 1, 'IT5A': 10, 'PT5B': 11}
	lfpDataTEST_fullElecs = getSumLFP(dataFile=dataFile, pops=includePops, elecs=False, timeRange=timeRange, showFig=True)
	# lfpDataTEST = getSumLFP(dataFile=dataFile, popElecDict=popElecDict, timeRange=timeRange, showFig=False)


lfpPSD = 0 ## fix bug !! see nb for error 
if lfpPSD: 
	psdData = getPSDdata(dataFile=dataFile, inputData = lfpDataTEST_fullElecs['sum'])	# inputData = lfpDataTEST['sum'])
	plotPSD(psdData)



### CSD ###
csdPSD = 0

dataFile = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/v34_batch67_CINECA_0_0_data.pkl'
betaOscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': -2270.0, 'left': 42993, 'right': 46654, 'w2':1098}
oscEventInfo = betaOscEventInfo

if csdPSD:
	includePops=['CT6', 'IT6']		#['ITS4']			#['ITS4', 'ITP4', 'IT5A']
	for pop in includePops:
		csdDataDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=betaOscEventInfo, pop=pop)  # thetaOscEventInfo
		csdData = csdDataDict['csdDuring'] 
		print('FOR POPULATION: ' + pop + ' , psdData')
		psdData = getPSDdata(dataFile=dataFile, inputData=csdData, inputDataType='timeSeries', minFreq=1, maxFreq=50, stepFreq=0.1)
		# plotPSD(psdData)
		### Unclear on significance of these two sections -- figure that out / look back in nb 
		spectDict = getCSDdata(dataFile=dataFile, outputType=['spectrogram'], oscEventInfo=betaOscEventInfo, pop=pop, maxFreq=40)  # thetaOscEventInfo,
		print('FOR POPULATION: ' + pop + ' , psdDataSpect')
		psdDataSpect = getPSDdata(dataFile=dataFile, inputData=spectDict, inputDataType='spectrogram', duringOsc=1, minFreq=1, maxFreq=40, stepFreq=1)
		#plotPSD(psdDataSpect)



### SUMMED CSD (FROM MULTIPLE POPS) ### 
csdPSD_multiple = 0

if csdPSD_multiple:
	includePops=['ITS4'] 				# ['ITS4', ITP4', 'IT5A'] # ECortPops 
	csdPopData = {}
	for pop in includePops:
		csdDataDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=thetaOscEventInfo, pop=pop)
		csdData = csdDataDict['csdDuring'] 
		csdPopData[pop] = csdData

	csdSummedData =  np.zeros(shape=csdPopData[includePops[0]].shape)
	for pop in includePops:
		csdSummedData += csdPopData[pop]

	psdSummedData = getPSDdata(dataFile=dataFile, inputData=csdSummedData, inputDataType='timeSeries', minFreq=1, maxFreq=100, stepFreq=0.25)
	plotPSD(psdSummedData)


### ENTIRE CSD (DURING OSC EVENT, AT SPECIFIED CHANNEL) ### 
csdPSD_wholeCSD = 0

if csdPSD_wholeCSD:
	maxFreq = 110 
	pop = 'IT5A' #None

	csdDataDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=thetaOscEventInfo, pop=pop)
	csdData = csdDataDict['csdDuring'] 

	psdData = getPSDdata(dataFile=dataFile, inputData=csdData, inputDataType='timeSeries', minFreq=1, maxFreq=maxFreq, stepFreq=0.25)
	plotPSD(psdData)



### LOOK AT MAX POWER FOR EACH POP -- TO DO: MAKE THIS INTO A FUNCTION !!!! ### 
PSDbyPop = 0

dataFile = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/v34_batch67_CINECA_0_0_data.pkl'
betaOscEventInfo = {'chan': 19, 'minT': 2149.6607483037415, 'maxT': 2332.7116635583175, 'alignoffset': -2270.0, 'left': 42993, 'right': 46654, 'w2':1098}
oscEventInfo = betaOscEventInfo


if PSDbyPop:
	includePops = ['CT6', 'IT6']		#ECortPops    					# ['IT3'] 	# thalPops 		# AllCortPops
	maxPowerByPop = {}

	for pop in includePops:
		# Get the max power frequency for each pop and put it in a dict
		csdDataDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=betaOscEventInfo, pop=pop)  # thetaOscEventInfo
		csdDataPop = csdDataDict['csdDuring']

		psdDataPop = getPSDdata(dataFile=dataFile, inputData=csdDataPop, inputDataType='timeSeries', minFreq=1, maxFreq=100, stepFreq=0.25)
		maxPowerByPop[pop] = psdDataPop['maxPowerFrequency']

	{k: v for k, v in sorted(maxPowerByPop.items(), key=lambda item: item[1])}




############################
###### load.py checks ######
############################

#### plot blob ####
dataFile = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/v34_batch67_CINECA_0_0_data.pkl'

plotBlobs = 0
if plotBlobs:
	sim.load(dataFile, instantiate=False)

	# Calculate sampling rate from simulation timestep 
	dt = sim.cfg.recordStep  	# or should I divide by 1000.0 up here, and then just do 1.0/dt below for sampr? 
	sampr = 1.0/(dt/1000.0) 	# divide by 1000.0 to turn denominator from units of ms to s

	spacing_um = 100
	norm=True
	lfpData = np.array(sim.allSimData['LFP'])

	inputData = getCSDa1dat(lfps=lfpData,sampr=sampr, spacing_um=spacing_um, minf=0.05, maxf=300,norm=norm,vaknin=True)
	inputData_chan19 = inputData[19, :]
	peakFData = getPeakF(dataFile=dataFile,inputData=inputData_chan19, csdAllChans=inputData, chan=19, plotTest=True, plotNorm=0)




#### peakF calculations ####
peakF = 0
if peakF:
	maxFreq = 20 #110 #10
	plotNorm = 0
	chan=8
	pop = None #'ITP4'

	timeSeriesDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=thetaOscEventInfo, pop=pop, maxFreq=maxFreq)

	## CSD DATA
	csdData = timeSeriesDict['csdData']   					## All chans, all timepoints 
	fullTimeRange = [0,(sim.cfg.duration/1000.0)] 
	dt = sim.cfg.recordStep / 1000.0  						# thus, dt also converted to seconds (from ms)
	tt = np.arange(fullTimeRange[0],fullTimeRange[1],dt) 	# tt = timeSeriesDict['tt']

	## timeRange so it's like load.py 	# put timeRange in seconds
	timeRange = [0,6] # [0,11.5]

	## Segment csdData and tt by timeRange
	dat = csdData[:, int(timeRange[0]/dt):int(timeRange[1]/dt)]
	datChan = dat[chan,:]
	tt = tt[int(timeRange[0]/dt):int(timeRange[1]/dt)]

	peakFData = getPeakF(dataFile=dataFile, inputData=datChan, csdAllChans=dat, timeData=tt, chan=chan, freqmax=maxFreq, plotTest=True, plotNorm=plotNorm)
	# peakFData = getPeakF(dataFile=dataFile, inputData=csdData_theta, csdAllChans=csdData_theta_allChans, timeData=tt, freqmax=maxFreq, plotTest=False, plotNorm=plotNorm)
	# peakFData = getPeakF(dataFile=dataFile, inputData=csdDuring, csdAllChans=csdDataAllChans, timeData=tt_During, freqmax=maxFreq, plotTest=True, plotNorm=plotNorm)
	# peakFData = getPeakF(dataFile=dataFile, inputData=csdOscChan_plusTimeBuffer, csdAllChans=csdAllChans_plusTimeBuffer, timeData=tt_plusTimeBuffer, 
				# freqmax=maxFreq, plotTest=True, plotNorm=plotNorm)

	imgpk = peakFData['imgpk']
	imgpk_nonNorm = peakFData['imgpk_nonNorm']
	lms = peakFData['lms']
	lsidx = peakFData['lsidx']
	leidx = peakFData['leidx']
	lmsnorm = peakFData['lmsnorm']
	lnoise = peakFData['lnoise']
	## lblob
	lblob_norm = peakFData['lblob_norm']
	lblob_nonNorm = peakFData['lblob_nonNorm']
	lblob = lblob_nonNorm
	## llevent
	llevent = peakFData['llevent']
	llevent = llevent[0]
	llevent_norm = peakFData['llevent_norm']
	llevent_norm = llevent_norm[0]

	peaks = np.where(imgpk==True)
	peaks_nonNorm = np.where(imgpk_nonNorm==True)


	### LOOK AT CANDIDATES
	print('Looking at llevent')
	for i in range(len(llevent)):
		if llevent[i].peakF > 4 and llevent[i].peakF < 7:
			print('index: ' + str(i) + ' --> peakF: ' + str(llevent[i].peakF) + ', peakT: ' + str(llevent[i].peakT))


	print('Looking at llevent_norm')
	for i in range(len(llevent_norm)):
		if llevent_norm[i].peakF > 4 and llevent_norm[i].peakF < 7:
			print('index: ' + str(i) + ' --> peakF: ' + str(llevent_norm[i].peakF) + ', peakT: ' + str(llevent_norm[i].peakT))


	print('Looking at lblob')
	for i in range(len(lblob)):
		if lblob[i].peakF > 4 and lblob[i].peakF < 7:
			print('index: ' + str(i) + ' --> peakF: ' + str(lblob[i].peakF) + ', peakT: ' + str(lblob[i].peakT))


	print('Looking at lblob_norm')
	for i in range(len(lblob_norm)):
		if lblob_norm[i].peakF > 4 and lblob_norm[i].peakF < 7:
			print('index: ' + str(i) + ' --> peakF: ' + str(lblob_norm[i].peakF) + ', peakT: ' + str(lblob_norm[i].peakT))


	# x = zip(np.arange(len(lms)),lsidx,lms,lmsnorm,lnoise)




#### getIEIstatsbyBand Testing ####
getIEIstatsbyBandTEST=0
if getIEIstatsbyBandTEST:
	maxFreq = 110 			# 10
	chan=8
	pop = None 				# 'ITP4'	

	timeSeriesDict = getCSDdata(dataFile=dataFile, outputType=['timeSeries'], oscEventInfo=thetaOscEventInfo, pop=pop, maxFreq=maxFreq)

	## CSD DATA
	csdData = timeSeriesDict['csdData']   					## All chans, all timepoints 
	fullTimeRange = [0,(sim.cfg.duration/1000.0)] 
	dt = sim.cfg.recordStep / 1000.0  						# thus, dt also converted to seconds (from ms)
	sampr = 1.0 / dt
	tt = np.arange(fullTimeRange[0],fullTimeRange[1],dt) 	# tt = timeSeriesDict['tt']

	## timeRange so it's like load.py
	timeRange = [0,6]					# in seconds 

	## Segment csdData and tt by timeRange
	dat = csdData[:, int(timeRange[0]/dt):int(timeRange[1]/dt)]
	datChan = dat[chan,:]
	tt = tt[int(timeRange[0]/dt):int(timeRange[1]/dt)]



	noiseampCSD = 200.0 / 10.0 # amplitude cutoff for CSD noise; was 200 before units fix
	noiseamp=noiseampCSD 
	winsz = 10
	medthresh = 4.0
	lchan = [chan]
	MUA = None


	dout = getIEIstatsbyBand2(inputData=datChan,winsz=winsz,sampr=sampr,freqmin=0.25,freqmax=maxFreq,freqstep=0.25,
		medthresh=medthresh,lchan=lchan,MUA=MUA,overlapth=0.5,getphase=True,savespec=True,
		threshfctr=2.0,useloglfreq=False,mspecwidth=7.0,noiseamp=noiseampCSD,endfctr=0.5,
		normop=mednorm)



####################################
##### COMPARING CINECA AND GCP ##### 
####################################

compareHPC = 0

dataFileCineca = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/v34_batch67_CINECA/data_pklFiles/v34_batch67_CINECA_0_0_data.pkl'
dataFileGCP = '/Users/ericagriffith/Desktop/NEUROSIM/A1/data/simDataFiles/spont/A1_v34_batch67_v34_batch67_0_0_data.pkl' #v34_batch57/v34_batch57_0_0_data.pkl'

if compareHPC:
	sim.load(dataFileCineca, instantiate=False)
	cineca_allData = sim.allSimData
	cineca_lfp = cineca_allData['LFP']			#sim.allSimData['LFP']
	## CHECK LFP POPS 
	cineca_lfpPops = cineca_allData['LFPPops']
	## CHECK SPIKE TIME
	cineca_spkt = cineca_allData['spkt']


	sim.load(dataFileGCP, instantiate=False)
	gcp_allData = sim.allSimData
	gcp_lfp = gcp_allData['LFP']				#sim.allSimData['LFP']
	## CHECK LFP POPS 
	#### gcp_lfpPops = gcp_allData['LFPPops'] <-- WILL NOT EXIST!!! DO THIS ONLY WITH THE batch65 batch67 0_0 0_1 files 
	## CHECK SPIKE TIME
	gcp_spkt = gcp_allData['spkt']

	## TRY PLOTTING SOMETHING EQUIVALENT


	## Differences between spike times:
	diff_inds = []
	diff_values = []
	for i in range(len(gcp_spkt)):
		if gcp_spkt[i]!= cineca_spkt[i]:
			diff_inds.append(i)
			diff_values.append([gcp_spkt[i], cineca_spkt[i]])














