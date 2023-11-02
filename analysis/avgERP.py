import matplotlib
from pylab import *

import numpy,scipy
from csd import *
from erp import *
from nhpdat import *

import os 
from netpyne import sim 



##############################


def drawDataAvgERP(usemodel, dataFile, spacing_um=100.0, swindowms=0, ewindowms=50, lfpSmooth=1, ERPdata='LFP', region='gran', chan=None, saveFigs=1, saveDir=None, dataFilesInfoDict=None): 
	### usemodel: bool 							--> 0: NHP, 1: SIM
	### dataFile: str 								--> FULL PATH TO DATA FILE 
	### spacing_um: spacing btwn electrodes 		--> DEFAULT: 100 uM
	### ewindowms: float / int 						--> DEFAULT: 50 ms
	### lfpSmooth: bool 							--> DEFAULT: 1 
	### ERPdata: str 								--> DEFAULT: 'LFP' (alt option: 'CSD')
	### region: str 								--> DEFAULT: 'gran' (alt options: 'supra', 'infra')
	### dataFileDict: dict


	if usemodel:
		# Get LFP data
		sim.load(dataFile, instantiate=False)
		LFP = np.array(sim.allSimData['LFP'])

		## Get sampling rate and spacing_um 
		dt = sim.cfg.recordStep   # dt is in ms by default 
		sampr = 1.0/(dt/1000.0)   # divide by 1000 to get sampr in Hz (1/s)


		# Get CSD
		CSD = getCSD(LFP,sampr,spacing_um=spacing_um, minf=1, maxf=110, vaknin=True)

		# Get trigger (stimulus input) times 
		trigtimes_ms = sim.cfg.ICThalInput['startTime']
		if type(trigtimes_ms) is not list:
			trigtimes_ms = [trigtimes_ms]
		trigtimes = []
		for trigtime_ms in trigtimes_ms:
			trigtime = ms2index(trigtime_ms, sampr)
			trigtimes.append(trigtime)


		if type(trigtimes) == list:
			tts = trigtimes
		else:
			tts = [trigtimes]

		# layer / cortical region label information 
		values = [0, 100, 160, 950, 950, 1250, 1250, 1334, 1550, 2000]  

		# set epoch params   ### <-- This should be the length of the stimulus (approx)!
		if swindowms is None:
			swindowms = 0
		if ewindowms is None:
			ewindowms = 50 

		### LFP SMOOTHING
		smooth = 80
		if lfpSmooth:
			LFP = getbandpass(LFP,sampr,1,smooth)
		else:
			LFP = LFP.T

		### CSD SMOOTHING
		smooth = 30
		if lfpSmooth:
			CSD = getbandpass(CSD.T,sampr,1,smooth)

		## Get channel(s) to plot
		chan = region

		if ERPdata == 'LFP':
			drawAvgERP(LFP, sampr, tts, swindowms, ewindowms, whichchan=chan)
			plt.ylabel('LFP (mV)')
		elif ERPdata == 'CSD':
			drawAvgERP(CSD, sampr, tts, swindowms, ewindowms, whichchan=chan, lw=2)
			plt.ylabel('CSD (mV/mm^2)')


		# Title 
		filenameStr = dataFile.split('/')[-1]
		if 'BBN' in dataFile: 
			stimType = 'BBN'
		if chan is None: 
			figTitle = ERPdata + ' avgERP for MODEL A1 with ' + stimType + ' input (n=%d)'%len(tts) + '\n' + filenameStr + ', smooth_' + str(smooth)
		else:
			figTitle = ERPdata + ' avgERP for MODEL A1 with ' + stimType + ' input (n=%d)'%len(tts) + '\n' + filenameStr + ', ' + str(region) + ', smooth_' + str(smooth)#chan ' + str(chan)
		plt.suptitle(figTitle)


		# Axes labels 
		plt.xlabel('Time (ms)')

		# Fig Saving 
		if saveFigs:
			print('Saving Figure')
			if chan is None:
				saveFilename = filenameStr.split('.pkl')[0] + '_avgERP_' + ERPdata + '.png'
			else:
				saveFilename = filenameStr.split('.pkl')[0] + '_avgERP_' + ERPdata + '_chan' + str(chan) + '_'  + str(swindowms) + '_' + str(ewindowms) + 'ms' + '_smooth_' + str(smooth)+ '.png'
			print('filenameStr: ' + filenameStr)
			print('saveFilename: ' + saveFilename)
			if saveDir is not None:
				saveFilePath = saveDir + saveFilename
				plt.savefig(saveFilePath)
				plt.close()
			else:
				print('need save dir!')
				plt.show()
		else:
			plt.show()

	else:
		fn = dataFile
		samprds = getdownsampr(fn)
		divby = getorigsampr(fn) / samprds
		trigtimes = [int(round(x)) for x in np.array(getTriggerTimes(fn)) / divby]   # NOTE: These are indices, not times
		trigIDs = getTriggerIDs(fn)


		sampr,LFP,dt,tt,CSD,MUA = loadfile(fn, samprds)
		stimFile = fn.split('/')[-1]
		
		chan = region


		### Calculate actual stimulus times ####
		stimTimes = []
		for idx in trigtimes:
			stimTimes.append(tt[idx])

		# set epoch params
		if swindowms is None:
			swindowms = 0
		if ewindowms is None:
			ewindowms = 50
		windowms = ewindowms - swindowms


		tts = trigtimes

		### LFP SMOOTHING 
		if lfpSmooth:
			LFP = getbandpass(LFP,sampr,1,80)#110)
		else:
			LFP = LFP.T

		### CSD SMOOTHING
		smooth = 30
		if lfpSmooth:
			CSD = getbandpass(CSD.T,sampr,1,smooth)
			#CSD = CSD.T

		plt.figure()

		if ERPdata == 'LFP':
			drawAvgERP(LFP, sampr, trigtimes, swindowms, ewindowms, whichchan=chan)
			plt.ylabel('LFP (mV)')
		elif ERPdata == 'CSD':
			drawAvgERP(CSD, sampr, trigtimes, swindowms, ewindowms, whichchan=chan, lw=2)
			plt.ylabel('CSD (mV/mm^2)')


		stimType = dataFilesInfoDict[stimFile]['stimType']
		figTitle = ERPdata + ' avgERP for NHP A1 with ' + stimType + ' input (n=%d)'%len(tts) + '\n' + stimFile + ', chan ' + str(chan) 
		plt.suptitle(figTitle)

		# Axes labels 
		plt.xlabel('Time (ms)')

		# Fig Saving 
		if saveFigs:
			print('Saving Figure')
			saveFilename = stimFile.split('.mat')[0]  + '_avgERP_' + ERPdata + '_chan ' + str(chan)+ '_'  + str(swindowms) + '_' + str(ewindowms) + 'ms' + '.png'
			if saveDir is not None:
				saveFilePath = saveDir + saveFilename
				plt.savefig(saveFilePath)
			else:
				print('Need saveDir!')
				plt.show()
		else:
			plt.show()







