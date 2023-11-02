"""
batchAnalysis.py 

Code to anlayse batch sim results

Contributors: salvadordura@gmail.com
"""

#import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers

from batchAnalysisFilter import *
from batchAnalysisPlotCombined import *


# Main code
if __name__ == '__main__': 
    dataFolder = '../data/'
    batchLabel = 'v34_batch7'  # 'v50_batch1' #
    #batchLabels = ['v103_batch3/gen_%d' % (i) for i in range(68)]
    loadAll = 1

    # ---------------------------------------------
    # Filtering wrapper funcs
    # ---------------------------------------------

    #df = filterDepolBlock(dataFolder, batchLabel, loadAll, gids=[])
    
    # df = applyFilterRates(dataFolder, batchLabel, loadAll, skipDepol=0)  
    
    # filterStimRates(dataFolder, batchLabel, load=loadAll)

    # var = [('simData','popRates')]
    # params, data = utils.readBatchData(dataFolder, batchLabel, loadAll=loadAll, saveAll=1-loadAll, vars=var, maxCombs=None)

    # Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'CT5B' , 'PT5B', 'IT6', 'CT6']  # all layers

    # Ipops = ['NGF1',                            # L1
    #         'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
    #         'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
    #         'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
    #         'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
    #         'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
    #         'PV6', 'SOM6', 'VIP6', 'NGF6']      # L6 

    # loadAll = 1
    # df = addFitness(dataFolder, batchLabel, loadAll, tranges=[[1500, 1750], [1750,2000], [2000,2250], [2250,2500]], Epops=Epops, Ipops=Ipops)
    


    # ---------------------------------------------
    # Combined sim (batch) funcs
    # ---------------------------------------------

    #allpops = ['NGF1', 'IT2', 'PV2', 'SOM2', 'VIP2', 'NGF2', 'IT3', 'SOM3', 'PV3', 'VIP3', 'NGF3', 'ITP4', 'ITS4', 'PV4', 'SOM4', 'VIP4', 'NGF4', 'IT5A', 'CT5A', 'PV5A', 'SOM5A', 'VIP5A', 'NGF5A', 'IT5B', 'PT5B', 'CT5B', 'PV5B', 'SOM5B', 'VIP5B', 'NGF5B', 'IT6', 'CT6', 'PV6', 'SOM6', 'VIP6', 'NGF6', 'TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI']
    
    fIAnalysis(dataFolder, batchLabel, loadAll)
    
    #df = popRateAnalysis(dataFolder, batchLabel, loadAll, pars=['weightBkgE', 'weightBkgI','rateBkg_exc', 'rateBkg_inh'], vals=['ITS4'], groupStat='first', plotLine=False) 

    # df = ihEPSPAnalysis(dataFolder, batchLabel, loadAll, pars=['groupWeight','ihGbar'], vals=['Vpeak_PTih'], zdcomp=0, plotLine=1)#, \
    # query = 'epas == 1.0 and groupWeight > 0.0003')# and axonNa==7 and gpas==0.65') #'ihLkcBasal == 0.01 and excTau2Factor==1.0') #, 'excTau2Factor', 'ihLkcBasal

    # df = ihEPSPAnalysis(dataFolder, batchLabel, loadAll, pars=['groupWeight', 'ihGbar'], vals=['Vpeak_PTih'], plotLine=True, \
    #    query = 'epas==1.0 and gpas==0.8') #, 'excTau2Factor', 'ihLkcBasal

    # extractRatePSD(dataFolder, batchLabel)

    #extractRates(dataFolder, batchLabel)

    # plotSeedsRates(dataFolder, batchLabel)

    #df = plotSimultLongRates(dataFolder, batchLabel)

    #fig, allSignal = freqAnalysis(dataFolder, batchLabel)
        

    # ---------------------------------------------
    # Comparing funcs
    # ---------------------------------------------

    #utils.compare(dataFolder+'v52_batch11/v52_batch11_0_0.json', dataFolder+'v52_batch12/v52_batch12_0_1_0.json', source_key='simConfig', target_key='simConfig')

    #utils.compare(dataFolder+'../sim/net1.json', dataFolder+'../sim/net2.json')#, source_key='net', target_key='net')

    



