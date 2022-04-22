'''
Macaque A1 model 
Local connectivity preprocessing script

- Loads experimental data from file
- Preprocesses so ready to use by NetPyNE high-level specs

'''
import numpy as np
import json
import csv

# ----------------------------------------------------------------------------------------------------------------
# Func to load data from published studies
# ----------------------------------------------------------------------------------------------------------------
def loadData():
    data = {'test': {'prob': 0.2, 'weight': 0.5}}  # oversimplified example -- will have data from different papers for differnet projections 
    
    # ----------------------------------------------------------------------------------------------------------------
    # load and pre-process Allen mouse V1 data (Billeh et al, 2019; https://www.dropbox.com/sh/xb7xasih3d8027u/AAAbKXe0Zmk86o3_y1iPVPCLa?dl=0)
    data['Allen_V1'] = {}
    
    ## load conn probs
    with open('../data/conn/Allen_V1_conn_probs.json', 'r') as f:
        data['Allen_V1']['connProb'] = json.load(f)

    ## func to calculate probConn at distance = 0um based on probConn at distance = 75 and lambda (length constant)
    ## adapted from Allen V1 code in Allen_V1_connect_cells.py

    for proj in data['Allen_V1']['connProb'].keys():

        # A_literature is different for every source-target pair and was estimated from the literature.
        A_literature = data['Allen_V1']['connProb'][proj]['A_literature']

        # R0 read from the dictionary, but setting it now at 75um for all cases but this allows us to change it
        R0 = data['Allen_V1']['connProb'][proj]['R0']

        # Sigma is measure from the literature or internally at the Allen Institute
        sigma = data['Allen_V1']['connProb'][proj]['sigma']

        # Gaussian equation was intergrated to and solved to calculate new A0 (amplitude at distance=0)
        A0 = A_literature / ((sigma / R0)** 2 * (1 - np.exp(-(R0 / sigma)** 2)))
        if A0 > 1.0: A0 = 1.0  # make max prob = 1.0
        data['Allen_V1']['connProb'][proj]['A0'] = A0


    ## load conn weights (soma PSP amp in mV)
    data['Allen_V1']['connWeight'] = {}
    with open('../data/conn/Allen_V1_conn_weights.csv', 'r') as f:
        csv_reader = csv.reader(f)
        for irow,row in enumerate(csv_reader):
            if irow == 0:
                headings = row
            else:
                for ih, h in enumerate(headings[1:]):
                    try:
                        data['Allen_V1']['connWeight'][row[0] + '-' + h] = float(row[ih + 1])
                    except:
                        data['Allen_V1']['connWeight'][row[0] + '-' + h] = 0.0
                
    
    # set correspondence between A1 pops and Allen V1 pops 
    data['Allen_V1']['pops'] = {
        'NGF1': 'i1H',                                                                              # L1
        'IT2': 'e2',                'PV2': 'i2P',   'SOM2': 'i2S',  'VIP2': 'i2P',  'NGF2': 'i2H', # L2
        'IT3': 'e2',                'PV3': 'i2P',   'SOM3': 'i2S',  'VIP3': 'i2P',  'NGF3': 'i2H',  # L3
        'ITP4': 'e4', 'ITS4': 'e4', 'PV4': 'i4P',   'SOM4': 'i4S',  'VIP4': 'i4P',  'NGF4': 'i4H',  # L4
        'IT5A': 'e5', 'CT5A': 'e5', 'PV5A': 'i5P',  'SOM5A': 'i5S', 'VIP5A': 'i5P', 'NGF5A': 'i5H', # L5A
        'IT5B': 'e5', 'PT5B': 'e5', 'CT5B': 'e5',   'PV5B': 'i5P',  'SOM5B': 'i5S', 'VIP5B': 'i5P', 'NGF5B': 'i5H', # L5B
        'IT6': 'e6',  'CT6': 'e6',  'PV6': 'i6P',   'SOM6': 'i6S',  'VIP6': 'i6P',  'NGF6': 'i6H'}  # L6


    # ----------------------------------------------------------------------------------------------------------------
    # load and pre-process BBP mouse S1 data (Markram et al, 2015; https://bbp.epfl.ch/nmc-portal/downloads
    data['BBP_S1'] = {}
    # field to use -> data['BBP_S1']['connProb'][projection]['connection_probability']
    # project format = '[pre pop]:[post pop]' e.g. 'L5_TTPC1:L23_SBC'
    with open('../data/conn/BBP_S1_pathways_anatomy_factsheets_simplified.json', 'r') as f:
        data['BBP_S1']['connProb'] = json.load(f) 

    # field to use -> data['BBP_S1']['connWeight'][projection]['epsp_mean']
    with open('../data/conn/BBP_S1_pathways_physiology_factsheets_simplified.json', 'r') as f:
        data['BBP_S1']['connWeight'] = json.load(f)
    
    # calculate A0 so can use in combination with Allen 
    ## note BBP connProb represent "avg probability within 100um"; can approximate as prob at 75um used in Allen (R0)
    ## need to calculate corresponding A0 (max prob) based on R0 (prob at 75um) for BBP
    for proj in data['BBP_S1']['connProb']:
        A_literature = data['BBP_S1']['connProb'][proj]['connection_probability'] / 100.
        sigma = 75
        A0 = A_literature / ((sigma / R0)** 2 * (1 - np.exp(-(R0 / sigma)** 2)))
        if A0 > 1.0: A0 = 1.0  # make max prob = 1.0
        data['BBP_S1']['connProb'][proj]['A0'] = A0

    # set correspondence between A1 pops and Allen V1 pops 
    data['BBP_S1']['pops'] = {
        'NGF1': 'L1_NGC-DA',                                                                                                                             # L1
        'IT2':  'L23_PC',                                             'PV2':  'L23_LBC',   'SOM2': 'L23_MC',  'VIP2': 'L23_BP', 'NGF2':  'L23_NGC-DA', # L2
        'IT3':  'L23_PC',                                             'PV3':  'L23_LBC',   'SOM3': 'L23_MC',  'VIP3': 'L23_BP', 'NGF3':  'L23_NGC-DA', # L3
        'ITP4': 'L4_PC',     'ITS4': 'L4_SS',                         'PV4':  'L4_LBC',    'SOM4': 'L4_MC',   'VIP4': 'L4_BP',  'NGF4':  'L4_NGC-DA',  # L4
        'IT5A': 'L5_UTPC',   'CT5A': 'L6_TPC_L4',                     'PV5A': 'L5_LBC',   'SOM5A': 'L5_MC',  'VIP5A': 'L5_BP',  'NGF5A': 'L5_NGC-DA',  # L5A
        'IT5B': 'L5_UTPC',   'CT5B': 'L6_TPC_L4', 'PT5B': 'L5_TTPC2', 'PV5B': 'L5_LBC',   'SOM5B': 'L5_MC',  'VIP5B': 'L5_BP',  'NGF5B': 'L5_NGC-DA',  # L5B
        'IT6':  'L6_TPC_L1', 'CT6':  'L6_TPC_L4',                     'PV6':  'L6_LBC',    'SOM6': 'L6_MC',   'VIP6': 'L6_BP',  'NGF6':  'L6_NGC-DA'}  # L6


    # ------------------------------------------------------------------------------------------------------------------
    # load Thal -> A1 (Ji et al 2016)
    # mouse A1 MGB -> layers/cell types
    # no distinction core vs matrix (but most input from MGBv=core; same as allen)
    # adjusted amplitude pA ~= strength = prob * weight
    # % innervated ~= proxy for probability; many presyn axons innervated

    data['TC_Ji2016'] = {'innvervated': {}, 'amplitude': {}}

    data['TC_Ji2016']['innervated'] = {
        'L1': 12. / 22,
        'L23_Pyr': 11. / 15, 'L23_PV': 11./14, 'L23_SOM': 1./12, 'L23_VIP': 0/17,
        'L4_Pyr': 18. / 18, 'L4_PV': 13./13, 'L4_SOM': 5./14, 'L4_VIP': 0/10,
        'L5_Pyr': 12. / 12, 'L5_PV': 15./15, 'L5_SOM': 0/10, 'L5_VIP': 0/6,
        'L6_Pyr': 8. / 10, 'L6_PV': 8./8., 'L6_SOM': 0/6, 'L6_VIP': 0/6
        }

    data['TC_Ji2016']['amplitude'] = {
        'L1': 425,
        'L23_Pyr': 129, 'L23_PV': 269, 'L23_SOM': 0, 'L23_VIP': 0,
        'L4_Pyr': 418, 'L4_PV': 962, 'L4_SOM': 20, 'L4_VIP': 24,
        'L5_Pyr': 195, 'L5_PV': 426, 'L5_SOM': 0, 'L5_VIP': 0,
        'L6_Pyr': 132, 'L6_PV': 426, 'L6_SOM': 0, 'L6_VIP': 0
        }


    # ------------------------------------------------------------------------------------------------------------------
    # load Thal -> A1 (Constantinople & Bruno, 2006,2015)
    # rat S1 (TC->L5/6 pmat+wmat)
    # prob of conn between thalamus and cortical layers / cell type
    # weight distribution: 0.2 - 1.2 mV
    
    data['TC_Cons2006_2015'] = {'prob': {}}
    data['TC_Cons2006_2015']['prob']['L23'] = 0.0
    data['TC_Cons2006_2015']['prob']['L4'] = 0.43
    data['TC_Cons2006_2015']['prob']['L5_IT'] = 0.17
    data['TC_Cons2006_2015']['prob']['L5_PT'] = 0.44
    data['TC_Cons2006_2015']['prob']['L6'] = 0.09

   # ------------------------------------------------------------------------------------------------------------------
    # Cruishcanck 2010 (mouse VB/TRN <-> Barrel cortex L4/5 RS/FS/LTS)
    # strong T -> RS and PV but not SOM
    # C -> VB and ITRN; ITRN -> VB (but not -> TRN) i.e. "CT activation produced strong inhibition in VB but not in TRN cells"
    # "CT projections outnumber TC projections and that CT synapses provide major input to thalamic neurons"

    # Data from fig 6B (avg EPSP amplitudes to optogen stim): 
    data['TC_Crui2010'] = {'epsp': {}, 'prob': {}}
    data['TC_Crui2010']['epsp']['L4_RS'] = 11
    data['TC_Crui2010']['epsp']['L4_FS'] = 20
    data['TC_Crui2010']['epsp']['L4_LTS'] = 2
    data['TC_Crui2010']['epsp']['L5_RS'] = 15
    data['TC_Crui2010']['epsp']['L5_FS'] = 25 
    data['TC_Crui2010']['epsp']['L5_LTS'] = 1

    data['TC_Crui2010']['prob']['TRN_TRN'] = 0.028
    data['TC_Crui2010']['prob']['TRN_VB'] = 0.08

    return data

# ----------------------------------------------------------------------------------------------------------------
# Params
# ----------------------------------------------------------------------------------------------------------------
Etypes = ['IT', 'ITS4', 'PT', 'CT']
Epops = ['IT2', 'IT3', 'ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']  # all layers

Itypes = ['PV', 'SOM', 'VIP', 'NGF']
Ipops = ['NGF1',                            # L1
        'PV2', 'SOM2', 'VIP2', 'NGF2',      # L2
        'PV3', 'SOM3', 'VIP3', 'NGF3',      # L3
        'PV4', 'SOM4', 'VIP4', 'NGF4',      # L4
        'PV5A', 'SOM5A', 'VIP5A', 'NGF5A',  # L5A  
        'PV5B', 'SOM5B', 'VIP5B', 'NGF5B',  # L5B
        'PV6', 'SOM6', 'VIP6', 'NGF6']  # L6 
        
Tpops = ['TC', 'TCM', 'HTC', 'IRE', 'IREM', 'TI', 'TIM']

layer = {'1': [0.00, 0.05], '2': [0.05, 0.08], '3': [0.08, 0.475], '4': [0.475,0.625], '5A': [0.625,0.667], '5B': [0.667,0.775], '6': [0.775,1], 'thal': [1.2, 1.4]} # 


# initialize prob and weight matrices
# format: pmat[presynaptic_pop][postsynaptic_pop] 
pmat = {}  # probability of connection matrix
lmat = {}  # length constant (lambda) for exp decaying prob conn (um) matrix
wmat = {}  # connection weight matrix = unitary conn somatic PSP (mV)
secmat = {}  # target section matrix

for p in Epops + Ipops + Tpops:
    pmat[p] = {}
    lmat[p] = {}
    wmat[p] = {}
    
# Load exp data
data = loadData()

# Set source of conn data
connDataSource = {}
connDataSource['E->E/I'] = 'Allen_BBP' #'Allen_V1' #'BBP_S1'  # 'Allen_V1' 
connDataSource['I->E/I'] = 'Allen_custom' #'custom_A1' #'BBP_S1'  # 'Allen_V1' 


# --------------------------------------------------
## E -> E/I 
# --------------------------------------------------

# --------------------------------------------------
## Probabilities, length constants (lambda), weights (=unitary conn somatic PSP amplitude), and target sections

'''
Probabilities are distance dependent based on:
A0 * np.exp(- (intersomatic_distance / lambda) ** 2)

where A0 is the probability of connection and distance 0um
and lambda is the length constant
'''

# start with base data from Allen V1
if connDataSource['E->E/I'] ==  'Allen_V1': 
    for pre in Epops:
        for post in Epops+Ipops:
            proj = '%s-%s' % (data['Allen_V1']['pops'][pre], data['Allen_V1']['pops'][post])
            pmat[pre][post] = data['Allen_V1']['connProb'][proj]['A0']
            lmat[pre][post] = data['Allen_V1']['connProb'][proj]['sigma']
            wmat[pre][post] = data['Allen_V1']['connWeight'][proj]

# use BBP S1 instead (has more cell-type specificity)
elif connDataSource['E->E/I'] == 'BBP_S1': 
    for pre in Epops:
        for post in Epops+Ipops:
            proj = '%s:%s' % (data['BBP_S1']['pops'][pre], data['BBP_S1']['pops'][post])
            if proj in data['BBP_S1']['connProb']:
                pmat[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability']/100.
                wmat[pre][post] = data['BBP_S1']['connWeight'][proj]['epsp_mean']
            else:
                pmat[pre][post] = 0.
                wmat[pre][post] = 0.

# use Allen but update with BBP cell-type specificity
elif connDataSource['E->E/I'] ==  'Allen_BBP': 
    for pre in Epops:
        for post in Epops+Ipops:
            proj = '%s-%s' % (data['Allen_V1']['pops'][pre], data['Allen_V1']['pops'][post])
            pmat[pre][post] = data['Allen_V1']['connProb'][proj]['A0']
            lmat[pre][post] = data['Allen_V1']['connProb'][proj]['sigma']
            wmat[pre][post] = data['Allen_V1']['connWeight'][proj]

    # List of pops to update based on BBP ratios of cell subtypes (key = pop to update; value = pop to use as reference)
    updatePopConnUsingBBP = {'VIP2': 'PV2', 'VIP3': 'PV3', 'VIP4': 'PV4', 'VIP5A': 'PV5A', 'VIP5B': 'PV5A', 'VIP6': 'PV6',
                            'ITS4': 'ITP4',
                            'CT5A': 'IT5A',
                            'CT5B': 'IT5B',
                            'PT5B': 'IT5B',
                            'CT6': 'IT6'}
    
    basedOnBBPCT6 = [] #'CT5A', 'CT5B']  # pops based on BBP CT6 (since BBP doesn't have CT5A and CT5B) so treated differently (NOT USED)

    fixVerbose = 1 #False  # print info messages
    
    ## update all fix pop (e.g. VIP) by making proportional to ref pop (e.g. PV): e.g. VIP_Allen = (VIP_BBP/PV_BBP) * PV_Allen
    for fixpop, refpop in updatePopConnUsingBBP.items():
        if fixVerbose:
            print('\nUpdating conn probability of pop %s using as reference BBP conn probability ratio of %s:%s ...' % (fixpop, fixpop, refpop))
        
        # E -> fixpop
        for pre in Epops:
            if pre not in updatePopConnUsingBBP:
                projAllen_ref = '%s-%s' % (data['Allen_V1']['pops'][pre], data['Allen_V1']['pops'][refpop])
                if fixpop in basedOnBBPCT6:
                    ## Make CT5A <-> L5A E/I and CT5B <-> L5B E/I == CT6 <-> L6 E/I (so based on local conn) (NOT USED)
                    projPre = pre.replace('5A', '6').replace('5B', '6').replace('PT6', 'PT5B')
                else:
                    projPre = pre
                projBBP_ref = '%s:%s' % (data['BBP_S1']['pops'][projPre], data['BBP_S1']['pops'][refpop])
                projBBP_fix = '%s:%s' % (data['BBP_S1']['pops'][projPre], data['BBP_S1']['pops'][fixpop])

            # conn probs 
            ref_Allen = data['Allen_V1']['connProb'][projAllen_ref]['A0'] if projAllen_ref in data['Allen_V1']['connProb'] else 0.
            ref_BBP = data['BBP_S1']['connProb'][projBBP_ref]['A0'] if projBBP_ref in data['BBP_S1']['connProb'] else 0.
            fix_BBP = data['BBP_S1']['connProb'][projBBP_fix]['A0'] if projBBP_fix in data['BBP_S1']['connProb'] else 0.
            if ref_BBP > 0. and fix_BBP > 0.:
                if fixVerbose:
                    print(' Prob %s->%s:'%(pre, fixpop), 'ref_BBP: %.2f'%(ref_BBP), 'fix_BBP: %.2f'%(fix_BBP), 'ref_Allen: %.2f'%(ref_Allen), 'fix_Allen: %.2f'%((fix_BBP/ref_BBP) * ref_Allen))
                pmat[pre][fixpop] = (fix_BBP / ref_BBP) * ref_Allen
 
        # fixpop -> E/I
        for post in Epops+Ipops:
            projAllen_ref = '%s-%s' % (data['Allen_V1']['pops'][refpop], data['Allen_V1']['pops'][post])
            if fixpop in basedOnBBPCT6:
                ## Make CT5A <-> L5A E/I and CT5B <-> L5B E/I == CT6 <-> L6 E/I (so based on local conn)  (NOT USED)
                projPost = post.replace('5A', '6').replace('5B', '6').replace('PT6', 'PT5B')
            else:
                projPost = post
            projBBP_ref = '%s:%s' % (data['BBP_S1']['pops'][refpop], data['BBP_S1']['pops'][projPost])
            projBBP_fix = '%s:%s' % (data['BBP_S1']['pops'][fixpop], data['BBP_S1']['pops'][projPost])

            # conn probs 
            ref_Allen = data['Allen_V1']['connProb'][projAllen_ref]['A0'] if projAllen_ref in data['Allen_V1']['connProb'] else 0.
            ref_BBP = data['BBP_S1']['connProb'][projBBP_ref]['A0'] if projBBP_ref in data['BBP_S1']['connProb'] else 0.
            fix_BBP = data['BBP_S1']['connProb'][projBBP_fix]['A0'] if projBBP_fix in data['BBP_S1']['connProb'] else 0.
            if ref_BBP > 0. and fix_BBP > 0.:
                if fixVerbose:
                    print(' Prob %s->%s:'%(fixpop,post), 'ref_BBP: %.2f'%(ref_BBP), 'fix_BBP: %.2f'%(fix_BBP), 'ref_Allen: %.2f'%(ref_Allen), 'fix_Allen: %.2f'%((fix_BBP/ref_BBP) * ref_Allen))
                pmat[fixpop][post] = (fix_BBP / ref_BBP) * ref_Allen
        



# --------------------------------------------------
## I -> E
# --------------------------------------------------

bins = {}
bins['inh'] = [[0.0, 0.37], [0.37, 0.8], [0.8,1.0]]

# --------------------------------------------------
## Probabilities 

'''
- I->E/I connections data from:
-- Tremblay, 2016, mouse S1;
-- Budinger, 2018, multispecies A1; 
-- Naka et al 2016, mouse multiple regions; 
-- Kato et al 2017, mouse A1;

- Local, intralaminar only; all-to-all but distance-based; high weights
- Although evidence for L2/3,4,6 -> L5A/B, strongest is intralaminar (Naka16)
- Consistent with Allen V1, except Allen doesn't show strong L2/3 I -> L5 E
'''

if connDataSource['I->E/I'] == 'custom_A1': 
    # I->E particularities:
    ## inh cells target apical dends (not strictly local to layer) 
    ## L1 + L2/3 NGF -> L2/3 E (prox apic) + L5 E (tuft) -- specify target dend in subConnParams
    ## upper layer SOM, VIP, NGF project strongly to deeper E cells (Kato 2017) with exp-decay distance-dep

    for pre in ['SOM', 'VIP', 'NGF']:
        pmat[pre] = {}
        pmat[pre]['E'] = np.array([[1.0, 1.0, 1.0],  # from L1+L2/3+L4 to all layers 
                                    [0.25, 1.0, 1.0],  # from L5A+L5B to all layers
                                    [0.25, 0.25, 1.0]])  # from L6 to all layers

    ## upper layer PV project weakly to deeper E cells (mostly intralaminar) (Kato 2017) with exp-decay distance-dep
    ## although Naka 2016 shows L2/3 PV -> L5 Pyr
    pmat['PV'] = {}
    pmat['PV']['E'] = np.array([[1.0, 0.5, 0.25],  # from L1+L2/3+L4 to all layers 
                                [0.25, 1.0, 0.5],  # from L5A+L5B to all layers
                                [0.1, 0.25, 1.0]])  # from L6 to all layers

    # VIP -> E (very low; 3/42; Pi et al 2013) 
    pmat['VIP']['E'] *= 0.1

    # --------------------------------------------------
    ## Weights  (=unitary conn somatic PSP amplitude)
    IEweight = 1.0
    for pre in Ipops:
        for post in Ipops:
            wmat[pre][post] = IEweight

# Allen V1
elif connDataSource['I->E/I'] == 'Allen_V1': 
    for pre in Ipops:
        for post in Epops:
            proj = '%s-%s' % (data['Allen_V1']['pops'][pre], data['Allen_V1']['pops'][post])
            pmat[pre][post] = data['Allen_V1']['connProb'][proj]['A0']
            lmat[pre][post] = data['Allen_V1']['connProb'][proj]['sigma']
            wmat[pre][post] = data['Allen_V1']['connWeight'][proj]

# use BBP S1 instead? (has more cell-type specificity)
elif connDataSource['I->E/I'] ==  'BBP_S1': 
    for pre in Ipops:
        for post in Epops:
            proj = '%s:%s' % (data['BBP_S1']['pops'][pre], data['BBP_S1']['pops'][post])
            if proj in data['BBP_S1']['connProb']:
                pmat[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability']/100.0
                wmat[pre][post] = data['BBP_S1']['connWeight'][proj]['epsp_mean']
            else:
                pmat[pre][post] = 0.
                wmat[pre][post] = 0.

# Allen V1
elif connDataSource['I->E/I'] == 'Allen_custom': 
    for pre in Ipops:
        for post in Epops:
            proj = '%s-%s' % (data['Allen_V1']['pops'][pre], data['Allen_V1']['pops'][post])
            pmat[pre][post] = data['Allen_V1']['connProb'][proj]['A0']
            lmat[pre][post] = data['Allen_V1']['connProb'][proj]['sigma']
            wmat[pre][post] = data['Allen_V1']['connWeight'][proj]


## Update L2/3 SOM,PV -> deeper E depth-dependent (Kato 2017; Fig 8B); supported by L2/3 PV -> L5 E (Naka 2016, Fig 2)
## L2 and L3 combined in Kato 2017 fig; also L2 very thin in A1
SOM23_E = {'L1': 0.92, 'L23u': 0.91, 'L23l': 0.84, 'L4': 0.67, 'L5A': 0.48, 'L5B': 0.28, 'L6': 0.08}
PV23_E = {'L1': 0.69, 'L23u': 0.89, 'L23l': 0.72, 'L4': 0.26, 'L5A': 0.05, 'L5B': 0.03, 'L6': 0.01}

SOM23_E23_modelToKato2017Ratio = pmat['SOM3']['IT3'] / SOM23_E['L23u']
PV23_E23_modelToKato2017Ratio = pmat['PV3']['IT3'] / PV23_E['L23u']

verbose = False
for post in ['ITP4', 'ITS4', 'IT5A', 'CT5A', 'IT5B', 'PT5B', 'CT5B', 'IT6', 'CT6']:
        layer = post[-1] if post[-1] in ['4', '6'] else post[-2:]  # '4', '5A', '5B', or '6'
        if verbose:
            print('Before %s -> %s = %.2f' % ('SOM3', post, pmat['SOM3'][post]))
            print('Before %s -> %s = %.2f' % ('PV3', post, pmat['PV3'][post]))
        pmat['SOM2'][post] = pmat['SOM3'][post] = SOM23_E23_modelToKato2017Ratio * SOM23_E['L'+layer]
        pmat['PV2'][post] = pmat['PV3'][post] = PV23_E23_modelToKato2017Ratio * PV23_E['L' + layer]
        if verbose:
            print('After %s -> %s = %.2f' % ('SOM3', post, pmat['SOM3'][post]))
            print('After %s -> %s = %.2f' % ('PV3', post, pmat['PV3'][post]))


## Update strong L1 NGF -> L5; currently at 0.148, so moderately strong already

## Update VIP -> E very low (3/42; Pi et al 2013); but L2/3 VIP -> L5 (Naka 2016, Fig 2)
## make same as SOM but multiply by ration of VIP3->E5 (Pi 2013) / SOM3->E5 (Kato 2017)
Pi2013_VIP_E5 = 3. / 42
ratio_VIP_SOM = Pi2013_VIP_E5 / SOM23_E['L5A']
verbose = 0
for prelayer in ['2', '3', '4', '5A', '5B', '6']:
    for post in Epops:
        pmat['VIP' + prelayer][post] = pmat['SOM' + prelayer][post] * ratio_VIP_SOM
        if verbose: print('%s -> %s = %.2f' % ('VIP'+prelayer, post, pmat['VIP'+prelayer][post]))

# --------------------------------------------------
## I -> I
# --------------------------------------------------

# --------------------------------------------------
## Probabilities 

'''
- NGF -> all I, local/intralaminar (+ distance-dep)
- VIP -> SOM (strong; 14/18), PV (weak; 4/15), VIP (very weak -- remove?) (Pi 2013)
- SOM -> FS+VIP (strong); SOM (weak -- remove?)  (Naka et al 2016;Tremblay, 2016; Sohn, 2016)
- PV  -> PV (strong); SOM+VIP (weak -- remove?) (Naka et al 2016;Tremblay, 2016; Sohn, 2016)
- Generally consistent with the more detailed Allen V1 I->I
'''

if connDataSource['I->E/I'] == 'custom_A1':
    
### I->I particularities
    for pre in Itypes:
        pmat[pre]
        for post in Itypes:
            pmat[pre][post] = 1.0

    # NGF -> all I, local/intralaminar (+ distance-dep)
    # no change required

    strong = 1.0
    weak = 0.35  # = ~ normalized strong/weak = (4/15) / (14/18)
    veryweak = 0.1

    # VIP -> SOM (strong; 14/18), PV (weak; 4/15), VIP (very weak -- remove?) (Pi 2013)
    pmat['VIP']['SOM'] = strong
    pmat['VIP']['PV'] = weak
    pmat['VIP']['VIP'] = veryweak
    pmat['VIP']['NGF'] = weak  # unknown; assume weak

    # SOM -> FS+VIP (strong); SOM (weak -- remove?)  (Naka et al 2016;Tremblay, 2016; Sohn, 2016)
    pmat['SOM']['PV'] = strong
    pmat['SOM']['VIP'] = strong
    pmat['SOM']['SOM'] = weak
    pmat['SOM']['NGF'] = weak  # unknown; assume weak
    
    # PV  -> PV (strong); SOM+VIP (weak -- remove?) (Naka et al 2016;Tremblay, 2016; Sohn, 2016)
    pmat['PV']['PV'] = strong
    pmat['PV']['SOM'] = weak
    pmat['PV']['VIP'] = weak
    pmat['PV']['NGF'] = weak  # unknown; assume weak

    # NGF -> I; unknown, assume weak since most data focuses on NGF -> E
    pmat['NGF']['PV'] = weak
    pmat['NGF']['SOM'] = weak
    pmat['NGF']['VIP'] = weak
    pmat['NGF']['NGF'] = weak  

    # --------------------------------------------------
    ## Weights  (=unitary conn somatic PSP amplitude)
    IIweight = 1.0
    for pre in Ipops:
        for post in Ipops:
            wmat[pre][post] = IIweight


# Allen V1
elif connDataSource['I->E/I'] ==  'Allen_V1': 
    for pre in Ipops:
        for post in Ipops:
            proj = '%s-%s' % (data['Allen_V1']['pops'][pre], data['Allen_V1']['pops'][post])
            pmat[pre][post] = data['Allen_V1']['connProb'][proj]['A0']
            lmat[pre][post] = data['Allen_V1']['connProb'][proj]['sigma']
            wmat[pre][post] = data['Allen_V1']['connWeight'][proj]

# use BBP S1 instead? (has more cell-type specificity)
elif connDataSource['I->E/I'] ==  'BBP_S1': 
    for pre in Ipops:
        for post in Ipops:
            proj = '%s:%s' % (data['BBP_S1']['pops'][pre], data['BBP_S1']['pops'][post])
            if proj in data['BBP_S1']['connProb']:
                pmat[pre][post] = data['BBP_S1']['connProb'][proj]['connection_probability']/100.0
                wmat[pre][post] = data['BBP_S1']['connWeight'][proj]['epsp_mean']
            else:
                pmat[pre][post] = 0.
                wmat[pre][post] = 0.

# Allen V1 customized
elif connDataSource['I->E/I'] ==  'Allen_custom': 
    for pre in Ipops:
        for post in Ipops:
            proj = '%s-%s' % (data['Allen_V1']['pops'][pre], data['Allen_V1']['pops'][post])
            pmat[pre][post] = data['Allen_V1']['connProb'][proj]['A0']
            lmat[pre][post] = data['Allen_V1']['connProb'][proj]['sigma']
            wmat[pre][post] = data['Allen_V1']['connWeight'][proj]

    # VIP uses by default PV weights; made the following changes:
    # VIP -> SOM = strong (but PV -> SOM = weak)    
    # VIP -> PV = weak (but PV -> PV = strong)   
    # VIP -> VIP = weak/veryweak (but PV -> PV = strong) 
    verbose = 0
    for prelayer in ['2', '3', '4', '5A', '5B', '6']:
        for post in Ipops:
            if 'SOM' in post:
                pmat['VIP' + prelayer][post] = pmat['PV' + prelayer][post.replace('SOM', 'PV')] # replace VIP->SOM with VIP->PV
            elif 'PV' in post:
                pmat['VIP' + prelayer][post] = pmat['PV' + prelayer][post.replace('PV', 'SOM')] # replace VIP->PV with VIP->SOM
            elif 'VIP' in post:
                pmat['VIP' + prelayer][post] = pmat['PV' + prelayer][post.replace('VIP', 'SOM')] # replace VIP->VIP with VIP->SOM
                
            if verbose: print('%s -> %s = %.2f' % ('VIP'+prelayer, post, pmat['VIP' + prelayer][post]))

# --------------------------------------------------
# Delays
## Make distance-dependent for now


# --------------------------------------------------
## INTRATHALAMIC (from old model; based on Bonj12, Bazhenov https://www.jneurosci.org/content/32/15/5250.full Lakatos observations; not checked)
# --------------------------------------------------

# --------------------------------------------------
## Probabilities 
## note: missing values mean 0 probability
pmat['TC']['TC'] =	    0.0 # 0.1
pmat['HTC']['HTC'] =	0.0 #0.1
pmat['TC']['HTC'] =	    0.0 # 0.1
pmat['HTC']['TC'] =	    0.0 # 0.1
pmat['TCM']['TCM'] =	0.0 #0.1
pmat['IRE']['IRE'] =	0.1  # > data['TC_Crui2010']['prob']['TRN_TRN'] = 0.028
pmat['IREM']['IREM'] =	0.1  # > data['TC_Crui2010']['prob']['TRN_TRN'] = 0.028
pmat['IRE']['IREM'] =	0.1  # > data['TC_Crui2010']['prob']['TRN_TRN'] = 0.028
pmat['IREM']['IRE'] =	0.1  # > data['TC_Crui2010']['prob']['TRN_TRN'] = 0.028
pmat['TC']['IREM'] =	0.2
pmat['HTC']['IREM'] =	0.2
pmat['IREM']['TC'] =	0.1  # ~= data['TC_Crui2010']['prob']['TRN_VB'] = 0.08
pmat['IREM']['HTC'] =	0.1  # ~= data['TC_Crui2010']['prob']['TRN_VB'] = 0.08
pmat['TCM']['IRE'] =	0.2
pmat['IRE']['TCM'] =	0.1
pmat['TC']['IRE'] =	    0.4  
pmat['HTC']['IRE'] =	0.4  
pmat['IRE']['TC'] =	    0.3  # > data['TC_Crui2010']['prob']['TRN_VB'] = 0.08
pmat['IRE']['HTC'] =	0.3  # > data['TC_Crui2010']['prob']['TRN_VB'] = 0.08
pmat['TCM']['IREM'] =	0.4
pmat['IREM']['TCM'] =	0.3  # > data['TC_Crui2010']['prob']['TRN_VB'] = 0.08
pmat['TI']['TC']    =   0.21  # TI values from Serkov 1996
pmat['TI']['HTC']   =   0.21
pmat['TI']['TCM']   =   0.21
pmat['TC']['TI']    =   0.01
pmat['HTC']['TI']   =   0.01
pmat['TCM']['TI']   =   0.09
pmat['IRE']['TI']   =   0.09
pmat['IREM']['TI']  =	0.3 
pmat['TI']['TI']    =  0.53
pmat['TIM']['TC']   =   0.21  # TIM values from Serkov 1996
pmat['TIM']['HTC']   =   0.21
pmat['TIM']['TCM']   =   0.21
pmat['TC']['TIM']    =   0.01
pmat['HTC']['TIM']   =   0.01
pmat['TCM']['TIM']   =   0.09
pmat['IRE']['TIM']   =   0.09
pmat['IREM']['TIM']  =	0.3 
pmat['TIM']['TIM']    =   0.53

# --------------------------------------------------
## Weights  (=unitary conn somatic PSP amplitude)
wmat['HTC']['HTC'] =    0.1
wmat['HTC']['TC'] =     0.1
wmat['TC']['HTC'] =     0.1
wmat['TC']['TC'] =  	0.1
wmat['TCM']['TCM'] =    0.1
wmat['IRE']['IRE'] =    1.5
wmat['IREM']['IREM'] =  1.5
wmat['IRE']['IREM'] =   1.5
wmat['IREM']['IRE'] =   1.5
wmat['TC']['IREM'] =    0.23
wmat['HTC']['IREM'] =   0.123
wmat['IREM']['TC'] =    0.83
wmat['IREM']['HTC'] =   0.83
wmat['TCM']['IRE'] =    0.2
wmat['IRE']['TCM'] =    0.83
wmat['TC']['IRE'] = 	0.2
wmat['HTC']['IRE'] =    0.2
wmat['IRE']['TC'] =     0.83
wmat['IRE']['HTC'] =    0.83
wmat['TCM']['IREM'] =   0.2
wmat['IREM']['TCM'] =   0.83
wmat['TI']['TC']    =   0.83  # for TI using same values as RE
wmat['TI']['HTC']   =   0.83
wmat['TI']['TCM']   =   0.83
wmat['TC']['TI']    =   0.2
wmat['HTC']['TI']   =   0.2
wmat['TCM']['TI']   =   0.2
wmat['TCM']['TI']   =   0.2
wmat['IRE']['TI']   =   0.5
wmat['IREM']['TI']  =	0.5 
wmat['TI']['TI']    =   0.5
wmat['TIM']['TC']    =   0.83  # for TIM using same values as RE
wmat['TIM']['HTC']   =   0.83
wmat['TIM']['TCM']   =   0.83
wmat['TC']['TIM']    =   0.2
wmat['HTC']['TIM']   =   0.2
wmat['TCM']['TIM']   =   0.2
wmat['TCM']['TIM']   =   0.2
wmat['IRE']['TIM']   =   0.5
wmat['IREM']['TIM']  =   0.5 
wmat['TIM']['TIM']    =   0.5


# --------------------------------------------------
## CORTICOTHALAMIC (from old model; partly from Bonj12, Bazhenov https://www.jneurosci.org/content/32/15/5250.full and discuss with Lakatos)
# --------------------------------------------------

# --------------------------------------------------
## Probabilities 
pmat['CT5A']['TC']	= 0.1
pmat['CT5A']['HTC']	= 0.1
pmat['CT5A']['IRE']	= 0.1
pmat['CT5A']['TI']	= 0.05
pmat['CT5B']['TC']	= 0.1
pmat['CT5B']['HTC']	= 0.1
pmat['CT5B']['IRE'] = 0.1
pmat['CT5B']['TI']  = 0.05
pmat['CT6']['TC']	= 0.1
pmat['CT6']['HTC']	= 0.1
pmat['CT6']['IRE']	= 0.1
pmat['CT6']['TI']	= 0.05

pmat['IT5B']['TCM']	= 0.1
pmat['PT5B']['TCM']	= 0.1
pmat['IT5B']['IREM'] = 0.1
pmat['PT5B']['IREM'] = 0.1
pmat['IT5B']['TIM']	= 0.05
pmat['PT5B']['TIM']	= 0.05

# --------------------------------------------------
## Weights  (=unitary conn somatic PSP amplitude)
wmat['CT5A']['TC']	= 0.7
wmat['CT5A']['HTC']	= 0.7
wmat['CT5A']['IRE']	= 0.23
wmat['CT5A']['TI']	= 0.23
wmat['CT5B']['TC']	= 0.7
wmat['CT5B']['HTC']	= 0.7
wmat['CT5B']['IRE']	= 0.23
wmat['CT5B']['IRE']	= 0.23
wmat['CT5B']['TI']	= 0.23
wmat['CT6']['TC']	= 0.7
wmat['CT6']['HTC']	= 0.7
wmat['CT6']['IRE']	= 0.23
wmat['CT6']['TI']	= 0.23

wmat['IT5B']['TCM']	= 0.7
wmat['PT5B']['TCM']	= 0.7
wmat['IT5B']['IREM'] = 0.23
wmat['PT5B']['IREM'] = 0.23
wmat['IT5B']['TIM']	= 0.23
wmat['PT5B']['TIM']	= 0.23

# --------------------------------------------------
## CORE THALAMOCORTICAL (from old model; partly from Bonj12, Bazhenov https://www.jneurosci.org/content/32/15/5250.full and discuss with Lakatos)
# --------------------------------------------------

# factor to convert amplitudes to probability; normalized max value in  Ji 2016 (amplitude) based on max value in Cons2006 (probability)
normProb =  data['TC_Cons2006_2015']['prob']['L4'] / data['TC_Ji2016']['amplitude']['L4_PV']

# Use conn data from Ji et al 2016; layer and cell-type specific (consistent with Constantinople)
# - use 'amplitude' data which represents conn strength (prob * unitary conn) and map to model prob conn, keeping weight fixed
# - use SOM values for NGF 
# - added inputs to IT3 and PV3 following Ji2016 data and Lakatos suggestion (didn't include SOM, VIP, NGF since very low)

pmat['TC']['IT3']       = data['TC_Ji2016']['amplitude']['L23_Pyr'] * normProb  # orig value: 0.1   
pmat['HTC']['IT3']      = data['TC_Ji2016']['amplitude']['L23_Pyr'] * normProb  # orig value: 0.1   
pmat['TC']['ITP4']      = data['TC_Ji2016']['amplitude']['L4_Pyr'] * normProb  # orig value: 0.25
pmat['HTC']['ITP4']     = data['TC_Ji2016']['amplitude']['L4_Pyr'] * normProb  # orig value: 0.25
pmat['TC']['ITS4']      = data['TC_Ji2016']['amplitude']['L4_Pyr'] * normProb  # orig value: 0.25
pmat['HTC']['ITS4']     = data['TC_Ji2016']['amplitude']['L4_Pyr'] * normProb  # orig value: 0.25
pmat['TC']['PT5B']      = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.1   #thalfctr
pmat['HTC']['PT5B']     = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.1   
pmat['TC']['IT5A']      = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.1   
pmat['HTC']['IT5A']     = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.1   
pmat['TC']['IT5B']      = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.1   
pmat['HTC']['IT5B']     = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.1   
pmat['TC']['IT6']       = data['TC_Ji2016']['amplitude']['L6_Pyr'] * normProb  # orig value: 0.15  
pmat['HTC']['IT6']      = data['TC_Ji2016']['amplitude']['L6_Pyr'] * normProb  # orig value: 0.15  
pmat['TC']['CT5A']      = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.15  
pmat['HTC']['CT5A']     = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.15  
pmat['TC']['CT5B']      = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.15  
pmat['HTC']['CT5B']     = data['TC_Ji2016']['amplitude']['L4_Pyr'] * normProb  # orig value: 0.15  
pmat['TC']['CT6']       = data['TC_Ji2016']['amplitude']['L6_Pyr'] * normProb  # orig value: 0.15  
pmat['HTC']['CT6']      = data['TC_Ji2016']['amplitude']['L6_Pyr'] * normProb  # orig value: 0.15  
pmat['TC']['PV3']       = data['TC_Ji2016']['amplitude']['L23_PV'] * normProb  # orig value: 0.25
pmat['HTC']['PV3']      = data['TC_Ji2016']['amplitude']['L23_PV'] * normProb  # orig value: 0.25
pmat['TC']['PV4']       = data['TC_Ji2016']['amplitude']['L4_PV'] * normProb  # orig value: 0.25
pmat['HTC']['PV4']      = data['TC_Ji2016']['amplitude']['L4_PV'] * normProb  # orig value: 0.25
pmat['TC']['SOM4']      = data['TC_Ji2016']['amplitude']['L4_SOM'] * normProb  # orig value: 0.25
pmat['HTC']['SOM4']     = data['TC_Ji2016']['amplitude']['L4_SOM'] * normProb  # orig value: 0.25
pmat['TC']['NGF4']      = data['TC_Ji2016']['amplitude']['L4_SOM'] * normProb  # orig value: 0.25
pmat['HTC']['NGF4']     = data['TC_Ji2016']['amplitude']['L4_SOM'] * normProb  # orig value: 0.25	
pmat['TC']['PV5A']      = data['TC_Ji2016']['amplitude']['L5_PV'] * normProb  # orig value: 0.1   
pmat['HTC']['PV5A']     = data['TC_Ji2016']['amplitude']['L5_PV'] * normProb  # orig value: 0.1   
pmat['TC']['SOM5A']     = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb  # orig value: 0.1   
pmat['HTC']['SOM5A']    = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb  # orig value: 0.1  
pmat['TC']['NGF5A']     = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb	# orig value: 0.1 #	
pmat['HTC']['NGF5A']    = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb	# orig value: 0.1 #	
pmat['TC']['PV5B']      = data['TC_Ji2016']['amplitude']['L5_PV'] * normProb  # orig value: 0.1   
pmat['HTC']['PV5B']     = data['TC_Ji2016']['amplitude']['L5_PV'] * normProb  # orig value: 0.1   
pmat['TC']['SOM5B']     = data['TC_Ji2016']['amplitude']['L5_PV'] * normProb  # orig value: 0.1   
pmat['HTC']['SOM5B']    = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb  # orig value: 0.1  
pmat['TC']['NGF5B']     = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb  # orig value: 0.1 #	
pmat['HTC']['NGF5B']    = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb  # orig value: 0.1 #	
pmat['TC']['PV6']       = data['TC_Ji2016']['amplitude']['L6_PV'] * normProb  # orig value: 0.15  
pmat['HTC']['PV6']      = data['TC_Ji2016']['amplitude']['L6_PV'] * normProb  # orig value: 0.15  
pmat['TC']['SOM6']      = data['TC_Ji2016']['amplitude']['L6_SOM'] * normProb  # orig value: 0.15  
pmat['HTC']['SOM6']     = data['TC_Ji2016']['amplitude']['L6_SOM'] * normProb  # orig value: 0.15  
pmat['TC']['NGF6']      = data['TC_Ji2016']['amplitude']['L6_SOM'] * normProb  # orig value: 0.15 #	
pmat['HTC']['NGF6']     = data['TC_Ji2016']['amplitude']['L6_SOM'] * normProb  # orig value: 0.15 #	


# --------------------------------------------------
# Weights  (=unitary conn somatic PSP amplitude)

# don't have specific exp data so use 0.5mV, which is approx avg in Constantinople, 2016  
TCweight = 0.5

wmat['TC']['IT3']      = TCweight  #	
wmat['HTC']['IT3']     = TCweight  #	
wmat['TC']['ITP4']      = TCweight  # 0.6
wmat['HTC']['ITP4']     = TCweight  # 0.6
wmat['TC']['ITS4']      = TCweight  # 0.6
wmat['HTC']['ITS4']     = TCweight  # 0.6
wmat['TC']['PT5B']      = TCweight  #	0.6  # [TC][E5B] / pmat[TC][E4]	
wmat['HTC']['PT5B']     = TCweight  # 0.6  #* pmat[TC][E5B] / pmat[TC][E4]
wmat['TC']['IT5A']      = TCweight  #	0.6  #* pmat[TC][E5R] / pmat[TC][E4]	
wmat['HTC']['IT5A']     = TCweight  #	0.6  #* pmat[TC][E5R] / pmat[TC][E4]	
wmat['TC']['IT5B']      = TCweight  #	0.6  #* pmat[TC][E5R] / pmat[TC][E4]	
wmat['HTC']['IT5B']     = TCweight  #	0.6  #* pmat[TC][E5R] / pmat[TC][E4]	
wmat['TC']['IT6']       = TCweight  #	0.6  #* pmat[TC][E6] / pmat[TC][E4]	
wmat['HTC']['IT6']      = TCweight  # 0.6  #* pmat[TC][E6] / pmat[TC][E4]	
wmat['TC']['CT5A']      = TCweight  #	0.6  #* pmat[TC][E6] / pmat[TC][E4]
wmat['HTC']['CT5A']     = TCweight  #	0.6  #* pmat[TC][E6] / pmat[TC][E4]	
wmat['TC']['CT5B']      = TCweight  #	0.6  #* pmat[TC][E6] / pmat[TC][E4]
wmat['HTC']['CT5B']     = TCweight  #	0.6  #* pmat[TC][E6] / pmat[TC][E4]	
wmat['TC']['CT6']       = TCweight  #	0.6  #* pmat[TC][E6] / pmat[TC][E4]
wmat['HTC']['CT6']      = TCweight  #	0.6  #* pmat[TC][E6] / pmat[TC][E4]	
wmat['TC']['PV3']       = TCweight  #	0.23 #	
wmat['HTC']['PV3']      = TCweight  #	0.23 #	
wmat['TC']['PV4']       = TCweight  #	0.23 #	
wmat['HTC']['PV4']      = TCweight  #	0.23 #	
wmat['TC']['SOM4']      = TCweight  #	0.23 #	
wmat['HTC']['SOM4']     = TCweight  #	0.23 #	
wmat['TC']['NGF4']      = TCweight  #	0.23 #	
wmat['HTC']['NGF4']     = TCweight  #	0.23 #	
wmat['TC']['PV5A']      = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['HTC']['PV5A']     = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['TC']['SOM5A']     = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['HTC']['SOM5A']    = TCweight  # 0.23  # * pmat[TC][I5] / pmat[TC][I4]	
wmat['TC']['NGF5A']     = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['HTC']['NGF5A']    = TCweight  # 0.23  # * pmat[TC][I5] / pmat[TC][I4]	
wmat['TC']['PV5B']      = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['HTC']['PV5B']     = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['TC']['SOM5B']     = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['HTC']['SOM5B']    = TCweight  # 0.23  # * pmat[TC][I5] / pmat[TC][I4]	
wmat['TC']['NGF5B']     = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['HTC']['NGF5B']    = TCweight  #	0.23 # * pmat[TC][I5] / pmat[TC][I4]	
wmat['TC']['PV6']       = TCweight  #	0.23 # * pmat[TC][I6] / pmat[TC][I4]	
wmat['HTC']['PV6']      = TCweight  #	0.23 # * pmat[TC][I6] / pmat[TC][I4]	
wmat['TC']['SOM6']      = TCweight  #	0.23 # * pmat[TC][I6] / pmat[TC][I4]	
wmat['HTC']['SOM6']     = TCweight  #	0.23 # * pmat[TC][I6] / pmat[TC][I4]
wmat['TC']['NGF6']      = TCweight  #	0.23 # * pmat[TC][I6] / pmat[TC][I4]	
wmat['HTC']['NGF6']     = TCweight  #	0.23 # * pmat[TC][I6] / pmat[TC][I4]

# --------------------------------------------------
# Length constant (matrix more broadly tunes spatially than core)
lmat['TC'] = {}
lmat['HTC'] = {}

for post in wmat['TC'].keys():
    lmat['TC'][post] = 500  # um - arbitrary value

for post in wmat['HTC'].keys():
    lmat['HTC'][post] = 1000  # um - arbitrary value

# --------------------------------------------------
## MATRIX THALAMOCORTICAL (from old model; partly from Bonj12, Bazhenov https://www.jneurosci.org/content/32/15/5250.full and discuss with Lakatos)
# --------------------------------------------------

# --------------------------------------------------
## Probabilities 

# Use conn data from Ji et al 2016; layer and cell-type specific
# use SOM values for NGF
pmat['TCM']['IT2']	    = data['TC_Ji2016']['amplitude']['L23_Pyr'] * normProb  # orig value: 0.25
pmat['TCM']['IT3']	    = data['TC_Ji2016']['amplitude']['L23_Pyr'] * normProb  # orig value: 0.25
pmat['TCM']['IT5A']	    = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.15  
pmat['TCM']['IT5B']	    = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.15  
pmat['TCM']['PT5B']	    = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.15  
pmat['TCM']['IT6']	    = data['TC_Ji2016']['amplitude']['L6_Pyr'] * normProb  # orig value: 0.05  
pmat['TCM']['CT5A']     = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.05  
pmat['TCM']['CT5B']     = data['TC_Ji2016']['amplitude']['L5_Pyr'] * normProb  # orig value: 0.05  
pmat['TCM']['CT6']      = data['TC_Ji2016']['amplitude']['L6_Pyr'] * normProb  # orig value: 0.05  

pmat['TCM']['NGF1']	    = data['TC_Ji2016']['amplitude']['L1'] * normProb  # orig value: 0.25
pmat['TCM']['PV2']	    = data['TC_Ji2016']['amplitude']['L23_PV'] * normProb  # orig value: 0.25
pmat['TCM']['SOM2']	    = data['TC_Ji2016']['amplitude']['L23_SOM'] * normProb # orig value: 0.25
pmat['TCM']['NGF2']	    = data['TC_Ji2016']['amplitude']['L23_SOM'] * normProb # orig value: 0.25
pmat['TCM']['PV3']	    = data['TC_Ji2016']['amplitude']['L23_PV'] * normProb # orig value: 0.25
pmat['TCM']['SOM3']	    = data['TC_Ji2016']['amplitude']['L23_SOM'] * normProb # orig value: 0.25
pmat['TCM']['NGF3']	    = data['TC_Ji2016']['amplitude']['L23_SOM'] * normProb # orig value: 0.25

pmat['TCM']['PV5A']	    = data['TC_Ji2016']['amplitude']['L5_PV'] * normProb  # orig value: 0.15  
pmat['TCM']['SOM5A']    = data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb  # orig value: 0.15  
pmat['TCM']['SOM5B']	= data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb # orig value: 0.15  
pmat['TCM']['PV5B']	    = data['TC_Ji2016']['amplitude']['L5_PV'] * normProb # orig value: 0.15  
pmat['TCM']['SOM5B']	= data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb  # orig value: 0.15  
pmat['TCM']['NGF5B']	= data['TC_Ji2016']['amplitude']['L5_SOM'] * normProb# orig value: 0.15  

pmat['TCM']['PV6']	    = data['TC_Ji2016']['amplitude']['L6_PV'] * normProb # orig value: 0.05  
pmat['TCM']['SOM6']	    = data['TC_Ji2016']['amplitude']['L6_SOM'] * normProb # orig value: 0.05  
pmat['TCM']['NGF6']	    = data['TC_Ji2016']['amplitude']['L6_SOM'] * normProb # orig value: 0.05  



# --------------------------------------------------
## Weights  (=unitary conn somatic PSP amplitude)

# don't have specific exp data so use 0.5mV, which is approx avg in Constantinople, 2016  
TCMweight = 0.5

wmat['TCM']['IT2']	    = TCMweight  # 0.6
wmat['TCM']['IT3']	    = TCMweight  # 0.6
wmat['TCM']['IT5A']	    = TCMweight  # 0.6  #* thalfctr
wmat['TCM']['IT5B']	    = TCMweight  # 0.6  #* thalfctr
wmat['TCM']['PT5B']	    = TCMweight  # 0.6  #* thalfctr
wmat['TCM']['IT6']	    = TCMweight  # 0.6  #* thalfctr
wmat['TCM']['CT5A']     = TCMweight  # 0.6  #* thalfctr
wmat['TCM']['CT5B']     = TCMweight  # 0.6  #* thalfctr
wmat['TCM']['CT6']      = TCMweight  # 0.6  #* thalfctr

wmat['TCM']['NGF1']	    = TCMweight  # 0.25
wmat['TCM']['PV2']	    = TCMweight  # 0.25
wmat['TCM']['SOM2']	    = TCMweight  # 0.25
wmat['TCM']['NGF2']	    = TCMweight  # 0.25
wmat['TCM']['PV3']	    = TCMweight  # 0.25
wmat['TCM']['SOM3']	    = TCMweight  # 0.25
wmat['TCM']['NGF3']	    = TCMweight  # 0.25

wmat['TCM']['PV5A']     = TCMweight  # 0.25 #* thalfctr
wmat['TCM']['SOM5A']    = TCMweight  # 0.25 #* thalfctr
wmat['TCM']['SOM5B']    = TCMweight  # 0.25 #* thalfctr
wmat['TCM']['PV5B']	    = TCMweight  # 0.25 #* thalfctr
wmat['TCM']['SOM5B']    = TCMweight  # 0.25 #* thalfctr
wmat['TCM']['NGF5B']    = TCMweight  # 0.25 #* thalfctr

wmat['TCM']['PV6']	    = TCMweight  # 0.25  #* thalfctr
wmat['TCM']['SOM6']	    = TCMweight  # 0.25  #* thalfctr
wmat['TCM']['NGF6']	    = TCMweight  # 0.25  #* thalfctr


# --------------------------------------------------
# Length constant (matrix more broadly tunes spatially than core)
lmat['TCM'] = {}
for post in wmat['TCM'].keys():
    lmat['TCM'][post] = 1000  # um - arbitrary value

# --------------------------------------------------
# Save data to pkl file
savePickle = 1

if savePickle:
    import pickle
    with open('conn.pkl', 'wb') as f:
        pickle.dump({'pmat': pmat, 'lmat': lmat, 'wmat': wmat, 'bins': bins, 'connDataSource': connDataSource}, f)