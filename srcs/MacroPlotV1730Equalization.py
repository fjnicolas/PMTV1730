import numpy as np
from matplotlib import pyplot as plt
import argparse
import pandas as pd
from PlotUtilsV1730 import *


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--Filepath", help="Input file name",  default="none")
parser.add_argument("-o", "--Option", help="Input option", type=int, default=1)
parser.add_argument("-n", "--NEv", help="Max Events", type=int, default=1e6)
parser.add_argument("-timeStamp", "--UseTimeStamp", help="Use time stamp", type=int, default=1)
parser.add_argument("-chSkip", "--ChSkip", type=int, action='append', help="Channels to skip", default=[])
parser.add_argument("-b", "--BoardID", type=int, action='append', help="Boards ID to analyze", default=[])
parser.add_argument("-debug", "--Debug", type=int, default=0)
parserargs = parser.parse_args()

fColor1 = "midnightblue"
fColor2 = "maroon"

# Setup 1
fDACStart = 0.8; fDACStep = 0.02
# Setup 2
fDACStart = 0.83; fDACStep = 0.02
# Setup 3
#fDACStart = 0.78; fDACStep = 0.04
# Setup 4
#fDACStart = 0.85; fDACStep = 0.01

fMaxDAC=2**16 - 1
fDACMaxCali = fMaxDAC-0.96*fMaxDAC
fTargetADC = 0.95*2**14
fChListWildcard = [160, 190, 200, 210]
fChListWildcard = [193, 205, 207]
print("MaxDAC for calibration", fDACMaxCali)

# read the TTree
ChTempDict = {}
EventIDDict = {} 
RunIDDict = {} 
TimeStamp_V = {} 
WvMeanChDict = {}
WvRMSChDict = {}
RunIDDict, EventIDDict, TimeStampDict, WvMeanChDict, WvRMSChDict, ChTempDict = ReadFromAnaROOT(parserargs.Filepath, "all", parserargs.ChSkip, parserargs.BoardID, parserargs.Debug, parserargs.NEv)

# get run numbers
RunSetV = [item for sublist in list(RunIDDict.values()) for item in sublist]
RunSetV = sorted (set( RunSetV ))

# get channel IDs
ChSetV = sorted (set( RunIDDict.keys() ))

# initialize variables
DataDict = {}
DCOffset_V = []
PedMean = {}
PedMeanErr = {}
for ch in ChSetV:
    PedMean[ch] = []
    PedMeanErr[ch] = []
runCounter = 0
fRMSMin=1e6; fRMSMax=-1e6;
fPedMin=1e6; fPedMax=-1e6;

# loop over different runs
for runID in RunSetV:
    print(runID)
    dacValue = fMaxDAC - (fDACStart + fDACStep*runCounter) * fMaxDAC
    runCounter+=1

    df = pd.DataFrame(columns=['Ch','ChPed','ChPedErr','ChRMS','ChRMSErr'])

    for ch in ChSetV:
        # initialize
        ped_mean=ped_err=rms_mean=rms_err=0

        # fill
        runV = RunIDDict[ch]
        pedV = WvMeanChDict[ch]
        rmsV = WvRMSChDict[ch]
        filterRun = (runV == runID)
        ped_mean = np.mean( pedV[filterRun] )
        ped_err = np.std( pedV[filterRun] )
        rms_mean = np.mean( rmsV[filterRun] )
        rms_err = np.std( rmsV[filterRun] )
       
        if(ped_mean>fPedMax): fPedMax=ped_mean
        if(ped_mean<fPedMin): fPedMin=ped_mean
        if(rms_mean>fRMSMax): fRMSMax=rms_mean
        if(rms_mean<fRMSMin): fRMSMin=rms_mean

        df.loc[ch] = [ch, ped_mean, ped_err, rms_mean, rms_err]
       
        if( dacValue>=fDACMaxCali ):
            PedMean[ch].append(ped_mean)
            PedMeanErr[ch].append(ped_err)
    
    DataDict[runID] = df
    if( dacValue>=fDACMaxCali ):
        DCOffset_V.append( dacValue )

print("Run data", DataDict)


params = {'figure.figsize': (10, 8),
        'axes.labelsize': 'x-large',   
        'axes.titlesize': 'x-large',      
        'xtick.labelsize':'large',
        'ytick.labelsize':'x-large'}

plt.rcParams.update(params)

fig2, axs = plt.subplots(2, 2, constrained_layout=True)
fig2.subplots_adjust(left=0.13, bottom=0.082, right=0.99, top=0.95, wspace=0.5, hspace=0.45)

# get color list
cList=GetColorList(len(ChSetV))
cDict = {}
for ch_ix, ch in enumerate(ChSetV):
    cDict[ch] = cList[ch_ix]

minCh = min(ChSetV)
maxCh = max(ChSetV)

FitSlopes = []
FitIntercepts = []
FitSlopesErr = []
FitInterceptsErr = []
AdHocDCOffsets = {}

print("Target ADC is ", fTargetADC)
for ch in ChSetV:
    # get linear fit paramters
    print("PO", DCOffset_V )
    print("PO",  np.array(PedMean[ch]) )
    pars, cov = np.polyfit(DCOffset_V, np.array(PedMean[ch]), 1, cov=True)
    pars_err = np.sqrt(np.diag(cov))
    FitSlopes.append( pars[0] )
    FitIntercepts.append( pars[1] )
    FitSlopesErr.append( pars_err[0] )
    FitInterceptsErr.append( pars_err[1] )

    adHocDC = (fTargetADC-pars[1])/pars[0]

    print(ch, adHocDC)


    AdHocDCOffsets[ch] = adHocDC

    equalizedFhiclLabel = "daq.fragment_receiver.channelPedestal"+str(ch)+": "+str(int(adHocDC))
    print(equalizedFhiclLabel)

    # plot
    xfit_plot = np.linspace( min(DCOffset_V), max(DCOffset_V), 50)
    axs[0,0].errorbar(DCOffset_V, PedMean[ch], PedMeanErr[ch], ls='none', marker="o", c=cDict[ch], label="Ch "+str(ch))
    axs[0,0].plot(xfit_plot, pars[0]*xfit_plot+pars[1], c=cDict[ch])

    if(ch in fChListWildcard):
        axs[0,1].errorbar(DCOffset_V, PedMean[ch], PedMeanErr[ch], ls='none', marker="o", c=cDict[ch], label="Ch "+str(ch))
        axs[0,1].plot(xfit_plot, pars[0]*xfit_plot+pars[1], c=cDict[ch])

axs[0,0].set_xlabel(r"${\rm DAC}_{\rm DCoffset}$ ")
axs[0,0].set_ylabel(r"$\overline{\rm B}$ [ADC]")
axs[0,0].grid()

divider = make_axes_locatable(axs[0][0])
cax = divider.new_vertical(size = '2.5%', pad = 0.25)
fig2.add_axes(cax)
myCmap = LinearSegmentedColormap.from_list("ChannelColorMap", list(cDict.values()), N=len(cDict.values()))
cbar = plt.colorbar(cm.ScalarMappable(cmap=myCmap), cax=cax, orientation = 'horizontal', ticks=np.linspace(0, 1, 4))
cbar.ax.set_xticklabels(np.linspace(minCh, maxCh-1, 4).astype(int).tolist(), fontsize=10)
cbar.ax.set_ylabel("Ch", fontsize=10)

axs[0,1].legend()
axs[0,1].set_xlabel(r"${\rm DAC}_{\rm DCoffset}$ ")
axs[0,1].set_ylabel(r"$\overline{\rm B}$ [ADC]")
axs[0,1].grid()


axs[1,0].errorbar(ChSetV, FitSlopes, FitSlopesErr, ls='none', marker="o", c=fColor1)
axs[1,1].errorbar(ChSetV, FitIntercepts, FitInterceptsErr, ls='none', marker="o", c=fColor2)

axs[1,0].set_xlabel("Channel")
axs[1,0].set_ylabel(r"$\alpha$ [ADC]")
axs[1,0].grid()

axs[1,1].set_xlabel("Channel")
axs[1,1].set_ylabel(r"$\beta$ [ADC]")
axs[1,1].grid()


for ch in ChSetV:
    adHocDC = AdHocDCOffsets[ch]
    equalizedFhiclLabel = "daq.fragment_receiver.channelPedestal"+str(ch)+": "+str(int(adHocDC))
    print(equalizedFhiclLabel)

for ch in ChSetV:
    if(ch%16==0): print("\n\n")
    adHocDC = AdHocDCOffsets[ch]
    equalizedFhiclLabel = "daq.fragment_receiver.channelPedestal"+str(ch%16)+": "+str(int(adHocDC))
    print(equalizedFhiclLabel)


outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/equalization/"
# create directory if it doesn't exist
if not os.path.exists(outputFilepath):
    os.makedirs(outputFilepath)
os.makedirs(outputFilepath+"png/", exist_ok=True)
os.makedirs(outputFilepath+"pdf/", exist_ok=True)

SaveSubplots(fig2, axs, outputFilepath, "equalization" )




plt.show()

