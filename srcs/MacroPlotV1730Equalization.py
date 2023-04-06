import numpy as np
from matplotlib import pyplot as plt
#%matplotlib inline
from scipy.signal import convolve
import random
from numpy import loadtxt
plt.rcParams['figure.figsize'] = [12, 7]
from scipy.fft import fft, ifft
import argparse
import glob
import ROOT
import pandas as pd
from pandas.plotting import table
from PlotUtilsV1730 import *

params = {'legend.fontsize': 'small',
          'figure.figsize': (12, 8),
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium'}
plt.rcParams.update(params)


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--Filepath", help="Input file name",  default="none")
parser.add_argument("-o", "--Option", help="Input option", type=int, default=1)
parser.add_argument("-n", "--NEv", help="Max Events", type=int, default=1e6)
parser.add_argument("-chSkip", "--ChSkip", type=int, action='append', help="Channels to skip", default=[])
parserargs = parser.parse_args()
fBaseline = 0
fNChannels = 16
fWfSize = 5000

fTargetADC = 8000
fChListWildcard = [0, 4, 8, 12]

InputFileList={}
for filepath in glob.iglob(parserargs.Filepath):
    print("File", filepath)
    if( filepath.find(".root")!=-1 ):
        ix=filepath.find(".root")-1
        setupLabel = filepath[ix]
        print("Setup label:", setupLabel)
        InputFileList[setupLabel] = filepath
   
   



fRMSMin=1e6; fRMSMax=-1e6;
fPedMin=1e6; fPedMax=-1e6;

fileix=-1
SetupLabelV=[]
file_counter = -1
DataDict = {}

for fileix, setup in enumerate(InputFileList):
    filepath=InputFileList[setup]
    file_counter+=1
    print(ix, " Analyzing file", filepath, " with setup", setup)

    setupLabel = setup
    SetupLabelV.append(setupLabel)
    
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromAnaROOT(filepath)

    df = pd.DataFrame(columns=['Ch','ChPed','ChPedErr','ChRMS','ChRMSErr'])

    for ch in range(fNChannels):
        # initialize
        ped_mean=ped_err=rms_mean=rms_err=0

        # fill
        ped_mean = np.mean(WvMeanChDict[ch])
        ped_err = np.std(WvMeanChDict[ch])
        rms_mean = np.mean(WvRMSChDict[ch])
        rms_err = np.std(WvRMSChDict[ch])
       
        if(ped_mean>fPedMax): fPedMax=ped_mean
        if(ped_mean<fPedMin): fPedMin=ped_mean
        if(rms_mean>fRMSMax): fRMSMax=rms_mean
        if(rms_mean<fRMSMin): fRMSMin=rms_mean

        df.loc[ch] = [ch, ped_mean, ped_err, rms_mean, rms_err]
    
    DataDict[setupLabel] = df

    print(df)


DCOffset_V = np.array( [32768, 30768, 34768, 36768] )

PedMean = {}
PedMeanErr = {}
for ch in range(0, fNChannels):
    PedMean[ch] = []
    PedMeanErr[ch] = []

for fileix, dflabel in enumerate(DataDict):
    df = DataDict[dflabel]
    print("Dataaa\n", df)

    for ch in range(0, fNChannels):
        PedMean[ch].append(df["ChPed"][ch])
        PedMeanErr[ch].append(df["ChPedErr"][ch])


#plt.set_cmap("YlOrBr")
#fRMSMax = 2.65
fig2, axs = plt.subplots(2, 2)
fig2.subplots_adjust(left=0.075, bottom=0.082, right=0.99, top=0.95, wspace=0.3, hspace=0.45)
fColorList = GetColorList(fNChannels)

FitSlopes = []
FitIntercepts = []
FitSlopesErr = []
FitInterceptsErr = []

AdHocDCOffsets = {}

for ch in range(0, fNChannels):
    # get linear fit paramters
    pars, cov = np.polyfit(DCOffset_V, np.array(PedMean[ch]), 1, cov=True)
    pars_err = np.sqrt(np.diag(cov))
    FitSlopes.append( pars[0] )
    FitIntercepts.append( pars[1] )
    FitSlopesErr.append( pars_err[0] )
    FitInterceptsErr.append( pars_err[1] )

    adHocDC = AdHocDCOffsets[ch] = (fTargetADC-pars[1])/pars[0]

    equalizedFhiclLabel = "daq.fragment_receiver.channelPedestal"+str(ch)+": "+str(int(adHocDC))
    print(equalizedFhiclLabel)

    # plot
    xfit_plot = np.linspace( min(DCOffset_V), max(DCOffset_V), 50)
    axs[0,0].errorbar(DCOffset_V, PedMean[ch], PedMeanErr[ch], ls='none', marker="o", c=fColorList[ch], label="Ch "+str(ch))
    axs[0,0].plot(xfit_plot, pars[0]*xfit_plot+pars[1], c=fColorList[ch])

    if(ch in fChListWildcard):
        axs[0,1].errorbar(DCOffset_V, PedMean[ch], PedMeanErr[ch], ls='none', marker="o", c=fColorList[ch], label="Ch "+str(ch))
        axs[0,1].plot(xfit_plot, pars[0]*xfit_plot+pars[1], c=fColorList[ch])

    

axs[0,0].legend()
axs[0,0].set_xlabel(r"DCOffset [$\tilde{\rm ADC}$]")
axs[0,0].set_ylabel("Pedestal mean [ADC]")
axs[0,0].grid()

axs[0,1].legend()
axs[0,1].set_xlabel(r"DCOffset [$\tilde{\rm ADC}$]")
axs[0,1].set_ylabel("Pedestal mean [ADC]")
axs[0,1].grid()

XChannelsV = np.arange(0, fNChannels, 1)
axs[1,0].errorbar(XChannelsV, FitSlopes, FitSlopesErr, ls='none', marker="o")
axs[1,1].errorbar(XChannelsV, FitIntercepts, FitInterceptsErr, ls='none', marker="o")

axs[1,0].legend()
axs[1,0].set_xlabel("Channel ID")
axs[1,0].set_ylabel(r"Fit slope [ADC/$\tilde{\rm ADC}$]")
axs[1,0].grid()


axs[1,1].legend()
axs[1,1].set_xlabel("Channel ID")
axs[1,1].set_ylabel("Fit intercept [ADC]")
axs[1,1].grid()

print("AdHocs DC offsets:\n", AdHocDCOffsets)


plt.show()