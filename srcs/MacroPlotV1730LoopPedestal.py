import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = [12, 7]
import argparse
import glob
import pandas as pd
from PlotUtilsV1730 import *

# Macro to obtaing the correlation matrix

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

fChSkip = parserargs.ChSkip

InputFileList={}
for filepath in glob.iglob(parserargs.Filepath):
    print(filepath)
    setupLabel = ""
    for ch in range(fNChannels):
        if( filepath.find("ch"+str(ch)+".")!=-1 ): 
            InputFileList[ch]=filepath


#Data=[ [for j in range(fNChannels)] for i in range(len(InputFileList)) ]

DataRMS = np.empty((len(InputFileList), fNChannels))
DataPed = np.empty((len(InputFileList), fNChannels))

fRMSMin=1e6; fRMSMax=-1e6;
fPedMin=1e6; fPedMax=-1e6;

fileix=-1
SetupLabelV=[]
AveragePed=[]
AverageRMS=[]
AveragePedErr=[]
AverageRMSErr=[]
for ch in range(fNChannels):
    if( (ch in InputFileList)==False): continue
    filepath=InputFileList[ch]
    fileix+=1
    print(filepath)
    chToSkip = [ch]
    setupLabel = "S_in_"+str(ch)
    SetupLabelV.append(setupLabel)
    
    print("Skip Channels: ", chToSkip)
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromAnaROOT(filepath, chToSkip)

    df = pd.DataFrame(columns=['Ch','ChPed','ChPedErr','ChRMS','ChRMSErr'])

    for ch in range(fNChannels):
        ped_mean=ped_err=rms_mean=rms_err=0
        if(ch in chToSkip): 
            DataRMS[fileix][ch]=None
            DataPed[fileix][ch]=None
        else:
            ped_mean = np.mean(WvMeanChDict[ch])
            ped_err = np.std(WvMeanChDict[ch])
            rms_mean = np.mean(WvRMSChDict[ch])
            rms_err = np.std(WvRMSChDict[ch])
            DataRMS[fileix][ch]=rms_mean
            DataPed[fileix][ch]=ped_mean
            if(ped_mean>fPedMax): fPedMax=ped_mean
            if(ped_mean<fPedMin): fPedMin=ped_mean
            if(rms_mean>fRMSMax): fRMSMax=rms_mean
            if(rms_mean<fRMSMin): fRMSMin=rms_mean

        df.loc[ch] = [ch, ped_mean, ped_err, rms_mean, rms_err]

    print(df)

    AverageRMS.append(np.mean(df['ChRMS']))
    AveragePed.append(np.mean(df['ChPed']))
    AverageRMSErr.append(0)#np.std(df['ChRMS']))
    AveragePedErr.append(0)#np.std(df['ChPed']))
    



print(AverageRMS, AverageRMSErr, AveragePedErr)
fig2, axs = plt.subplots(2, 2)
fig2.subplots_adjust(left=0.075, bottom=0.082, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

PlotBaseline = axs[0,0].imshow(DataPed.T, vmin=fPedMin, vmax=fPedMax)
axs[0,0].set_xticks(np.arange(len(InputFileList)))
axs[0,0].set_xticklabels(SetupLabelV, rotation='vertical')
barBaseline = plt.colorbar(PlotBaseline, ax=axs[0,0])
barBaseline.ax.set_title('Pedestal [ADC]')
axs[0,0].set_ylabel('ChID')

PlotRMS = axs[1,0].imshow(DataRMS.T, vmin=fRMSMin, vmax=fRMSMax)
axs[1,0].set_xticks(np.arange(len(InputFileList)))
axs[1,0].set_xticklabels(SetupLabelV, rotation='vertical')
barRMS = plt.colorbar(PlotRMS, ax=axs[1,0])
barRMS.ax.set_title('RMS [ADC]')
axs[1,0].set_ylabel('ChID')


axs[0,1].errorbar(np.arange(len(InputFileList)), AveragePed, AveragePedErr, ls='none', marker="o")
axs[0,1].set_xticks(np.arange(len(InputFileList)))
axs[0,1].set_xticklabels(SetupLabelV, rotation='vertical')
axs[0,1].set_ylabel(r"$\langle $ Pedestal $\rangle$ [ADC]")

axs[1,1].errorbar(np.arange(len(InputFileList)), AverageRMS, AverageRMSErr, ls='none', marker="o")
axs[1,1].set_xticks(np.arange(len(InputFileList)))
axs[1,1].set_xticklabels(SetupLabelV, rotation='vertical')
axs[1,1].set_ylabel(r"$\langle $ RMS $\rangle$ [ADC]")

plt.show()

