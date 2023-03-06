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
#plt.set_cmap("YlOrBr")
#fRMSMax = 2.65
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





    

"""
fig2, axs = plt.subplots(2, 3)
fig2.subplots_adjust(left=0.075, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

ChPed={}
ChPedErr={}
ChRMS={}
ChRMSErr={}


df = pd.DataFrame(columns=['Ch','ChPed','ChPedErr','ChRMS','ChRMSErr'], index=np.arange(0, fNChannels, 1))
dfString = pd.DataFrame(columns=["Ch", "Ped", "RMS"], index=np.arange(0, fNChannels, 1))

for ch in range(fNChannels):
    if(ch in fChSkip): continue
    ped_mean = np.mean(WvMeanChDict[ch])
    ped_err = np.std(WvMeanChDict[ch])
    rms_mean = np.mean(WvRMSChDict[ch])
    rms_err = np.std(WvRMSChDict[ch])

    df.loc[ch] = [ch,ped_mean, ped_err, rms_mean, rms_err]

    dfString.loc[ch] = [
        str(ch),
        "{:.1f}".format(ped_mean)+r"$\pm$"+"{:.1f}".format(ped_err),
        "{:.2f}".format(rms_mean)+r"$\pm$"+"{:.2f}".format(rms_err)
    ]


print(df)

for ch in range(fNChannels):
    if(ch in fChSkip): continue
    axs[0][0].hist(WvMeanChDict[ch], label="Ch"+str(ch), histtype="step")
axs[0][0].set_xlabel("Baseline mean"); axs[0][0].set_ylabel("# entries"); 
axs[0][0].grid()
axs[0][0].legend()

binsRMS=np.arange(0, 5, 0.005)
for ch in range(fNChannels):
    if(ch in fChSkip): continue
    axs[0][1].hist(WvRMSChDict[ch], bins=binsRMS, label="Ch"+str(ch), histtype="step")
axs[0][1].set_xlabel("Baseline RMS"); axs[0][1].set_ylabel("# entries"); 
axs[0][1].grid()
axs[0][1].legend()


axs[1][0].table(cellText=dfString.values, colLabels=dfString.columns, loc='center')
axs[1][0].axis('off')

axs[1][1].errorbar(df.Ch, df.ChRMS, df.ChRMSErr, ls='none', marker="o")
axs[1][1].grid()
axs[1][1].set_xlabel("Channel"); axs[1][1].set_ylabel("Pedestal RMS"); 

axs[1][2].errorbar(df.Ch, df.ChPed, df.ChPedErr, ls='none', marker="o")
axs[1][2].grid()
axs[1][2].set_xlabel("Channel"); axs[1][2].set_ylabel("Pedestal Mean"); 

plt.show()

"""