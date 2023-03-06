import numpy as np
from matplotlib import pyplot as plt
#%matplotlib inline
from scipy.signal import convolve
import random
from numpy import loadtxt
plt.rcParams['figure.figsize'] = [12, 7]
from scipy.fft import fft, ifft
import argparse
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

if(parserargs.Option == 1):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromAnaROOT(parserargs.Filepath)
elif(parserargs.Option == 2):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromROOT(parserargs.Filepath, fChSkip)
elif(parserargs.Option == 3):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromTxt(parserargs.Filepath)


fig, axs = plt.subplots(2, 1)
fig.subplots_adjust(left=0.075, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

print("Total processed events", len(EventID_V))

for ch in range(fNChannels):
    if(ch in fChSkip): continue
    #print(ch, len(WvMeanChDict[ch]), len(WvRMSChDict[ch]), len(EventID_V))
    axs[0].scatter(EventID_V, WvMeanChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
    
axs[0].set_xlabel("event ID"); axs[0].set_ylabel("Pedestal mean [ADC]"); 
axs[0].legend(loc='right')
axs[0].grid()


print(WvRMSChDict[0])
for ch in range(fNChannels):
    if(ch in fChSkip): continue
    axs[1].scatter(EventID_V, WvRMSChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
axs[1].set_xlabel("event ID"); axs[1].set_ylabel("Pedestal RMS [ADC]"); 
axs[1].legend(loc='right')
axs[1].grid()



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
axs[1][1].set_xlabel("Channel"); axs[1][1].set_ylabel("Pedestal RMS [ADC]"); 
axs[1][1].set_xticks(np.arange(fNChannels))

axs[1][2].errorbar(df.Ch, df.ChPed, df.ChPedErr, ls='none', marker="o")
axs[1][2].grid()
axs[1][2].set_xlabel("Channel"); axs[1][2].set_ylabel("Pedestal mean [ADC]");
axs[1][2].set_xticks(np.arange(fNChannels)) 

plt.show()

