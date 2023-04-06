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
parser.add_argument("-tHandle", "--TimeHandle", help="Use time handle", type=int, default=0)
parser.add_argument("-chSkip", "--ChSkip", type=int, action='append', help="Channels to skip", default=[])
parserargs = parser.parse_args()


# compute pedestal in two zones
fPreBins = [0, 800]
fPostBins = [1500, -1]
fSelectedBins=[0, 5000]
fSelectedBins = fPostBins

fChSkip = parserargs.ChSkip

if(parserargs.Option == 1):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromAnaROOT(parserargs.Filepath)
elif(parserargs.Option == 2):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromROOT(parserargs.Filepath, fChSkip, waveformRange=fSelectedBins, maxEvents=parserargs.NEv)
elif(parserargs.Option == 3):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromTxt(parserargs.Filepath)



print("MAX", max(EventID_V), "MIN", min(EventID_V) )
print(EventID_V)



if(parserargs.TimeHandle==1):
    # run 03/23
    #timeHandle = TimeXAxisHandle( datetime.datetime(2023, 3, 23, 19, 42, 0), 0.05)
    # run 03/24
    timeHandle = TimeXAxisHandle( datetime.datetime(2023, 3, 24, 16, 56, 0), 0.05)
    # run 02/28
    #timeHandle = TimeXAxisHandle( datetime.datetime(2023, 2, 28, 17, 18, 0), 0.04)
    
    PlotAverageBaseline(EventID_V, WvMeanChDict, WvRMSChDict, timeAxisHandle=timeHandle)
else:
    PlotAverageBaseline(EventID_V, WvMeanChDict, WvRMSChDict)



PlotStatisticsBaseline(EventID_V, WvMeanChDict, WvRMSChDict)


plt.show()