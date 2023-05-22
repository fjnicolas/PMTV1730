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
parser.add_argument("-timeStamp", "--UseTimeStamp", help="Use time stamp", type=int, default=1)
parser.add_argument("-chSkip", "--ChSkip", type=int, action='append', help="Channels to skip", default=[])
parser.add_argument("-b", "--BoardID", type=int, action='append', help="Boards ID to analyze", default=[])
parser.add_argument("-debug", "--Debug", type=int, default=0)
parser.add_argument("-stats", "--MakeStats", type=int, default=0)
parserargs = parser.parse_args()


# compute pedestal in two zones
fPreBins = [0, 800]
fPostBins = [1500, -1]
fSelectedBins=[0, 5000]
fSelectedBins = fPostBins

fChSkip = parserargs.ChSkip

ChTempDict = {}
EventID_V = {} 
TimeStamp_V = {} 
WvMeanChDict = {}
WvRMSChDict = {}

if(parserargs.Option == 1):
    EventIDDict, TimeStampDict, WvMeanChDict, WvRMSChDict, ChTempDict = ReadFromAnaROOT(parserargs.Filepath, parserargs.ChSkip, parserargs.BoardID, parserargs.Debug, parserargs.NEv)
elif(parserargs.Option == 2):
    EventIDDict, WvMeanChDict, WvRMSChDict = ReadFromROOT(parserargs.Filepath, fChSkip, waveformRange=fSelectedBins, maxEvents=parserargs.NEv)
elif(parserargs.Option == 3):
    EventIDDict, WvMeanChDict, WvRMSChDict = ReadFromTxt(parserargs.Filepath)



PlotAverageBaseline(EventIDDict, TimeStampDict, WvMeanChDict, WvRMSChDict,  ChTempDict, parserargs.ChSkip, parserargs.UseTimeStamp)
if(parserargs.MakeStats==1):
    PlotStatisticsBaseline(EventIDDict, WvMeanChDict, WvRMSChDict, ChTempDict)


plt.show()