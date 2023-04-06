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

# execute as 
#py /Users/franciscojaviernicolas/Work/PMTV1730/PMTV1730/srcs/MacroComparePedestalStats.py -s basline_stats_prebins.csv -s basline_stats_postbins.csv -l "Noise only" -l "Pulse all channels"
#py /Users/franciscojaviernicolas/Work/PMTV1730/PMTV1730/srcs/MacroComparePedestalStats.py -s ../run_022823/basline_stats.csv -s basline_stats_postbins.csv -l "Noise only" -l "All channels (post pulse)"

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--Filepath", type=str, action='append', help="Input dataframes", default=[])
parser.add_argument("-l", "--Label", type=str, action='append', help="Input label", default=[])
parserargs = parser.parse_args()

print(parserargs.Filepath)
if(len(parserargs.Label)!=len(parserargs.Filepath)):
    print("Warning: no labels specified")
    parserargs.Label = ["" for _ in range(len(parserargs.Filepath))]

fig, axs = plt.subplots(2, 2)
fig.subplots_adjust(left=0.075, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

fLS = ["solid", "dashed"]
fMS = ["o", "^"]


dfV = []
for ix, path in enumerate(parserargs.Filepath):
    print(path)
    df = pd.read_csv(path) 
    dfV.append( df )
    Label = parserargs.Label[ix]

    eb1=axs[0][0].errorbar(df.Ch, df.ChRMS, df.ChRMSErr, ls="none", marker=fMS[ix], label=Label, capsize=10, alpha=0.9)
    eb2=axs[0][1].errorbar(df.Ch, df.ChPed, df.ChPedErr, ls="none", marker=fMS[ix], label=Label, capsize=10, alpha=0.9)
    eb1[-1][0].set_linestyle(fLS[ix])
    eb2[-1][0].set_linestyle(fLS[ix])


fNChannels  = 16
axs[0][0].grid()
axs[0][0].legend()
axs[0][0].set_xlabel("Channel"); 
axs[0][0].set_ylabel("Pedestal RMS [ADC]"); 
axs[0][0].set_xticks(np.arange(fNChannels))
axs[0][0].set_ylim(2.4, 2.8)

axs[0][1].grid()
axs[0][1].legend()
axs[0][1].set_xlabel("Channel"); 
axs[0][1].set_ylabel("Pedestal mean [ADC]");
axs[0][1].set_xticks(np.arange(fNChannels)) 

fig.savefig("ComparisonPedestal.pdf")

plt.show()