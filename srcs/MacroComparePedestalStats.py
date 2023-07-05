import numpy as np
from matplotlib import pyplot as plt
#%matplotlib inline
plt.rcParams['figure.figsize'] = [12, 7]
import argparse
import pandas as pd


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

fMS = ['o', 's', 'D', '^', 'v', '>', '<', 'p', '*', 'h', 'H', 'x']
fLS = ['-', '--', '-.', ':'] * len(fMS)

fFiles = parserargs.Filepath
fLabels = parserargs.Label

"""fFiles.append("V1730/freqscna10Hz/baseline_stats.csv")
fLabels.append("10 Hz")
fFiles.append("V1730/freqscna5Hz/baseline_stats.csv")
fLabels.append("5 Hz")
fFiles.append("V1730/freqscna1Hz/baseline_stats.csv")
fLabels.append("1 Hz")
fFiles.append("V1730/freqscna500mHz/baseline_stats.csv")
fLabels.append("0.5 Hz")"""

fFiles.append("Comparisons/Run070523Fan/baseline_stats_nofan.csv")
fLabels.append("woFan")
fFiles.append("Comparisons/Run070523Fan/baseline_stats_fan.csv")
fLabels.append("wFan")


dfV = []
for ix, path in enumerate(fFiles):
    print(path)
    df = pd.read_csv(path) 
    dfV.append( df )
    Label = fLabels[ix]

    eb1=axs[0][0].errorbar(df.Ch, df.ChRMS, df.ChRMSErr, ls="none", marker=fMS[ix], label=Label, capsize=10, alpha=0.9)
    eb2=axs[0][1].errorbar(df.Ch, df.ChPed, df.ChPedErr, ls="none", marker=fMS[ix], label=Label, capsize=10, alpha=0.9)
    eb3=axs[1][0].errorbar(df.Ch, df.ChTemp, df.ChTempErr, ls="none", marker=fMS[ix], label=Label, capsize=10, alpha=0.9)
    eb1[-1][0].set_linestyle(fLS[ix])
    eb2[-1][0].set_linestyle(fLS[ix])


axs[0][0].grid()
axs[0][0].legend()
axs[0][0].set_xlabel("Channel"); 
axs[0][0].set_ylabel(r"$\sigma_{\rm B}$ [ADC]"); 
axs[0][0].set_ylim(2.4, 2.8)

axs[0][1].grid()
axs[0][1].legend()
axs[0][1].set_xlabel("Channel"); 
axs[0][1].set_ylabel("B [ADC]");

axs[1][0].grid()
axs[1][0].legend()
axs[1][0].set_xlabel("Channel"); 
axs[1][0].set_ylabel(r"T [$^\circ$ C]")
axs[1][0].set_xticks(np.arange( min(dfV[0].Ch), max(dfV[0].Ch), 16 ))
for index, xmaj in enumerate(axs[1][0].xaxis.get_majorticklocs()):
    axs[1][0].axvline(x=xmaj, ls="--", linewidth = 2., color = 'dimgrey')

fig.savefig("ComparisonPedestal.pdf")

if(len(fFiles)==2):
    fig2, axs = plt.subplots(2, 2)
    fig2.subplots_adjust(left=0.075, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)


    df1 = pd.read_csv(fFiles[0]) 
    df2 = pd.read_csv(fFiles[1]) 

    axs[0][0].plot(df1.Ch, np.array(df2.ChRMS)/np.array(df1.ChRMS), ls="none", marker=fMS[0])
    axs[0][1].plot(df1.Ch, np.array(df2.ChPed)/np.array(df1.ChPed), ls="none", marker=fMS[0])
    axs[1][0].plot(df1.Ch, df2.ChTemp/df1.ChTemp, ls="none", marker=fMS[0])

    axs[0][0].grid()
    axs[0][0].legend()
    axs[0][0].set_xlabel("Channel"); 
    axs[0][0].set_ylabel(r"$\sigma_{\rm B_1}/\sigma_{\rm B_2}$");

    axs[0][1].grid()
    axs[0][1].legend()
    axs[0][1].set_xlabel("Channel"); 
    axs[0][1].set_ylabel(r"${\rm B}_1/{\rm B}_2$");

    axs[1][0].grid()
    axs[1][0].legend()
    axs[1][0].set_xlabel("Channel"); 
    axs[1][0].set_ylabel(r"${\rm T}_1/{\rm T}_2$")
    axs[1][0].set_xticks(np.arange( min(df1.Ch), max(df1.Ch), 16 ))
    for index, xmaj in enumerate(axs[1][0].xaxis.get_majorticklocs()):
        axs[1][0].axvline(x=xmaj, ls="--", linewidth = 2.0, color = 'dimgrey')
plt.show()