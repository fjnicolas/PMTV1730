import numpy as np
from matplotlib import pyplot as plt
import argparse
import pandas as pd
from PlotUtilsV1730 import *

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


params = {'figure.figsize': (10, 8),
        'axes.labelsize': 'x-large',   
        'axes.titlesize': 'x-large',      
        'xtick.labelsize':'large',
        'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

fig, axs = plt.subplots(2, 2)
fig.subplots_adjust(left=0.13, bottom=0.082, right=0.99, top=0.95, wspace=0.5, hspace=0.45)
fig2, axs2 = plt.subplots(2, 2)
fig2.subplots_adjust(left=0.13, bottom=0.082, right=0.99, top=0.95, wspace=0.5, hspace=0.45)
fig3, axs3 = plt.subplots(2, 2)
fig3.subplots_adjust(left=0.13, bottom=0.082, right=0.99, top=0.95, wspace=0.5, hspace=0.45)


fMS = ['o', 's', 'D', '^', 'v', '>', '<', 'p', '*', 'h', 'H', 'x']
fLS = ['-', '--', '-.', ':'] * len(fMS)
fColor1 = "midnightblue"
fColor2 = "maroon"
fCList = ["midnightblue", "maroon", "forestgreen", "slategrey", "darkmagenta",]
fCapSize = 5.
fMarkerSize = 5.
fAlpha=0.8
fLineWidth=2.

fBoardNumberDict = {160:"S3",176:"S5",192:"S7",208:"S9",
                    224:"S11",240:"S13",256:"S16",272:"S18"}

fChList = list(range(160, 290))

fOption=5
fFiles = parserargs.Filepath
fLabels = parserargs.Label
if(fOption==1):   
    """fFiles.append("Comparisons/FreqScan/freqscna500mHz/baseline_stats.csv")
    fLabels.append("0.5 Hz")"""
    """fFiles.append("Comparisons/FreqScan/freqscna1Hz/baseline_stats.csv")
    fLabels.append("1 Hz")"""
    fFiles.append("Comparisons/FreqScan/freqscna2_5Hz/baseline_stats.csv")
    fLabels.append("2.5 Hz")
    fFiles.append("Comparisons/FreqScan/freqscna5Hz/baseline_stats.csv")
    fLabels.append("5 Hz")
    fFiles.append("Comparisons/FreqScan/freqscna10Hz/baseline_stats.csv")
    fLabels.append("10 Hz")
elif(fOption==2):
    fFiles.append("Comparisons/Run070523Fan/baseline_stats_nofan.csv"); fLabels.append("woFan")
    fFiles.append("Comparisons/Run070523Fan/baseline_stats_fan.csv"); fLabels.append("wFan")
    #fFiles.append("Comparisons/Run070523Fan/baseline_stats_fanoff.csv"); fLabels.append("wFanOff")
elif(fOption==3):
    fFiles.append("Comparisons/Signals/baseline_stats_signal.csv"); fLabels.append("wSignal")
    fFiles.append("Comparisons/Signals/baseline_stats_nosignal.csv"); fLabels.append("woSignal")
    fBoard1ChList = [160, 162, 164, 166, 168, 170, 172, 174]
    fBoard2ChList = list(range(176, 183))
    fBoard3ChList = list(range(192, 208))
    fChList = fBoard1ChList + fBoard2ChList +fBoard3ChList
elif(fOption==4):
    fFiles.append("Comparisons/QuasiEqualization/baseline_stats_eq.csv"); fLabels.append("Setup 1")
    fFiles.append("Comparisons/QuasiEqualization/baseline_stats_quasi.csv"); fLabels.append("Setup 2")
elif(fOption==5):
    fFiles.append("Comparisons/AfterEqua/baseline_stats_preequa.csv"); fLabels.append("PreEq")
    fFiles.append("Comparisons/AfterEqua/baseline_stats_afterequa.csv"); fLabels.append("PostEq")


def AddSlotNumbers(ax, minCh, maxCh):
    ax.set_xticks(np.arange( minCh, maxCh, 16 ))
    for index, xmaj in enumerate(ax.xaxis.get_majorticklocs()):
        if(index<len(ax.xaxis.get_majorticklocs())):
            ax.axvline(x=xmaj, ls="--", linewidth = 2., color = 'dimgrey')
            ax.text(xmaj, 1.05, fBoardNumberDict[xmaj], ha='center', va='center', transform=ax.get_xaxis_transform(), color= 'dimgrey', rotation="45")

dfV = []
for ix, path in enumerate(fFiles):
    print(path)
    df = pd.read_csv(path) 
    Label = fLabels[ix]
    
    chFilter = df["Ch"].astype("int").isin(fChList)
    df = df[chFilter]
    dfV.append( df )

    eb1=axs[0][0].errorbar(df.Ch, df.ChRMS, df.ChRMSErr, ls="none", marker=fMS[ix], label=Label, capsize=fCapSize, alpha=fAlpha, c=fCList[ix], ms=fMarkerSize)
    eb2=axs[0][1].errorbar(df.Ch, df.ChPed, df.ChPedErr, ls="none", marker=fMS[ix], label=Label, capsize=fCapSize, alpha=fAlpha, c=fCList[ix], ms=fMarkerSize)
    eb3=axs[1][0].errorbar(df.Ch, df.ChTemp, df.ChTempErr, ls="none", marker=fMS[ix], label=Label, capsize=fCapSize, alpha=fAlpha, c=fCList[ix], ms=fMarkerSize)
    eb1[-1][0].set_linestyle(fLS[ix])
    eb2[-1][0].set_linestyle(fLS[ix])
    eb3[-1][0].set_linestyle(fLS[ix])
    
    statsL = " (" + GetStatsStr(df.ChRMS, entries=False, prec="{:.2f}") + ")"
    binsRMS=np.linspace( min(df.ChRMS), max(df.ChRMS)+0.6*(max(df.ChRMS)-min(df.ChRMS)), 15)
    axs3[0][0].hist(df.ChRMS, bins=binsRMS, histtype="step", label=Label+statsL, ls=fLS[ix], color=fCList[ix], lw=fLineWidth)

    statsL = " ("+GetStatsStr(df.ChPed, entries=False, prec="{:.0f}")+")"
    binsPed=np.linspace( min(df.ChPed), max(df.ChPed)+0.6*(max(df.ChPed)-min(df.ChPed)), 15)
    axs3[0][1].hist(df.ChPed, bins=binsPed, histtype="step", label=Label+statsL, ls=fLS[ix], color=fCList[ix],lw=fLineWidth)

    statsL = " ("+GetStatsStr(df.ChTemp, entries=False, prec="{:.0f}")+")"
    binsTemp=np.linspace( min(df.ChTemp), max(df.ChTemp)+0.6*(max(df.ChTemp)-min(df.ChTemp)), 15)
    axs3[1][0].hist(df.ChTemp, binsTemp, histtype="step", label=Label+statsL, ls=fLS[ix], color=fCList[ix],lw=fLineWidth)

    if(ix!=0):

        df1 = dfV[0]
        df2 = dfV[ix]
        l1=fLabels[0]
        l2=fLabels[ix]

        sigmaLabel = r"$\overline{\ \sigma}_{\rm B_{\rm"+l2+r"}}/\overline{\ \sigma}_{\rm B_{\rm"+l1+r"}}$"
        baselineLabel = r"$\overline{\rm B}_{\rm"+l2+r"}/\overline{\rm B}_{\rm"+l1+r"}$"
        axs2[0][0].plot(df1.Ch, df2.ChRMS/df1.ChRMS, ls="none", marker=fMS[ix], c=fCList[ix], label=sigmaLabel, ms=fMarkerSize)
        axs2[0][1].plot(df1.Ch, df2.ChPed/df1.ChPed, ls="none", marker=fMS[ix], c=fCList[ix], label=baselineLabel, ms=fMarkerSize)
        tempLabel = r"$\overline{\rm T}_{\rm"+l2+r"}/\overline{\rm T}_{\rm"+l1+r"}$"
        axs2[1][0].plot(df1.Ch, df2.ChTemp/df1.ChTemp, ls="none", marker=fMS[ix], c=fCList[ix], label=tempLabel, ms=fMarkerSize)

        
        

minCh = min(dfV[0].Ch)
maxCh = max(dfV[0].Ch)

# figure 1 axis labels
##---## RMS
axs[0][0].grid()
axs[0][0].legend()
axs[0][0].set_xlabel("Channel"); 
axs[0][0].set_ylabel(r"$\overline{\sigma}_{\rm B}$ [ADC]"); 
#axs[0][0].set_ylim(2.4, 2.8)
AddSlotNumbers(axs[0,0], minCh, maxCh)
##---## Baseline
axs[0][1].grid()
axs[0][1].legend()
axs[0][1].set_xlabel("Channel")
axs[0][1].set_ylabel(r"$\overline{\rm B}$ [ADC]")
AddSlotNumbers(axs[0,1], minCh, maxCh)
##---## Temperature
axs[1][0].grid()
axs[1][0].legend()
axs[1][0].set_xlabel("Channel"); 
axs[1][0].set_ylabel(r"$\overline{\rm T}$ [$^\circ$ C]")
axs[1][0].set_xticks(np.arange( min(dfV[0].Ch), max(dfV[0].Ch), 16 ))
AddSlotNumbers(axs[1,0], minCh, maxCh)

# figure 2 axis labels
##---## RMS
axs2[0][0].grid()
axs2[0][0].legend()
axs2[0][0].set_xlabel("Channel");
axs2[0][0].set_ylabel( r"$\overline{\ \sigma}_{\rm B}$ ratio" )
AddSlotNumbers(axs2[0,0], minCh, maxCh)
##---## Baseline
axs2[0][1].grid()
axs2[0][1].legend()
axs2[0][1].set_xlabel("Channel"); 
axs2[0][1].set_ylabel(r"$\overline{\rm B}$ ratio")
AddSlotNumbers(axs2[0,1], minCh, maxCh)
##---## Temperature
axs2[1][0].grid()
axs2[1][0].legend()
axs2[1][0].set_xlabel("Channel"); 
axs2[1][0].set_ylabel(r"$\overline{\rm T}$ ratio")
AddSlotNumbers(axs2[1,0], minCh, maxCh)


# figure 3 axis labels
##---## RMS
axs3[0][0].legend()
axs3[0][0].set_ylabel("# entries");
axs3[0][0].set_xlabel( r"$\overline{\ \sigma}_{\rm B}$" )
##---## Baseline
axs3[0][1].legend()
axs3[0][1].set_ylabel("# entries"); 
axs3[0][1].set_xlabel(r"$\overline{\rm B}$")
##---## Temperature
axs3[1][0].legend()
axs3[1][0].set_ylabel("# entries"); 
axs3[1][0].set_xlabel(r"$\overline{\rm T}$")
        


outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/data_comparison/"
# create directory if it doesn't exist
if not os.path.exists(outputFilepath):
    os.makedirs(outputFilepath)
os.makedirs(outputFilepath+"png/", exist_ok=True)
os.makedirs(outputFilepath+"pdf/", exist_ok=True)

SaveSubplots(fig, axs, outputFilepath, "comp_baseline" )
SaveSubplots(fig2, axs2, outputFilepath, "comp_ratio" )
SaveSubplots(fig3, axs3, outputFilepath, "comp_1dhist" )
plt.show()