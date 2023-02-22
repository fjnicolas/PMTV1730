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

params = {'legend.fontsize': 'small',
          'figure.figsize': (12, 8),
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium'}
plt.rcParams.update(params)


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--Filepath", help="Input file name",  default="none")
parser.add_argument("-fft", "--FFT", help="Plot FFT", type=int, default=0)
parser.add_argument("-n", "--NEv", help="Max Events", type=int, default=2)
parserargs = parser.parse_args()
fBaseline = 0

print("FFFF", parserargs.Filepath, parserargs.FFT)

##########################################
file = ROOT.TFile.Open(parserargs.Filepath)
tree =  file.Get("caenv1730dump/events")
print ("Tree Entries: ", tree.GetEntries())
NTreeEntries=tree.GetEntries(); NTreeEntries=75
##########################################

evCounter=0
for tree_entry in range( NTreeEntries ):
    tree.GetEntry(tree_entry)
    evCounter+=1
    if(evCounter>parserargs.NEv): continue
    #Event IDs
    eventID= tree.fEvent
    runID= tree.fRun
    #Waveforms
    #fTicksVec=tree.fTicksVec
    fWvfmsVec=tree.fWvfmsVec

    fig, axs = plt.subplots(4, 4, num="Event ID: "+str(eventID))
    fig.subplots_adjust(left=0.05, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)
    
    print("Event ID: ", eventID)
    for ch, wf in enumerate(fWvfmsVec):
        print(ch, len(wf))
        wf=np.array(wf)-fBaseline
        if(parserargs.FFT==1):
            print("Making FFT")
            wf_fft = fft(np.array(wf))
            axs[ch//4, ch%4].plot(np.abs(wf_fft))
            axs[ch//4, ch%4].set_yscale("log")
    
        else:
            ch_stddev = np.std(wf)
            labName = "StdDev="+"{:.1f}".format(ch_stddev)+" ADC"
            axs[ch//4, ch%4].plot(wf, label=labName)
            axs[ch//4, ch%4].legend()

        axs[ch//4, ch%4].set_title("Ch="+str(ch))
        axs[ch//4, ch%4].set_xlabel("Time Tick [2 ns]")
        axs[ch//4, ch%4].set_ylabel("[ADC]")
        axs[ch//4, ch%4].grid()

        

    plt.show()
wfsize=48
