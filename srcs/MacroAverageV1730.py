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
parserargs = parser.parse_args()
fBaseline = 0



##########################################
file = ROOT.TFile.Open(parserargs.Filepath)
tree =  file.Get("caenv1730dump/events")
print ("Tree Entries: ", tree.GetEntries())
##########################################

eventCounter = 0

fNChannels = 16
WvMeanChDict = {}
WvRMSChDict = {}
EventID_V = []

for ch in range(fNChannels):
    WvMeanChDict[ch]=[]
    WvRMSChDict[ch]=[]


for tree_entry in range( tree.GetEntries() ):
    tree.GetEntry(tree_entry)
    eventCounter+=1
    #Event IDs
    eventID= tree.fEvent
    runID= tree.fRun
    print(eventCounter, "Event ID: ", eventID)
    #Waveforms
    #fTicksVec=tree.fTicksVec
    fWvfmsVec=tree.fWvfmsVec


    for ch, wf in enumerate(fWvfmsVec):
        wf=np.array(wf)

        ch_mean = np.mean(wf)   
        ch_stddev = np.std(wf)
     
        WvMeanChDict[ch].append(ch_mean)
        WvRMSChDict[ch].append(ch_stddev)

    EventID_V.append(eventID)
   


        


fig, axs = plt.subplots(2, 1)
fig.subplots_adjust(left=0.05, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

print("Total processed events", eventCounter)

for ch in range(fNChannels):
    axs[0].scatter(EventID_V, WvMeanChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
axs[0].set_xlabel("event ID"); axs[0].set_ylabel("Pedestal mean [ADC]"); 
axs[0].legend()
axs[0].grid()


print(WvRMSChDict[0])
for ch in range(fNChannels):
    axs[1].scatter(EventID_V, WvRMSChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
axs[1].set_xlabel("event ID"); axs[1].set_ylabel("Pedestal RMS [ADC]"); 
axs[1].legend()
axs[1].grid()

plt.show()

