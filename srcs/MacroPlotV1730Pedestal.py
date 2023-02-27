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
parser.add_argument("-o", "--Option", help="Input option", type=int, default=1)
parserargs = parser.parse_args()
fBaseline = 0
fNChannels = 16
fWfSize = 5000



#### Read from ROOT txt file ####
##########################################
def ReadFromROOT(filename):
    file = ROOT.TFile.Open(filename)
    tree =  file.Get("caenv1730dump/events")
    print ("Tree Entries: ", tree.GetEntries())

    ## Initialize dictionaries
    eventCounter = 0
    wvMeanChDict = {}
    wvRMSChDict = {}
    eventID_V = []
    for ch in range(fNChannels):
        wvMeanChDict[ch]=[]
        wvRMSChDict[ch]=[]

    for tree_entry in range( tree.GetEntries() ):
        tree.GetEntry(tree_entry)
        eventCounter+=1
        #Event IDs
        eventID= tree.fEvent
        runID= tree.fRun
        print(eventCounter, "Event ID: ", eventID)
        #Waveforms
        #fTicksVec=tree.fTicksVec
        WvfmsVec=tree.fWvfmsVec


        for ch, wf in enumerate(WvfmsVec):
            wf=np.array(wf)

            ch_mean = np.mean(wf)   
            ch_stddev = np.std(wf)
        
            wvMeanChDict[ch].append(ch_mean)
            wvRMSChDict[ch].append(ch_stddev)

        eventID_V.append(eventID)

    return eventID_V, wvMeanChDict, wvRMSChDict
##########################################


#### Read from ROOT file ####
##########################################
def ReadFromTxt(filename):

    ##########################################
    import pandas as pd
    DF = pd.read_table(parserargs.Filepath, delim_whitespace=True, header=None)
    ##########################################

    ## Initialize dictionaries
    eventCounter = 0
    wvMeanChDict = {}
    wvRMSChDict = {}
    eventID_V = []
    for ch in range(fNChannels):
        wvMeanChDict[ch]=[]
        wvRMSChDict[ch]=[]

    for ixStep in range( 0, len(DF), fWfSize ):
        print(ixStep, ixStep+fWfSize)
        data = DF[ixStep:ixStep+fWfSize]

        eventCounter+=1
        #Event IDs
        eventID= eventCounter
        print(eventCounter, "Event ID: ", eventID)

        for (ch, wf) in data.iteritems():
            if(ch==0): continue
            chIx=ch-1
            wf = np.array(wf.values)-fBaseline
            print('Plotting channel : ', ch, " Length : ", len(wf))

            ch_mean = np.mean(wf)   
            ch_stddev = np.std(wf)
        
            wvMeanChDict[chIx].append(ch_mean)
            wvRMSChDict[chIx].append(ch_stddev)

        eventID_V.append(eventID)

    return eventID_V, wvMeanChDict, wvRMSChDict
##########################################

if(parserargs.Option == 1):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromROOT(parserargs.Filepath)
elif(parserargs.Option == 2):
    EventID_V, WvMeanChDict, WvRMSChDict = ReadFromTxt(parserargs.Filepath)




fig, axs = plt.subplots(2, 1)
fig.subplots_adjust(left=0.075, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

print("Total processed events", len(EventID_V))

for ch in range(fNChannels):
    axs[0].scatter(EventID_V, WvMeanChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
axs[0].set_xlabel("event ID"); axs[0].set_ylabel("Pedestal mean [ADC]"); 
axs[0].legend(loc='right')
axs[0].grid()


print(WvRMSChDict[0])
for ch in range(fNChannels):
    axs[1].scatter(EventID_V, WvRMSChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
axs[1].set_xlabel("event ID"); axs[1].set_ylabel("Pedestal RMS [ADC]"); 
axs[1].legend(loc='right')
axs[1].grid()

plt.show()

