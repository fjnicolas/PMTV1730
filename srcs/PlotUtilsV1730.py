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


fBaseline = 0
fNChannels = 16
fWfSize = 5000

#### Read from ROOT txt file ####
##########################################
def ReadFromROOT(filename, ch_skip=[]):
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
        if(eventCounter>parserargs.NEv): continue
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
            if(ch in ch_skip): continue
            wf=np.array(wf)

            ch_mean = np.mean(wf)   
            ch_stddev = np.std(wf)
        
            wvMeanChDict[ch].append(ch_mean)
            wvRMSChDict[ch].append(ch_stddev)

        eventID_V.append(eventID)

    return eventID_V, wvMeanChDict, wvRMSChDict
##########################################


#### Read from AnaROOT txt file ####
##########################################
def ReadFromAnaROOT(filename, ch_skip=[]):
    from ROOT import gROOT

    file = ROOT.TFile.Open(filename)
    tuple =  file.Get("caenv1730ana/nt_wvfm;1")
    #names = [b.GetName() for b in tuple.GetListOfBranches()]
    #print(names)
    ## Initialize dictionaries
    
    eventCounter = 0
    wvMeanChDict = {}
    wvRMSChDict = {}
    eventID_V = set()
    for ch in range(fNChannels):
        wvMeanChDict[ch]=[]
        wvRMSChDict[ch]=[]

    for tree_entry in range( tuple.GetEntries() ):
        tuple.GetEntry(tree_entry)
        #print(int(tuple.ch), tuple.ped, tuple.rms)
        if(int(tuple.ch) in ch_skip): continue
        wvMeanChDict[int(tuple.ch)].append(tuple.ped)
        wvRMSChDict[int(tuple.ch)].append(tuple.rms)
        eventID_V.add(int(tuple.art_ev))

    return list(eventID_V), wvMeanChDict, wvRMSChDict
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

        if(ixStep>parserargs.NEv): continue

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
