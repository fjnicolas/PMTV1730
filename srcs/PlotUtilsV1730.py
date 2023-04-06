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
import datetime
import matplotlib.dates as mdates

fBaseline = 0
fNChannels = 16
fWfSize = 5000


##### COLOR LIST #####
#derived from https://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
def GetColorList(N):
    NUM_COLORS = N
    LINE_STYLES = ['solid']#, 'dashed', 'dashdot', 'dotted']
    NUM_STYLES = len(LINE_STYLES)

    ColorList = []
    cm = plt.get_cmap('gist_ncar')
    for i in range(NUM_COLORS):
        ColorList.append(cm(i//NUM_STYLES*float(NUM_STYLES)/NUM_COLORS))
        #lines[0].set_linestyle(LINE_STYLES[i%NUM_STYLES])

    return ColorList

#### Read from ROOT txt file ####
##########################################
def ReadFromROOT(filename, ch_skip=[], waveformRange=None, maxEvents=1e6):
    file = ROOT.TFile.Open(filename)
    tree =  file.Get("caenv1730dump/events")
    print ("Tree Entries: ", tree.GetEntries())
    
    startIx = 0
    endIx = fWfSize
    if(waveformRange):
        startIx = waveformRange[0]
        endIx = waveformRange[1]

    ## Initialize dictionaries
    eventCounter = 0
    wvMeanChDict = {}
    wvRMSChDict = {}
    eventID_V = []
    for ch in range(fNChannels):
        wvMeanChDict[ch]=[]
        wvRMSChDict[ch]=[]

    for entry,tree_entry in enumerate(tree):
        if(eventCounter>maxEvents): continue
        eventCounter+=1
        
        #Event IDs
        eventID= tree_entry.fEvent
        runID= tree_entry.fRun
        print(eventCounter, "Run ID", runID, "Event ID: ", eventID)
        
        #Waveforms
        #fTicksVec=tree.fTicksVec
        WvfmsVec=list(tree_entry.fWvfmsVec)

        
        for ch, wave in enumerate(WvfmsVec):
            if(ch in ch_skip): continue
            wf=list(wave)
            if(endIx>len(wf)):
                endIx = len(wf)
            wf=np.array(wf[startIx:endIx])

            ch_mean = np.mean(wf)   
            ch_stddev = np.std(wf)
        
            wvMeanChDict[ch].append(ch_mean)
            wvRMSChDict[ch].append(ch_stddev)

        eventID_V.append(eventID)

    return eventID_V, wvMeanChDict, wvRMSChDict
##########################################


#### Read from AnaROOT txt file ####
##########################################
def ReadFromAnaROOT(filename, ch_skip=[], fVerbose=0):

    print("Reading from AnaTree")
    file = ROOT.TFile.Open(filename)
    tuple =  file.Get("caenv1730ana/nt_wvfm;1")
    
    wvMeanChDict = {}
    wvRMSChDict = {}
    eventIDDict = {}
    for ch in range(fNChannels):
        wvMeanChDict[ch]=[]
        wvRMSChDict[ch]=[]
        eventIDDict[ch]=[]

    for tree_entry in range( tuple.GetEntries() ):
        tuple.GetEntry(tree_entry)
        if(fVerbose>=1): print(int(tuple.ch), tuple.ped, tuple.rms)
        if(int(tuple.ch) in ch_skip): continue
        wvMeanChDict[int(tuple.ch)].append(tuple.ped)
        wvRMSChDict[int(tuple.ch)].append(tuple.rms)
        eventIDDict[int(tuple.ch)].append(tuple.art_ev)
        if(fVerbose>=1): print(int(tuple.art_ev), int(tuple.caen_ev))

    print( list(eventIDDict.keys())[0] )
    eventID_V = eventIDDict[list(eventIDDict.keys())[0]]

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

##########################################
# class handle for baseline vs time plots
class TimeXAxisHandle:
    def __init__(self, startTime, frequency):
        self.startTime = startTime
        self.frequency = frequency
        self.period = 1./self.frequency
##########################################

##########################################
# plot baseline mean and RMS as a function of time
# input is EventID_V, WvMeanChDict, WvRMSDict
def PlotAverageBaseline(eventID_V, wvMeanChDict, wvRMSChDict, chSkip=[], timeAxisHandle=None):
    fig, axs = plt.subplots(2, 1)
    fig.subplots_adjust(left=0.075, bottom=0.14, right=0.97, top=0.95, wspace=0.3, hspace=0.45)

    print("Total processed events", len(eventID_V))

    dateTimesV =  []
    if(timeAxisHandle):
        print("Have handled events")
        for evID in eventID_V:
            dateTimesV.append(timeAxisHandle.startTime+datetime.timedelta(seconds=timeAxisHandle.period*evID ) )
    print(dateTimesV)

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        #print(ch, len(WvMeanChDict[ch]), len(WvRMSChDict[ch]), len(EventID_V))
        if(timeAxisHandle):
            axs[0].scatter(dateTimesV, wvMeanChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
        else:
            axs[0].scatter(eventID_V, wvMeanChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
        
    if(timeAxisHandle):
        fmt = mdates.DateFormatter('%D:%H:%M')
        axs[0].xaxis.set_major_formatter(fmt)
        axs[0].set_xlabel("Time");
        axs[0].tick_params(axis='x', labelrotation = 45)
    else:
        axs[0].set_xlabel("event ID");
    
    axs[0].set_ylabel("Pedestal mean [ADC]"); 
    axs[0].legend(loc='right')
    axs[0].grid()

    for ch in range(fNChannels):
        if(ch in chSkip): continue
        if(timeAxisHandle):
            axs[1].scatter(dateTimesV, wvRMSChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
        else:
            axs[1].scatter(eventID_V, wvRMSChDict[ch], label="Ch"+str(ch), marker='o', s=3.)
    if(timeAxisHandle):
        fmt = mdates.DateFormatter('%D:%H:%M')
        axs[1].xaxis.set_major_formatter(fmt)
        axs[1].set_xlabel("Time"); 
        axs[1].tick_params(axis='x', labelrotation = 45)
    else:
        axs[1].set_xlabel("event ID"); 
    
    axs[1].set_ylabel("Pedestal RMS [ADC]"); 
    axs[1].legend(loc='right')
    axs[1].grid()

    fig.savefig("average_pedestal.pdf")
##########################################



# plot baseline statistics
##########################################
def PlotStatisticsBaseline(eventID_V, wvMeanChDict, wvRMSChDict, chSkip=[]):
    fig2, axs = plt.subplots(2, 3)
    fig2.subplots_adjust(left=0.075, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

    ChPed={}
    ChPedErr={}
    ChRMS={}
    ChRMSErr={}


    df = pd.DataFrame(columns=['Ch','ChPed','ChPedErr','ChRMS','ChRMSErr'], index=np.arange(0, fNChannels, 1))
    dfString = pd.DataFrame(columns=["Ch", "Ped", "RMS"], index=np.arange(0, fNChannels, 1))

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        ped_mean = np.mean(wvMeanChDict[ch])
        ped_err = np.std(wvMeanChDict[ch])
        rms_mean = np.mean(wvRMSChDict[ch])
        rms_err = np.std(wvRMSChDict[ch])

        df.loc[ch] = [ch,ped_mean, ped_err, rms_mean, rms_err]

        dfString.loc[ch] = [
            str(ch),
            "{:.1f}".format(ped_mean)+r"$\pm$"+"{:.1f}".format(ped_err),
            "{:.2f}".format(rms_mean)+r"$\pm$"+"{:.2f}".format(rms_err)
        ]


    print("Baseline stats data frame", df)
    df.to_csv("baseline_stats.csv")

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        axs[0][0].hist(wvMeanChDict[ch], label="Ch"+str(ch), histtype="step")
    axs[0][0].set_xlabel("Baseline mean"); axs[0][0].set_ylabel("# entries"); 
    axs[0][0].grid()
    axs[0][0].legend()

    binsRMS=np.arange(0, 5, 0.005)
    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        axs[0][1].hist(wvRMSChDict[ch], bins=binsRMS, label="Ch"+str(ch), histtype="step")
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

    fig2.savefig("statistics_pedestal.pdf")
##########################################

