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
#import datetime
from datetime import datetime
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec


fBaseline = 0
fNChannels = 16
fNOpDet = 320
fWfSize = 5000


##### COLOR LIST #####
#derived from https://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
def GetColorList(N):
    NUM_COLORS = N
    LINE_STYLES = ['solid']#, 'dashed', 'dashdot', 'dotted']
    NUM_STYLES = len(LINE_STYLES)

    ColorList = []
    cm = plt.get_cmap('gist_rainbow')
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
def ReadFromAnaROOT(filename, ch_skip=[], fBoardIDList=[], fVerbose=0, fMaxEvents=-1):

    print("Reading from AnaTree")
    file = ROOT.TFile.Open(filename)
    file.Print()
    tree_path="caenv1730dump/nt_wvfm"#"caenv1730ana/nt_wvfm;1"
    tuple =  file.Get(tree_path)

    tuple.Print();
    
    wvMeanChDict = {}
    wvRMSChDict = {}
    chTempDict = {}
    eventIDDict = {}
    timeStampDict = {}

    # initialize dictionaries
    #for ch in range(fNChannels):
    for ch in range(fNOpDet):
        eventIDDict[ch]=[]
        timeStampDict[ch]=[]
        wvMeanChDict[ch]=[]
        wvRMSChDict[ch]=[]
        chTempDict[ch]=[]

    eventCounter = 0
    for tree_entry in range( tuple.GetEntries() ):
        if(fMaxEvents!=-1 and eventCounter>fMaxEvents): continue
        eventCounter+=1
        tuple.GetEntry(tree_entry)
        if(fVerbose>=1): print(int(tuple.ch), tuple.ped, tuple.rms, tuple.temp)
        if(int(tuple.ch) in ch_skip): continue
        if( (tuple.boardId in fBoardIDList)==False and len(fBoardIDList)>0): continue
        wvMeanChDict[int(tuple.ch)].append(tuple.ped)
        wvRMSChDict[int(tuple.ch)].append(tuple.rms)
        chTempDict[int(tuple.ch)].append(tuple.temp)

        eventIDDict[int(tuple.ch)].append(tuple.art_ev)
        timeStampDict[int(tuple.ch)].append(tuple.stamp_time)

        if(fVerbose>=1): print(int(tuple.art_ev), int(tuple.caen_ev), "Ch:",int(tuple.ch), tuple.stamp_time)

    eventIDDictFinal = {}
    timeStampDictFinal = {}
    wvMeanChDictFinal = {}
    wvRMSChDictFinal = {}
    chTempDictFinal = {}
    for ch in range(fNOpDet):
        if(len(eventIDDict[ch])>0):
            eventIDDictFinal[ch] = eventIDDict[ch]
            timeStampDictFinal[ch] = timeStampDict[ch]
            wvMeanChDictFinal[ch] = wvMeanChDict[ch]
            wvRMSChDictFinal[ch] = wvRMSChDict[ch]
            chTempDictFinal[ch] = chTempDict[ch]
            print(ch, len(eventIDDictFinal[ch]))
        
    

    return eventIDDictFinal, timeStampDictFinal, wvMeanChDictFinal, wvRMSChDictFinal, chTempDictFinal
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
def PlotAverageBaseline(eventIDDict, timeStampDict, wvMeanChDict, wvRMSChDict, chTempDict, chSkip=[], useTimeStamp=True):
    timeAxisHandle = len(timeStampDict)>0
    fig = plt.figure("PMT V1730")
    fig.subplots_adjust(left=0.075, bottom=0.15, right=0.97, top=0.95, wspace=0.3, hspace=0.45)
    gs = GridSpec(3, 8, hspace=0, wspace=0)

    # Get channel-color map
    #cList=GetColorList(fNChannels)
    cList=GetColorList(len(wvMeanChDict.keys()))
    cDict = {}
    for ch_ix, ch in enumerate(wvMeanChDict.keys()):
        cDict[ch] = cList[ch_ix]

    ax2 = fig.add_subplot(gs[2,0:7])
    ax1 = fig.add_subplot(gs[1,0:7], sharex=ax2)
    ax0 = fig.add_subplot(gs[0,0:7], sharex=ax1)
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    axLeg = fig.add_subplot(gs[0:2,7], sharex=ax1)
    axLeg.axis("off")
   

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        print(ch)
        #print(ch, len(WvMeanChDict[ch]), len(WvRMSChDict[ch]), len(EventID_V))
        labelStr=""
        if(ch%2==0):
            labelStr="Ch"+str(ch)
        if(useTimeStamp):
            dateTimesV =  []
            for ts in timeStampDict[ch]:
                #print(ts, datetime.fromtimestamp(ts).strftime('%D:%H:%M'))
                dateTimesV.append( datetime.fromtimestamp(ts) ) 
            ax0.scatter(dateTimesV, chTempDict[ch], c=cDict[ch], label=labelStr, marker='o', s=3.)
            ax1.scatter(dateTimesV, wvMeanChDict[ch], c=cDict[ch], label=labelStr, marker='o', s=3.)
            ax2.scatter(dateTimesV, wvRMSChDict[ch], c=cDict[ch], label=labelStr, marker='o', s=3.)
        else:
            ax0.scatter(eventIDDict[ch], chTempDict[ch], c=cDict[ch], label=labelStr, marker='o', s=3.)
            ax1.scatter(eventIDDict[ch], wvMeanChDict[ch], c=cDict[ch], label=labelStr, marker='o', s=3.)
            ax2.scatter(eventIDDict[ch], wvRMSChDict[ch], c=cDict[ch], label=labelStr, marker='o', s=3.) 

    ax1.set_ylabel("Pedestal mean [ADC]"); 
    ax1.grid()

    ax2.set_ylabel("Pedestal RMS [ADC]"); 
    ax2.grid()
    
    ax0.grid()
    ax0.set_ylabel(r"Temperature [$^\circ$C]"); 

    if(useTimeStamp):
        fmt = mdates.DateFormatter('%D:%H:%M:%S')
        ax2.xaxis.set_major_formatter(fmt)
        ax2.set_xlabel("Time");
        ax2.tick_params(axis='x', labelrotation = 45)
    else:
        ax2.set_xlabel("event ID");

    handles, labels = ax1.get_legend_handles_labels()
    axLeg.legend(handles, labels)


    fig.savefig("average_pedestal.pdf")
##########################################



# plot baseline statistics
##########################################
def PlotStatisticsBaseline(eventID_V, wvMeanChDict, wvRMSChDict, chTempDict, chSkip=[]):

    fig2, axs = plt.subplots(2, 3)
    fig2.subplots_adjust(left=0.075, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)

    df = pd.DataFrame(columns=['Ch','ChPed','ChPedErr','ChRMS','ChRMSErr', 'ChTemp', 'ChTempErr', 'ChLoc'])
    dfString = pd.DataFrame(columns=["Ch", "Ped", "RMS", "T"])

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        ped_mean = np.mean(wvMeanChDict[ch])
        ped_err = np.std(wvMeanChDict[ch])
        rms_mean = np.mean(wvRMSChDict[ch])
        rms_err = np.std(wvRMSChDict[ch])
        t_mean = np.mean(chTempDict[ch])
        t_err = np.std(chTempDict[ch])
        ch_loc = ch%2

        df.loc[ch] = [ch, ped_mean, ped_err, rms_mean, rms_err, t_mean, t_err, ch_loc]

        dfString.loc[ch] = [
            str(ch),
            "{:.1f}".format(ped_mean)+r"$\pm$"+"{:.1f}".format(ped_err),
            "{:.2f}".format(rms_mean)+r"$\pm$"+"{:.2f}".format(rms_err),
            "{:.2f}".format(t_mean)+r"$\pm$"+"{:.2f}".format(t_err)
        ]


    print("Baseline stats data frame", df)
    df.to_csv("baseline_stats.csv")

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        labelStr="Ch"+str(ch)
        axs[0][0].hist(wvMeanChDict[ch], label=labelStr, histtype="step")
    axs[0][0].set_xlabel("Baseline mean"); axs[0][0].set_ylabel("# entries"); 
    axs[0][0].grid()
    axs[0][0].legend()

    binsRMS=np.arange(0, 5, 0.005)
    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        labelStr="Ch"+str(ch)
        axs[1][0].hist(wvRMSChDict[ch], bins=binsRMS, label=labelStr, histtype="step")
    axs[1][0].set_xlabel("Baseline RMS"); axs[1][0].set_ylabel("# entries"); 
    axs[1][0].grid()
    axs[1][0].legend()


    axs[1][1].errorbar(df.Ch, df.ChRMS, df.ChRMSErr, ls='none', marker="o")
    axs[1][1].grid()
    axs[1][1].set_xlabel("Channel"); axs[1][1].set_ylabel("Pedestal RMS [ADC]"); 
    axs[1][1].set_xticks(np.arange( min(df.Ch), max(df.Ch), 2))

    axs[0][1].errorbar(df.Ch, df.ChPed, df.ChPedErr, ls='none', marker="o")
    axs[0][1].grid()
    axs[0][1].set_xlabel("Channel"); axs[0][1].set_ylabel("Pedestal mean [ADC]");
    axs[0][1].set_xticks(np.arange( min(df.Ch), max(df.Ch), 2)) 


    axs[1][2].errorbar(df.ChTemp[df.ChLoc==0], df.ChRMS[df.ChLoc==0], df.ChRMSErr[df.ChLoc==0], ls='none', marker="o", label="Even")
    axs[1][2].errorbar(df.ChTemp[df.ChLoc==1], df.ChRMS[df.ChLoc==1], df.ChRMSErr[df.ChLoc==1], ls='none', marker="o", label="Odd")
    axs[1][2].grid()
    axs[1][2].set_xlabel(r"Channel temperature [$^\circ$ C]");
    axs[1][2].set_ylabel("Pedestal RMS [ADC]"); 
    axs[1][2].legend()

    axs[0][2].errorbar(df.ChTemp[df.ChLoc==0], df.ChPed[df.ChLoc==0], df.ChPedErr[df.ChLoc==0], ls='none', marker="o", label="Even")
    axs[0][2].errorbar(df.ChTemp[df.ChLoc==1], df.ChPed[df.ChLoc==1], df.ChPedErr[df.ChLoc==1], ls='none', marker="o", label="Odd")
    axs[0][2].grid()
    axs[0][2].set_xlabel(r"Channel temperature [$^\circ$ C]");
    axs[0][2].set_ylabel("Pedestal mean [ADC]");
    axs[0][2].legend()

    fig2.savefig("statistics_pedestal.pdf")


    """fig3, axs = plt.subplots(1, 1)
    axs.table(cellText=dfString.values, colLabels=dfString.columns, loc='center')
    axs.axis('off')
    fig3.savefig("statistics_pedestal_table.pdf")"""
##########################################

