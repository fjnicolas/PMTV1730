import numpy as np
from matplotlib import pyplot as plt
from scipy.fft import fft
import argparse
import ROOT
import pandas as pd
import matplotlib.gridspec as gridspec
from PlotUtilsV1730 import *


# run as
# py MacroPlotV1730Waveform.py -s data/processed_data_result.root -o 1 -n 10
# if using ROOT output from CaenV1730Dump or
# py MacroPlotV1730Waveform.py -s data/PlotData.txt -o 2 -n 10
# if using CAEN wavedump txt output

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
parser.add_argument("-n", "--NEv", help="Max Events", type=int, default=1)
parser.add_argument("-o", "--Option", help="Input option", type=int, default=1)
parser.add_argument("-b", "--BoardID", type=int, action='append', help="Boards ID to analyze", default=[])
parserargs = parser.parse_args()

fX_min = 0
fX_max = 5000
fWfSize = 5000
fSamplingTime=2 #in ns
fDeltaF=1./(fWfSize*fSamplingTime)
fFMax=1./(fSamplingTime*2)
fBinsFrquencies=np.arange(0, fFMax, fDeltaF)

##########################################
def PlotBoardChannels(dat, eventID, fFFT=False):

    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(9, 6), sharex=True, sharey=True, num="Event ID: "+str(eventID))
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.99, top=0.95, wspace=0.04, hspace=0.04)
    print(dat)
    minCh = min(dat.keys())
    for (ch, wf) in dat.iteritems():
        
        if(ch==0): continue
        chIx=ch-minCh
        wf = np.array(wf.values)-fBaseline
        print('Plotting channel : ', ch, " Length : ", len(wf))
        ax =axs [chIx//4, chIx%4]
        if(fFFT==True):
            #subtract baseline
            wf=wf-np.mean(wf)
            wf_fft = fft(wf)
            ax.plot(fBinsFrquencies, np.abs(wf_fft)[0:int(fWfSize/2)])
            ax.set_yscale("log")
           
        
        else:
            chIx_stddev = np.std(wf)
            labName = "RMS="+"{:.1f}".format(chIx_stddev)+" ADC"
            ax.plot(wf, label=labName)
            ax.set_xlim(fX_min, fX_max)

        ax.legend(title="Ch="+str(ch))
    if(fFFT==True):
        fig.text(0.5, 0.03, 'Frequency [GHz]', ha='center', fontsize=14)
        fig.text(0.01, 0.5, 'Power (AU)', va='center', rotation='vertical', fontsize=14)
    else:
        fig.text(0.5, 0.03, 'Time Tick [2 ns]', ha='center', fontsize=14)
        fig.text(0.01, 0.5, '[ADC]', va='center', rotation='vertical', fontsize=14)

    outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/waveforms/"
    if not os.path.exists(outputFilepath):
        os.makedirs(outputFilepath)
    fig.savefig(outputFilepath+"waveform.pdf")
    plt.show()
##########################################


##########################################
def GetDataFromROOT(filepath, boardIdV = []):
    file = ROOT.TFile.Open(parserargs.Filepath)
    tree =  file.Get("caenv1730dump/events")
    print ("Tree Entries: ", tree.GetEntries())
   
    evCounter=0
    for tree_entry in range( tree.GetEntries() ):
        tree.GetEntry(tree_entry)
        #Event IDs
        eventID= tree.fEvent
        if(evCounter>=parserargs.NEv): continue
        print(evCounter, "Event ID: ", eventID)
        #Waveforms
        fWvfmsVec = tree.fWvfmsVec
        fChVec = tree.fWvfmsChVec
        fBoardId = tree.boardID
        if(len(boardIdV)>0 and (not fBoardId in boardIdV)): 
            print("BB", boardIdV, fBoardId)
            print(fBoardId in boardIdV)
            
            continue

        data = pd.DataFrame()
        for ix, wf in enumerate(fWvfmsVec):
            wf=np.array(wf)
            ch=fChVec[ix]
            data[ch] = wf

        print(data)
        
        PlotBoardChannels(data, eventID, parserargs.FFT)
        evCounter+=1
##########################################


##########################################
def GetDataFromTxt(filepath):
    DF = pd.read_table(filepath, delim_whitespace=True, header=None)
   
    for ixStep in range( 0, len(DF), fWfSize ):
        print(ixStep, ixStep+fWfSize)
        if(ixStep>parserargs.NEv): continue
        data = DF[ixStep:ixStep+fWfSize]
        PlotBoardChannels(data, ixStep, parserargs.FFT)
##########################################



if(parserargs.Option==1):
    GetDataFromROOT(parserargs.Filepath, parserargs.BoardID)
elif(parserargs.Option==2):
    GetDataFromTxt(parserargs.Filepath)
