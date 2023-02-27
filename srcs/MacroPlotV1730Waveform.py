import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = [12, 7]
from scipy.fft import fft, ifft
import argparse
import ROOT
import pandas as pd

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
parserargs = parser.parse_args()
fBaseline = 0

fX_min = 0
fX_max = 5000

fWfSize = 5000



##########################################
def PlotBoardChannels(dat, eventID, fFFT=False):

    fig, axs = plt.subplots(4, 4, num="Event ID: "+str(eventID))
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.95, wspace=0.3, hspace=0.45)
    print(dat)
    
    for (ch, wf) in dat.iteritems():

        if(ch==0): continue
        chIx=ch-1
        wf = np.array(wf.values)-fBaseline
        print('Plotting channel : ', ch, " Length : ", len(wf))

        if(fFFT==True):
            print("Making FFT")
            wf_fft = fft(np.array(wf))
            axs[chIx//4, chIx%4].plot(np.abs(wf_fft))
            axs[chIx//4, chIx%4].set_yscale("log")
            

        else:
            chIx_stddev = np.std(wf)
            labName = "StdDev="+"{:.1f}".format(chIx_stddev)+" ADC"
            axs[chIx//4, chIx%4].plot(wf, label=labName)
            axs[chIx//4, chIx%4].legend()
            axs[chIx//4, chIx%4].set_xlim(fX_min, fX_max)

        axs[chIx//4, chIx%4].set_title("Ch="+str(chIx))
        axs[chIx//4, chIx%4].set_xlabel("Time Tick [2 ns]")
        axs[chIx//4, chIx%4].set_ylabel("[ADC]")
        axs[chIx//4, chIx%4].grid()


    plt.show()
##########################################


##########################################
def GetDataFromROOT(filepath):
    file = ROOT.TFile.Open(parserargs.Filepath)
    tree =  file.Get("caenv1730dump/events")
    print ("Tree Entries: ", tree.GetEntries())
   
    evCounter=0
    for tree_entry in range( tree.GetEntries() ):
        tree.GetEntry(tree_entry)
        #Event IDs
        eventID= tree.fEvent
        evCounter+=1
        if(evCounter>parserargs.NEv): continue
        print(evCounter, "Event ID: ", eventID)
        #Waveforms
        #fTicksVec=tree.fTicksVec
        fWvfmsVec=tree.fWvfmsVec

        data = pd.DataFrame()
        for ch, wf in enumerate(fWvfmsVec):
            wf=np.array(wf)
            data[ch+1] = wf

        print(data)
        
        PlotBoardChannels(data, eventID, parserargs.FFT)
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
    GetDataFromROOT(parserargs.Filepath)
elif(parserargs.Option==2):
    GetDataFromTxt(parserargs.Filepath)