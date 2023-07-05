import numpy as np
from matplotlib import pyplot as plt
from scipy.fft import fft
import argparse
import ROOT
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
def PlotPowerSpectrum(dat):

    fig, axs = plt.subplots(nrows=4, ncols=4, figsize=(9, 6), sharex=True, sharey=True, num="Power Spectrum")
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.99, top=0.95, wspace=0.04, hspace=0.04)
    print("Data\n", dat)

    minCh = min(dat.keys())
    for ch in dat:
        chIx = ch - minCh
        wfFFT = dat[ch]
        
        print('Plotting channel : ', ch, " Length : ", len(wfFFT))
        ax =axs [chIx//4, chIx%4]
        
        
        ax.plot(fBinsFrquencies, np.abs(wf_fft)[0:int(fWfSize/2)])
        ax.set_yscale("log")
        ax.legend(title="Ch="+str(ch))
    
    fig.text(0.5, 0.03, 'Frequency [GHz]', ha='center', fontsize=14)
    fig.text(0.01, 0.5, 'Power (AU)', va='center', rotation='vertical', fontsize=14)
    
    outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/powerspectrum/"
    if not os.path.exists(outputFilepath):
        os.makedirs(outputFilepath)
    fig.savefig(outputFilepath+"waveform.pdf")
    plt.show()
##########################################



file = ROOT.TFile.Open(parserargs.Filepath)
tree =  file.Get("caenv1730dump/events")
print ("Tree Entries: ", tree.GetEntries())


# dictionaty to store average FFT
dataFFT = {}
numItersFFT = {}

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
   
    if(len(parserargs.BoardID)>0 and (not fBoardId in parserargs.BoardID)):     
        continue


    for ix, wf in enumerate(fWvfmsVec):
        wf = np.array(wf)
        wf_fft = np.abs( fft(wf) )
        ch=fChVec[ix]
        if(ch in numItersFFT):
            numItersFFT[ch]+=1
            dataFFT[ch]+=wf_fft
        else:
            numItersFFT[ch]=1
            dataFFT[ch]=wf_fft

for ch in numItersFFT:
    print(ch, numItersFFT[ch])
    dataFFT[ch]=dataFFT[ch]/numItersFFT[ch]


        


PlotPowerSpectrum(dataFFT)


