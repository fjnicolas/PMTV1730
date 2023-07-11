from matplotlib import pyplot as plt
import argparse
from PlotUtilsV1730 import *
from InputReaderV1730 import *

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--Filepath", help="Input file name",  default="none")
parser.add_argument("-r", "--Run", help="Select run number", type=int, default=-1)
parser.add_argument("-o", "--Option", help="Input option", type=int, default=1)
parser.add_argument("-n", "--NEv", help="Max Events", type=int, default=1e6)
parser.add_argument("-timeStamp", "--UseTimeStamp", help="Use time stamp", type=int, default=1)
parser.add_argument("-chSkip", "--ChSkip", type=int, action='append', help="Channels to skip", default=[])
parser.add_argument("-ch", "--Ch", type=int, help="Select channel", default=-1)
parser.add_argument("-b", "--BoardID", type=int, action='append', help="Boards ID to analyze", default=[])
parser.add_argument("-m", "--Mode", type=int, help="Read pedestal from all/start/end of the waveform", default=0)
parser.add_argument("-debug", "--Debug", type=int, default=0)
parser.add_argument("-pedOnly", "--PedOnly", type=int, default=0)
parser.add_argument("-stats", "--MakeStats", type=int, default=1)
parser.add_argument("-time", "--PlotTime", type=int, default=1)
parser.add_argument("-nn", "--NNeighboursSmooth", type=int, default=5)
parserargs = parser.parse_args()


fPedestalMode = "all"
if(parserargs.Mode==1):
    fPedestalMode = "start"
elif(parserargs.Mode==2):
    fPedestalMode = "end"

# compute pedestal in two zones
fPreBins = [0, 800]
fPostBins = [1500, -1]
fSelectedBins=[0, 5000]
fSelectedBins = fPostBins

fChSkip = parserargs.ChSkip

# decide what to plot in which subplot
# in this configuration: bottom plot->baseline mean,  middle plot->temperature, top plot-> baseline RMS
fPlotScheme={"RMS":2, "Temp":1, "B":0}
#fPlotScheme={"Temp":1, "B":0}
# plot only RMS
#fPlotScheme={"RMS":0}
#fPlotScheme={"B":0}
#fPlotScheme={"B":0, "RMS":1}

ChTempDict = {}
EventIDDict = {} 
RunIDDict = {} 
TimeStamp_V = {} 
WvMeanChDict = {}
WvRMSChDict = {}

if(parserargs.Option == 1):
    RunIDDict, EventIDDict, TimeStampDict, WvMeanChDict, WvRMSChDict, ChTempDict = ReadFromAnaROOT(parserargs.Filepath, fPedestalMode, parserargs.ChSkip, parserargs.BoardID, parserargs.Debug, parserargs.NEv, parserargs.Run, parserargs.Ch)
elif(parserargs.Option == 2):
    EventIDDict, WvMeanChDict, WvRMSChDict = ReadFromROOT(parserargs.Filepath, fChSkip, waveformRange=fSelectedBins, maxEvents=parserargs.NEv)
elif(parserargs.Option == 3):
    EventIDDict, WvMeanChDict, WvRMSChDict = ReadFromTxt(parserargs.Filepath)
print("Data...read")

# get run numbers
RunSetV = [item for sublist in list(RunIDDict.values()) for item in sublist]
RunSetV = sorted (set( RunSetV ))
print("Run numbers", RunSetV)

if(parserargs.PlotTime==1):
    PlotAverageBaseline(EventIDDict, TimeStampDict, WvMeanChDict, WvRMSChDict, ChTempDict, fPlotScheme, parserargs.ChSkip, parserargs.UseTimeStamp, parserargs.NNeighboursSmooth)
if(parserargs.MakeStats==1):
    PlotStatisticsBaseline(WvMeanChDict, WvRMSChDict, ChTempDict)

plt.show()