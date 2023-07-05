from matplotlib import pyplot as plt
import argparse
from PlotUtilsV1730 import *
from InputReaderV1730 import *


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--Filepath", help="Input file name",  default="none")
parser.add_argument("-r", "--Run", help="Select run number", type=int, default=-1)
parser.add_argument("-n", "--NEv", help="Max Events", type=int, default=1e6)
parser.add_argument("-timeStamp", "--UseTimeStamp", help="Use time stamp", type=int, default=1)
parser.add_argument("-chSkip", "--ChSkip", type=int, action='append', help="Channels to skip", default=[])
parser.add_argument("-ch", "--Ch", type=int, help="Select channel", default=-1)
parser.add_argument("-b", "--BoardID", type=int, action='append', help="Boards ID to analyze", default=[])
parser.add_argument("-m", "--Mode", type=int, help="Read pedestal from all/start/end of the waveform", default=0)
parser.add_argument("-debug", "--Debug", type=int, default=0)
parser.add_argument("-nn", "--NNeighboursSmooth", type=int, default=5)
parser.add_argument("-nplots", "--NPlots", type=int, default=2)
parserargs = parser.parse_args()


fPedestalMode = "all"
if(parserargs.Mode==1):
    fPedestalMode = "start"
elif(parserargs.Mode==2):
    fPedestalMode = "end"

fChSkip = parserargs.ChSkip

ChTempDict = {}
EventIDDict = {} 
RunIDDict = {} 
TimeStamp_V = {} 
WvMeanChDict = {}
WvRMSChDict = {}


RunIDDict, EventIDDict, TimeStampDict, WvMeanChDict, WvRMSChDict, ChTempDict = ReadFromAnaROOT(parserargs.Filepath, fPedestalMode, parserargs.ChSkip, parserargs.BoardID, parserargs.Debug, parserargs.NEv, parserargs.Run, parserargs.Ch)

print("Data...read")

# get run numbers
RunSetV = [item for sublist in list(RunIDDict.values()) for item in sublist]
RunSetV = sorted (set( RunSetV ))
print("Run numbers", RunSetV)


PlotBaselineTemperatureCorrelation(EventIDDict, TimeStampDict, WvMeanChDict, ChTempDict, parserargs.NPlots,parserargs.ChSkip, parserargs.UseTimeStamp, parserargs.NNeighboursSmooth)

plt.show()