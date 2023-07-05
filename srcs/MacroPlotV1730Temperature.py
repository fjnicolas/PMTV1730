from matplotlib import pyplot as plt
import argparse
import matplotlib.cm as cm
from PlotUtilsV1730 import *
from InputReaderV1730 import *

params = {'legend.fontsize': 'small',
          'figure.figsize': (12, 8),
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium'}
plt.rcParams.update(params)



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
parser.add_argument("-stats", "--MakeStats", type=int, default=1)
parser.add_argument("-time", "--PlotTime", type=int, default=1)
parser.add_argument("-nn", "--NNeighboursSmooth", type=int, default=5)
parser.add_argument("-grad", "--PlotTGradient", type=int, default=1)
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

fBoardWidth = 2. # cm
fBoardHeight = 30.
fChipDim = fBoardWidth/2


fBoardMap = {}
fBoardTempMap = {}

def PlotTemperatureMap(eventIDDict, chTempDict):
    print(chTempDict)
    nChannels = len(chTempDict.keys())
    steps = 16

    # Initializing variables for X and Y coordinates
    x_coord = 0
    y_coord = 0

    # Looping over channels
    for ch, ch_id in enumerate(chTempDict):
        # Checking if the channel is even
        if ch % 2 == 0:
            # Assigning coordinates for even channels
            x_coord = ((ch // steps) * steps )/16
            y_coord = ch % steps
            # Printing the channel and its assigned coordinates
            print(x_coord, y_coord)
            p = [2*fBoardWidth*x_coord, fBoardHeight*(1-y_coord/16)]
            print(p)
            fBoardMap[ch_id] = p
            fBoardTempMap[ch_id] = np.mean( chTempDict[ch_id] )

    minTemp = min(fBoardTempMap.values())
    maxTemp = max(fBoardTempMap.values())

    print("TTTTT", minTemp, maxTemp)

    fig = plt.figure("PMT V1730 Temperature Map")
    fig.subplots_adjust(left=0.075, bottom=0.175, right=0.97, top=0.95, wspace=0.3, hspace=0.45)
    ax = fig.add_subplot(1,1,1)

    # Create a color scale
    color_scale = cm.get_cmap('Reds')

    for ch in fBoardMap:
        if ch % 2 == 0:

            color = color_scale( (fBoardTempMap[ch]-minTemp)/(maxTemp-minTemp))
            rec = plt.Rectangle( (fBoardMap[ch][0], fBoardMap[ch][1]), fChipDim, fChipDim, facecolor=color, edgecolor='black')
            ax.add_patch(rec)


    # Set the limits of the plot
    ax.set_xlim(-fBoardWidth, 2*8*fBoardWidth)
    ax.set_ylim( -0.1*fBoardHeight, 1.2*fBoardHeight )

    # Add labels and title
    ax.set_xlabel('X [cm]')
    ax.set_ylabel('Y [cm]')

    cbar = plt.colorbar(cm.ScalarMappable(cmap=color_scale), ax=ax, ticks=np.linspace(0, 1, 2))
    cbar.ax.set_title(r"T [$^\circ$C]")
    cbar.ax.set_yticklabels([int(minTemp), int(maxTemp)])  # horizontal colorbar
    plt.show()
    

##########################################

PlotTemperatureMap(EventIDDict, ChTempDict)
