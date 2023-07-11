from matplotlib import pyplot as plt
import argparse
import matplotlib.cm as cm
from InputReaderV1730 import *
import matplotlib.patches as patches
import os


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--Filepath", help="Input file name",  default="none")
parser.add_argument("-r", "--Run", help="Select run number", type=int, default=-1)
parser.add_argument("-o", "--Option", help="Input option", type=int, default=1)
parser.add_argument("-n", "--NEv", help="Max Events", type=int, default=1e6)
parser.add_argument("-chSkip", "--ChSkip", type=int, action='append', help="Channels to skip", default=[])
parser.add_argument("-ch", "--Ch", type=int, help="Select channel", default=-1)
parser.add_argument("-b", "--BoardID", type=int, action='append', help="Boards ID to analyze", default=[])
parser.add_argument("-debug", "--Debug", type=int, default=0)

parserargs = parser.parse_args()

# compute pedestal in two zones
fPreBins = [0, 800]
fPostBins = [1500, -1]
fSelectedBins=[0, 5000]
fSelectedBins = fPostBins

fChSkip = parserargs.ChSkip

ChTempDict = {}

_, _, _, _, _, ChTempDict = ReadFromAnaROOT(parserargs.Filepath, "all", parserargs.ChSkip, parserargs.BoardID, parserargs.Debug, parserargs.NEv, parserargs.Run, parserargs.Ch)

fChipDimH = .6
fChipDimW = .9
fNChips = 8
fBoardNumberDict = {160:3,176:5,192:7,208:9,
                    224:11,240:13,256:16,272:18}

fBoardMap = {}
fBoardTempMap = {}

def PlotTemperatureMap(chTempDict):
    steps = 16
    boardIx = 0

    # Looping over channels
    for chIx, ch in enumerate(chTempDict):
        print(ch)
        if(ch in fBoardNumberDict.keys()):
            boardIx = fBoardNumberDict[ch]
        # Checking if the channel is even
        if ch % 2 == 0:
            p = [boardIx, fNChips-(ch % steps)/2 ]
            print(p)
            fBoardMap[ch] = p
            fBoardTempMap[ch] = np.mean( chTempDict[ch] )

    minTemp = min(fBoardTempMap.values())
    maxTemp = max(fBoardTempMap.values())

    params = {'figure.figsize': (10, 8),
        'axes.labelsize': 'x-large',   
        'axes.titlesize': 'x-large',      
        'xtick.labelsize':'xx-large',
        'ytick.labelsize':'xx-large'}
    plt.rcParams.update(params)
    fig = plt.figure("PMT V1730 Temperature Map")
    #fig.subplots_adjust(left=0.075, bottom=0.175, right=0.97, top=0.95, wspace=0.3, hspace=0.45)
    ax = fig.add_subplot(1,1,1)

    # Create a color scale
    color_scale = cm.get_cmap('Reds')

    for ch in fBoardMap:
        color = color_scale( (fBoardTempMap[ch]-minTemp)/(maxTemp-minTemp))
        rec = plt.Rectangle( (fBoardMap[ch][0]-fChipDimW/2, fBoardMap[ch][1]-fChipDimH/2), fChipDimW, fChipDimH, facecolor=color, edgecolor='black')
        ax.add_patch(rec)

    # Add labels and title
    ax.set_xlabel('Slot Index')
    ax.set_ylabel('ADC Chip Index')

    ax.set_xlim(0, 22)
    ax.set_ylim(0, 9)

    ax.set_xticks(np.arange( 1, 23, 2) )
    ax.set_xticklabels(np.arange( 1, 23, 2) )
    ax.set_yticks(np.arange( 1, 9, 1) )
    ax.set_yticklabels(np.arange( 1, 9, 1) )

    cbar = plt.colorbar(cm.ScalarMappable(cmap=color_scale), ax=ax, ticks=np.linspace(0, 1, 4))
    cbar.ax.set_title(r"$\overline{\rm T}$ [$^\circ$C]")
    cbar.ax.set_yticklabels( np.linspace( int(minTemp), int(maxTemp), 4).astype(int).tolist())

    # draw crate limits
    # Draw an empty box using Rectangle patch
    rect = patches.Rectangle((0.5, 0.5), 21, 8, linewidth=2, edgecolor='dimgrey', facecolor='none')

    # Add the box to the plot
    ax.add_patch(rect)

    ax.tick_params(axis='both', which='both', labelsize=20)
    ax.tick_params(axis='both', which='both', length=10, width=2)

    # Remove the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)


    outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/temperaturemap/"
    # create directory if it doesn't exist
    if not os.path.exists(outputFilepath):
        os.makedirs(outputFilepath)
    fig.savefig(outputFilepath+"temperaturemap.pdf")

    plt.show()
    

##########################################

PlotTemperatureMap(ChTempDict)
