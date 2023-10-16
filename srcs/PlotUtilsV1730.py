import numpy as np
from matplotlib import pyplot as plt
import ROOT
import pandas as pd
from datetime import datetime
import matplotlib.dates as mdates
import matplotlib.cm as cm
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
import os


plt.rcParams['figure.figsize'] = [13, 8]
from scipy.optimize import curve_fit


fColor1 = "midnightblue"
fColor2 = "maroon"

fBaseline = 0
fNChannels = 16
fNOpDet = 320
fWfSize = 5000
fPlotStackedHistograms=True

# Define Gaussian function
def gaussian(x, mu, sigma, A):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


def GetValuesFromGaussFit(data, binScheme):
    # Create histogram
    hist, bin_edges = np.histogram(data, bins=binScheme, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Fit Gaussian function to histogram
    p0 = [np.mean(data), np.std(data), np.max(hist)]
    popt, pcov = curve_fit(gaussian, bin_centers, hist, p0=p0)

    return popt[0], popt[1]


def calculate_local_variation(vector, n_neighbors):
    local_variations = []
    for i in range(len(vector)):
        start_index = max(0, i - n_neighbors)
        end_index = min(i + n_neighbors + 1, len(vector))
        local_values = vector[start_index:end_index]
        local_variation = np.max(local_values) - np.min(local_values)
        sign = np.sign(local_values[-1] - local_values[0])
        local_variations.append(sign * local_variation)
    return local_variations


def nearest_neighbors_smoothing(data, window_size):
    smoothed_data = []
    for i in range(len(data)):
        start_index = max(0, i - window_size + 1)
        end_index = min(i + window_size, len(data))
        neighbors = data[start_index:end_index]
        smoothed_value = sum(neighbors) / len(neighbors)
        smoothed_data.append(smoothed_value)
    return smoothed_data


def GetBiasAndStdDev(x, y, nbins, minentries_perbin):
    bin1x=min(x) 
    bin2x=max(x)
    hp=ROOT.TProfile("", ";", nbins, bin1x, bin2x, "s")
    for i in range(min(len(x), len(y))):
        hp.Fill( x[i], y[i] )

    X=[]; Mean=[]; StdDev=[];
    for k in range(hp.GetNbinsX()):
        if(hp.GetBinEntries(k)>=minentries_perbin):
            X.append(hp.GetBinCenter(k));
            Mean.append(hp.GetBinContent(k));
            StdDev.append(hp.GetBinError(k));
    return X, Mean, StdDev

# Make legend using a cmap
def DrawChannelColorMap(cax, colorDict):
    colorMap = LinearSegmentedColormap.from_list("ChannelColorMap", list(colorDict.values()), N=len(colorDict.values()))
    cbar = plt.colorbar(cm.ScalarMappable(cmap=colorMap), cax=cax, ticks=np.linspace(0, 1, 4), 
                        orientation="horizontal")
    cbar.ax.set_xticklabels(np.linspace( min(colorDict.keys()), max(colorDict.keys()), 4).astype(int).tolist())
    cbar.ax.set_ylabel("Ch", fontsize=10)


def AddSlotNumbers(ax, minCh, maxCh):

    fBoardNumberDict = {160:"S3",176:"S5",192:"S7",208:"S9",
                    224:"S11",240:"S13",256:"S16",272:"S18"}
    ax.set_xticks(np.arange( minCh, maxCh, 16 ))
    for index, xmaj in enumerate(ax.xaxis.get_majorticklocs()):
        if(index<len(ax.xaxis.get_majorticklocs())):
            ax.axvline(x=xmaj, ls="--", linewidth = 2., color = 'dimgrey')
            ax.text(xmaj, 1.05, fBoardNumberDict[xmaj], ha='center', va='center', transform=ax.get_xaxis_transform(), color= 'dimgrey', rotation="45")


def SaveSubplots(fig, axs, outputFilepath, filename, expFrac=1.05):
    # save full figure
    fig.savefig(outputFilepath+"pdf/"+filename+".pdf")
    fig.savefig(outputFilepath+"png/"+filename+".png")
    
    # Iterate over the subplots
    nrows = len(axs[0])
    ncols = len(axs)

    # save individual plots
    count = 0
    for i in range(nrows):
        for j in range(ncols):
            
            # format for cropping is [[x0, y0], [x1, y1]]
            dims = [ [i*1./nrows, j*1./ncols], [(i+1.)/nrows, (j+1)*1./ncols] ]
            #png
            fig.savefig(outputFilepath+"png/"+filename+f'_sp{count}'+".png", bbox_inches=mtransforms.Bbox(dims).
            transformed(fig.transFigure - fig.dpi_scale_trans))
            #pdf
            fig.savefig(outputFilepath+"pdf/"+filename+f'_sp{count}'+".pdf", bbox_inches=mtransforms.Bbox(dims).
            transformed(fig.transFigure - fig.dpi_scale_trans))
            count+=1
            # other method
            #fig.savefig(outputFilepath+"pdf/"+filename+f'_sp{i+1}'+".pdf", bbox_inches=extent.expanded(1.06, 1.22))


##### COLOR LIST #####
#derived from https://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
def GetColorList(N):
    NUM_COLORS = N
    LINE_STYLES = ['solid']#, 'dashed', 'dashdot', 'dotted']
    NUM_STYLES = len(LINE_STYLES)

    ColorList = []
    #cm = plt.get_cmap('gist_rainbow')
    #cm = plt.get_cmap('viridis')
    cm = plt.get_cmap('nipy_spectral')
    for i in range(NUM_COLORS):
        ColorList.append(cm(i//NUM_STYLES*float(NUM_STYLES)/NUM_COLORS))
        #lines[0].set_linestyle(LINE_STYLES[i%NUM_STYLES])

    return ColorList


#### Get stats string
def GetStatsStr(v, l0="", entries=True, mu=True, stddev=True, prec="{:.1f}"):
    if(l0!=""): l0+="\n"
    if(entries): l0+=r"Entries="+str(len(v))+"\n"
    if(mu): l0+=r"$\mu$="+prec.format(np.mean(v))
    if(stddev): l0+=r", $\sigma$="+prec.format(np.std(v))
    return l0

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
def PlotAverageBaseline(eventIDDict, timeStampDict, wvMeanChDict, wvRMSChDict, chTempDict, plotScheme={"RMS":2, "Temp":1, "B":0}, chSkip=[], useTimeStamp=True, nNeighboursSmooth=1):


    # Get channel-color map
    cList=GetColorList(len(wvMeanChDict.keys()))
    cDict = {}
    for ch_ix, ch in enumerate(wvMeanChDict.keys()):
        cDict[ch] = cList[ch_ix]
    fMarkerSize = 1.
    
    fNN = nNeighboursSmooth

    params = {'figure.figsize': (13, 8),
         'axes.labelsize': 'x-large',   
         'axes.titlesize': 'x-large',      
         'xtick.labelsize':'large',
         'ytick.labelsize':'x-large'}
    
    plt.rcParams.update(params)
    fig = plt.figure("PMT V1730")
    fig.subplots_adjust(left=0.1, bottom=0.175, right=0.97, top=0.95, wspace=0.3, hspace=0.45)
    gs = gridspec.GridSpec(3, 8, hspace=0.1, wspace=1.5)

    # define the subplots
    ax0 = fig.add_subplot(gs[2,0:7])
    axs = [ax0]

    if(len(plotScheme)>=2):
        ax1 = fig.add_subplot(gs[1,0:7], sharex=ax0)
        plt.setp(ax1.get_xticklabels(), visible=False)
        axs.append(ax1)
    if(len(plotScheme)==3):
        ax2 = fig.add_subplot(gs[0,0:7], sharex=ax1)
        plt.setp(ax2.get_xticklabels(), visible=False)
        axs.append(ax2)

    # define axis for legend
    axLeg = fig.add_subplot(gs[2,7])
    axLeg.axis("off")

    # channel counter
    channelIndex = 0

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        print("Using channel: ", ch)
    
        # xaxis
        xValues = eventIDDict[ch]
        if(useTimeStamp):
            xValues =  []
            for ts in timeStampDict[ch]:
                xValues.append( datetime.fromtimestamp(ts) ) 
        
        if("B" in plotScheme):
            pedForPlot = nearest_neighbors_smoothing(wvMeanChDict[ch], fNN)
            axs[ plotScheme["B"]].scatter(xValues, pedForPlot, color=cDict[ch], marker='o', s=fMarkerSize)
        if("RMS" in plotScheme):
            rmsForPlot = nearest_neighbors_smoothing(wvRMSChDict[ch], fNN)
            axs[ plotScheme["RMS"]].scatter(xValues, rmsForPlot, color=cDict[ch], marker='o', s=fMarkerSize)
        if("Temp" in plotScheme):
            tempForPlot = nearest_neighbors_smoothing(chTempDict[ch], fNN)
            axs[ plotScheme["Temp"]].scatter(xValues, tempForPlot, color=cDict[ch], marker='o', s=fMarkerSize)

        channelIndex+=1


    if("B" in plotScheme):
        axs[plotScheme["B"]].set_ylabel("B [ADC]"); 
        axs[plotScheme["B"]].grid()
        # force matplotliob to not use scientific notation for the baseline
        axs[plotScheme["B"]].yaxis.set_major_formatter(FormatStrFormatter('% 1i'))
    if("RMS" in plotScheme):
        axs[plotScheme["RMS"]].set_ylabel(r"$\sigma_{\rm B}$ [ADC]");     
        axs[plotScheme["RMS"]].grid()
    if("Temp" in plotScheme):
        axs[plotScheme["Temp"]].set_ylabel(r"T [$^\circ$C]"); 
        axs[plotScheme["Temp"]].grid() 
        

    if(useTimeStamp):
        #fmt = mdates.DateFormatter('%D:%H:%M:%S')
        fmt = mdates.DateFormatter('%b %d - %H:%M:%S')
        axs[0].xaxis.set_major_formatter(fmt)
        axs[0].tick_params(axis='x', labelrotation = 45) 
    else:
        axs[0].set_xlabel("event ID");
   
    # draw the color map legend
    colorMap = LinearSegmentedColormap.from_list("ChannelColorMap", list(cDict.values()), N=len(cDict.values()))
    cbar = plt.colorbar(cm.ScalarMappable(cmap=colorMap), ax=axLeg, ticks=np.linspace(0, 1, 4), location="left")
    cbar.ax.set_yticklabels(np.linspace( min(cDict.keys()), max(cDict.keys()), 4).astype(int).tolist())   
    cbar.ax.set_ylabel("Ch", fontsize=10)

   
    outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/baseline/"
    # create directory if it doesn't exist
    if not os.path.exists(outputFilepath):
        os.makedirs(outputFilepath)
    fig.savefig(outputFilepath+"average_pedestal.png", dpi=300)
##########################################



##########################################
# plot baseline mean and temperature
def PlotBaselineTemperatureCorrelation(eventIDDict, timeStampDict, wvMeanChDict, chTempDict, nPlots=2, chSkip=[], useTimeStamp=True, nNeighboursSmooth=1):

    
    # plotting style parameters
    fMarkerSize = .3
    params = {'figure.figsize': (13, 8),
         'axes.labelsize': 'x-large',   
         'axes.titlesize': 'x-large',      
         'xtick.labelsize':'large',
         'ytick.labelsize':'x-large'}
   
    # number of channels
    nChannels = len(wvMeanChDict.keys())

    fNN = nNeighboursSmooth
    
    # first figure: correlation baseline-temperature
    plt.rcParams.update(params)
    fig, AxGrad = plt.subplots(1, 1, figsize=(7, 7))
    fig.subplots_adjust(left=0.15, bottom=0.1, right=0.94, top=0.92)

    # second figure: baseline+temperature time series
    fig2 = plt.figure("PMT V1730: Baseline and temperature")
    fig2.subplots_adjust(left=0.1, bottom=0.175, right=0.91, top=0.91, wspace=0.6, hspace=0.75)
    
    # plot the time series only for a subsbample of channles
    chGradStep = int(nChannels/nPlots)
    gs2 = gridspec.GridSpec(nPlots, 1)
    axGrad = []
    for ix in range(nPlots):
        axGrad.append( fig2.add_subplot(gs2[ix,0]) )
    axGradCont=0

    # list to store the baseline and temperature derivatives
    tempGradV = []
    pedGradV = []

    # plot!
    channelIndex = 0
    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        print("Using channel: ", ch)

        # get smoothed time series
        tempForPlot = nearest_neighbors_smoothing(chTempDict[ch], fNN)
        pedForPlot = nearest_neighbors_smoothing(wvMeanChDict[ch], fNN)
        xValues = eventIDDict[ch]
        
        # xaxis
        if(useTimeStamp):
            xValues =  []
            for ts in timeStampDict[ch]:
                xValues.append( datetime.fromtimestamp(ts) ) 

        # get derivatives
        pedGradient = calculate_local_variation(pedForPlot, 2)
        tempGradient = calculate_local_variation(tempForPlot, 2)
        
        # store values for the all-channels plot
        pedGradV.extend(pedGradient)
        tempGradV.extend(tempGradient)
        
        # plot time series for selected channel
        if(ch%chGradStep==0):
            axTwin = axGrad[axGradCont].twinx()
            plot1 = axGrad[axGradCont].scatter(xValues, pedForPlot, color="navy", label="Pedestal", marker='o', s=fMarkerSize)
            plot2 = axTwin.scatter(xValues, tempForPlot, color="maroon", label="Temperature", marker='^', s=fMarkerSize)
            axGrad[axGradCont].set_ylabel("B [ADC]", color="navy")
            axGrad[axGradCont].tick_params('y', colors="navy")
            axGrad[axGradCont].ticklabel_format(useOffset=False, style='plain', axis='y')
            axTwin.set_ylabel(r"T [$^\circ$ C]", color="maroon")
            axTwin.tick_params('y', colors="maroon")
            
            if(useTimeStamp):
                fmt = mdates.DateFormatter('%D:%H:%M:%S')
                axGrad[axGradCont].xaxis.set_major_formatter(fmt)
                axGrad[axGradCont].tick_params(axis='x', labelrotation = 45)
            else:
                axGrad[axGradCont].ax2.set_xlabel("event ID");
    
            labels = [plot1.get_label(), plot2.get_label()]
            plt.legend([plot1, plot2], labels, title="Ch="+str(ch))

            axGradCont+=1

        channelIndex+=1

    # plot profile for baseline-temperature correlation
    gradX, gradY, gradYErr = GetBiasAndStdDev(tempGradV, pedGradV, 10, 10)
    #AxGrad.scatter(tempGradV, pedGradV, s=fMarkerSize, c="grey", alpha=0.5)
    AxGrad.errorbar(gradX, gradY, gradYErr, ls='none', marker="o", color="navy")
    AxGrad.set_xlabel(r"$\Delta_{\rm T} \ [^\circ C]$", fontsize=20)
    AxGrad.set_ylabel(r"$\Delta_{\rm B} \ [{\rm ADC}]$", fontsize=20)


    # save plots
    outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/baselinetempcorrelation/"
    # create directory if it doesn't exist
    if not os.path.exists(outputFilepath):
        os.makedirs(outputFilepath)
    fig.savefig(outputFilepath+"average_pedestal_gradientAll.pdf")
    fig2.savefig(outputFilepath+"average_pedestal_gradient.png", dpi=300)
##########################################

##########################################
# plot baseline statistics
def PlotStatisticsBaseline(wvMeanChDict, wvRMSChDict, chTempDict, chSkip=[], makeErrorsWithFit=False):

    fUseGaussFit = makeErrorsWithFit
    # data frame with channel information
    df = pd.DataFrame(columns=['Ch', 'ChLoc',
                               'ChPed','ChPedErr', 'ChPedMin', 'ChPedMax',
                               'ChRMS','ChRMSErr', 'ChRMSMin', 'ChRMSMax',
                               'ChTemp', 'ChTempErr', 'ChTempMin', 'ChTempMax'
                               ])
    # dataframe in string mode
    dfString = pd.DataFrame(columns=["Ch", "Ped", "RMS", "T"])

    cList=GetColorList(len(wvMeanChDict.keys()))
    cDict = {}
    for ch_ix, ch in enumerate(wvMeanChDict.keys()):
        cDict[ch] = cList[ch_ix]

    fBinSizeRMS = 0.001
    fBinSizePed = 0.25
    fBinSizeTemp = 1
    

    allChannelsPedestalMean = []
    allChannelsPedestalRMS = []
    allChannelsTemp = []
    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        ch_loc = ch%2
        allChannelsPedestalMean.extend(wvMeanChDict[ch])
        allChannelsPedestalRMS.extend(wvRMSChDict[ch])
        if(ch_loc==0):
            allChannelsTemp.extend(chTempDict[ch])
        # pedestal mean variables
        ped_mean, ped_err = 0, 0
        ped_min = np.min(wvMeanChDict[ch])
        ped_max = np.max(wvMeanChDict[ch])
        if(fUseGaussFit==1):
            ped_mean, ped_err = GetValuesFromGaussFit(wvMeanChDict[ch], np.arange(ped_min, ped_max, fBinSizePed))
        else:
            ped_mean = np.mean(wvMeanChDict[ch])
            ped_err = np.std(wvMeanChDict[ch])
        
        
        
        
        # pedestal RMS variables
        rms_min = np.min(wvRMSChDict[ch])
        rms_max = np.max(wvRMSChDict[ch])
        rms_mean, rms_err = 0, 0
        if(fUseGaussFit==1):
            rms_mean, rms_err = GetValuesFromGaussFit(wvRMSChDict[ch], np.arange(rms_min, rms_max, fBinSizeRMS))
        else:
            rms_mean = np.mean(wvRMSChDict[ch])
            rms_err = np.std(wvRMSChDict[ch])

        # temperature variables
        t_mean = np.mean(chTempDict[ch])
        t_err = np.std(chTempDict[ch])
        t_min = np.min(chTempDict[ch])
        t_max = np.max(chTempDict[ch])
        
        df.loc[int(ch)] = [int(ch), int(ch_loc), 
                      ped_mean, ped_err, ped_min, ped_max,
                      rms_mean, rms_err, rms_min, rms_max,
                      t_mean, t_err, t_min, t_max
                      ]

        dfString.loc[ch] = [
            str(ch),
            "{:.1f}".format(ped_mean)+r"$\pm$"+"{:.1f}".format(ped_err),
            "{:.2f}".format(rms_mean)+r"$\pm$"+"{:.2f}".format(rms_err),
            "{:.2f}".format(t_mean)+r"$\pm$"+"{:.2f}".format(t_err)
        ]
    
    
    fMarkerSize = 1.

    fBinSizeRMS = 0.05
    fBinSizePed = 2.5
    fBinSizeTemp = 1
    binsRMS=np.arange( df["ChRMSMin"].min(), df["ChRMSMax"].max(), fBinSizeRMS )
    binsPed=np.arange( df["ChPedMin"].min(), df["ChPedMax"].max(), fBinSizePed  )
    binsTemp=np.arange( df["ChTempMin"].min(), df["ChTempMin"].max(), fBinSizeTemp  )

    # Color map used for legend
    minCh = min(wvMeanChDict.keys())
    maxCh = max(wvMeanChDict.keys())
    myCmap = LinearSegmentedColormap.from_list("ChannelColorMap", list(cDict.values()), N=len(cDict.values()))

    params = {'figure.figsize': (13, 8),
        'axes.labelsize': 'x-large',   
        'axes.titlesize': 'large',      
        'xtick.labelsize':'large',
        'ytick.labelsize':'large'}
    
    plt.rcParams.update(params)


    #fig1, axs1 = plt.subplots(2, 3, constrained_layout=False)
    fig1 = plt.figure(constrained_layout=False)
    gs1 = gridspec.GridSpec(2, 3) 
    fig1.subplots_adjust(left=0.07, bottom=0.075, right=0.99, top=0.95, wspace=0.5, hspace=0.45)
    axs1 = [[plt.subplot(gs1[0]), plt.subplot(gs1[1]),plt.subplot(gs1[2])], 
            [plt.subplot(gs1[3]), plt.subplot(gs1[4]),plt.subplot(gs1[5])]]

    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        labelStr=""
        if(ch%4==0):
            labelStr="Ch"+str(ch)
        axs1[0][0].hist(wvMeanChDict[ch], bins=binsPed, label=labelStr, histtype="step", edgecolor=cDict[ch], stacked=True)
        axs1[1][0].hist(wvRMSChDict[ch], bins=binsRMS, label=labelStr, histtype="step", edgecolor=cDict[ch], stacked=True)
    
    axs1[0][0].set_xlabel("B [ADC]"); 
    axs1[0][0].set_ylabel("# entries"); 
    axs1[0][0].tick_params(axis='x', labelrotation = 15)
    axs1[0][0].grid()    
    divider = make_axes_locatable(axs1[0][0])
    cax = divider.new_vertical(size = '2.5%', pad = 0.25)
    fig1.add_axes(cax)
    DrawChannelColorMap(cax, cDict)
    """cbar = plt.colorbar(cm.ScalarMappable(cmap=myCmap), cax=cax, orientation = 'horizontal', ticks=np.linspace(0, 1, 4))
    cbar.ax.set_xticklabels(np.linspace(minCh, maxCh-1, 4).astype(int).tolist(), fontsize=10)
    cbar.ax.set_ylabel("Ch", fontsize=10)"""

    axs1[1][0].set_xlabel(r"$\sigma_{\rm B}$ [ADC]"); axs1[1][0].set_ylabel("# entries"); 
    axs1[1][0].grid()
    axs1[1][0].set_yscale("log")
    divider = make_axes_locatable(axs1[1][0])
    cax = divider.new_vertical(size = '2.5%', pad = 0.25)
    fig1.add_axes(cax)
    DrawChannelColorMap(cax, cDict)

    # baseline RMS as a function of the channel
    axs1[1][1].errorbar(df.Ch, df.ChRMS, df.ChRMSErr, ls='none', marker="o", c=fColor1)
    axs1[1][1].grid()
    axs1[1][1].set_xlabel("Channel");
    axs1[1][1].set_ylabel(r"$\overline{\sigma_{\rm B}}$ [ADC]"); 
    #axs1[1][1].set_xticks(np.arange( min(df.Ch), max(df.Ch), fNLabelsForCh))

    # baseline mean as a function of the channel
    gs1A = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[2,1.2], subplot_spec=gs1[0, 1], hspace=0)
    axs1A_1 = fig1.add_subplot(gs1A[0, 0]); 
    axs1A_2 = fig1.add_subplot(gs1A[1, 0], sharex = axs1A_1);
    plt.setp(axs1[0][1], visible=False) 
    axs1A_1.errorbar(df.Ch, df.ChPed, df.ChPedErr, ls='none', marker="o", c=fColor1)
    axs1A_2.plot(df.Ch, df.ChPedErr, ls='none', marker="o",  c=fColor2)
    axs1A_1.grid(); axs1A_2.grid()
    axs1A_2.set_xlabel("Channel");
    axs1A_1.set_ylabel(r"$\overline{\rmB}$ [ADC]");
    axs1A_2.set_ylabel(r"$\Delta\overline{\rm B}$ [ADC]");
    stats=GetStatsStr(allChannelsPedestalMean, entries=False)
    #axs1[0][1].set_xticks(np.arange( min(df.Ch), max(df.Ch), fNLabelsForCh))

    # all channels histogram
    axs1[0][2].hist(allChannelsPedestalMean, label="All ch", histtype="step", bins=binsPed, color=fColor1)
    axs1[0][2].grid()
    axs1[0][2].set_xlabel(r"B [ADC]");
    axs1[0][2].set_ylabel("# entries");
    axs1[0][2].tick_params(axis='x', labelrotation = 15)
    stats=GetStatsStr(allChannelsPedestalMean, entries=False)
    axs1[0][2].legend(title=stats)

    axs1[1][2].hist(allChannelsPedestalRMS, label="All ch", histtype="step", bins=binsRMS,  color=fColor1)
    axs1[1][2].grid()
    axs1[1][2].set_yscale("log")
    axs1[1][2].set_xlabel(r"$\sigma_{\rm B}$ [ADC]"); axs1[1][2].set_ylabel("# entries");
    stats=GetStatsStr(allChannelsPedestalRMS, entries=False, prec="{:.2f}" )
    axs1[1][2].legend(title=stats)




    fig2, axs2 = plt.subplots(2, 3)
    fig2.subplots_adjust(left=0.07, bottom=0.075, right=0.99, top=0.95, wspace=0.5, hspace=0.45)
    

    # RMS min and max
    axs2[0][0].plot(df.Ch, df.ChRMSMin, ls='none', marker="o", label="Min",  c=fColor1)
    axs2[0][0].plot(df.Ch, df.ChRMSMax, ls='none', marker="o", label="Max",  c=fColor2)
    axs2[0][0].grid()
    axs2[0][0].set_xlabel(r"Channel");
    axs2[0][0].set_ylabel(r"$\sigma_{\rm B}$ [ADC]");
    #ax1.set_xticks(np.arange( min(df.Ch), max(df.Ch), fNLabelsForCh))
    axs2[0][0].legend()

    # RMS Spread
    axs2[0][1].plot(df.ChRMSMax-df.ChRMSMin, ls='none', marker="o", label="Min",  c=fColor1)
    axs2[0][1].grid()
    axs2[0][1].set_xlabel(r"Channel");
    axs2[0][1].set_ylabel(r"${\rm D}_{\sigma_{\rm B}}$ [ADC]");


    # Baseline min and max
    axs2[1][0].plot(df.Ch, df.ChPedMin, ls='none', marker="o", label="Min",  c=fColor1)
    axs2[1][0].plot(df.Ch, df.ChPedMax, ls='none', marker="o", label="Max",  c=fColor2)
    axs2[1][0].grid()
    axs2[1][0].set_xlabel(r"Channel");
    axs2[1][0].set_ylabel("B [ADC]");
    #axs2[1][0].set_xticks(np.arange( min(df.Ch), max(df.Ch), fNLabelsForCh))
    axs2[1][0].legend()

    # Baseline spread
    axs2[1][1].plot(df.ChPedMax-df.ChPedMin, ls='none', marker="o", label="Min",  c=fColor1)
    axs2[1][1].grid()
    axs2[1][1].set_xlabel(r"Channel");
    axs2[1][1].set_ylabel(r"${\rm D}_{\rm B}$ [ADC]");



    for ch in wvMeanChDict.keys():
        if(ch in chSkip): continue
        labelStr=""
        if(ch%4==0):
            labelStr="Ch"+str(ch)
        axs2[0, 2].scatter(wvMeanChDict[ch], wvRMSChDict[ch], label=labelStr, color=cDict[ch], s=fMarkerSize)
    axs2[0, 2].set_xlabel("B [ADC]")
    axs2[0][2].tick_params(axis='x', labelrotation = 15)
    axs2[0, 2].set_ylabel(r"$\sigma_{\rm B}$ [ADC]")
    divider = make_axes_locatable(axs2[0][2])
    cax = divider.new_vertical(size = '2.5%', pad = 0.25)
    fig2.add_axes(cax)
    DrawChannelColorMap(cax, cDict)

    # table    
    """axs2[1, 2].table(cellText=dfString.values, colLabels=dfString.columns, loc='center')
    axs2[1, 2].axis('off')"""


    if(fPlotStackedHistograms):
        histMeanV = []
        histRMSV = []
        colorV = []
        labelV = [] 
        for ch in wvMeanChDict.keys():
            if(ch in chSkip): continue
            labelStr=""
            if(ch%4==0):
                labelStr="Ch"+str(ch)
            labelV.append(labelStr)
            colorV.append(cDict[ch])
            histMeanV.append(wvMeanChDict[ch])
            histRMSV.append(wvRMSChDict[ch])
        
        """axs2[0][2].hist(histMeanV, label=labelV, color=colorV, stacked=True)
        axs2[0][2].set_xlabel("B [ADC]"); axs[0][0].set_ylabel("# entries"); 
        axs2[0][2].grid()
        axs2[0][2].legend()"""
        binsRMS=np.arange( df["ChRMSMin"].min(), 4, fBinSizeRMS )
        axs2[1][2].hist(histRMSV, bins=binsRMS, label=labelV, color=colorV, stacked=True)
        axs2[1][2].set_xlabel(r"$\sigma_{\rm B}$ [ADC]");
        axs2[1][2].set_ylabel("# entries"); 
        axs2[1][2].grid()
        divider = make_axes_locatable(axs2[1][2])
        cax = divider.new_vertical(size = '2.5%', pad = 0.25)
        fig2.add_axes(cax)
        DrawChannelColorMap(cax, cDict)




    # temperature plots
    fig3, axs3 = plt.subplots(2, 3)
    fig3.subplots_adjust(left=0.07, bottom=0.075, right=0.99, top=0.95, wspace=0.5, hspace=0.45)
    axs3[0][0].errorbar(df.Ch, df.ChTemp, df.ChTempErr, ls='none', marker="o", label="Average",  c=fColor1)
    axs3[0][0].grid()
    axs3[0][0].set_xlabel(r"Channel");
    axs3[0][0].set_ylabel(r"$\overline{\rm T}$ [$^\circ$ C]");
    #axs3[0][0].set_xticks(np.arange( min(df.Ch), max(df.Ch), fNLabelsForCh))
    AddSlotNumbers(axs3[0,0], minCh, maxCh)
    axs3[0][0].legend()

    axs3[1][0].errorbar(df.Ch, df.ChTempMin, 1, ls='none', marker="o", label="Min",  c=fColor1)
    axs3[1][0].errorbar(df.Ch, df.ChTempMax, 1, ls='none', marker="o", label="Max",  c=fColor2)
    axs3[1][0].grid()
    axs3[1][0].set_xlabel(r"Channel");
    axs3[1][0].set_ylabel(r"T [$^\circ$ C]");
    #axs3[1][0].set_xticks(np.arange( min(df.Ch), max(df.Ch), fNLabelsForCh))
    AddSlotNumbers(axs3[1,0], minCh, maxCh)
    axs3[1][0].legend()

    # temperature as a function of the channel
    axs3[0][1].hist(allChannelsTemp, label="All ch", histtype="step", bins=binsTemp,  color=fColor1)
    axs3[0][1].grid()
    axs3[0][1].set_xlabel(r"T [$^\circ$ C]");
    axs3[0][1].set_ylabel("# entries");
    stats=GetStatsStr(allChannelsTemp, entries=False)
    axs3[0][1].legend(title=stats)


    # temperature dispersion
    axs3[1][1].plot(df.ChTempMax-df.ChTempMin, ls='none', marker="o", label="Min",  c=fColor1)
    axs3[1][1].grid()
    axs3[1][1].set_xlabel(r"Channel")
    axs3[1][1].set_ylabel(r"${\rm D}_{\rm T}$ [$^\circ$C]")
    axs3[1][1].set_xticks(np.arange( minCh, maxCh, 16 ))
    AddSlotNumbers(axs3[1,1], minCh, maxCh)
    axs3[1][1].legend()


    """axs3[0][2].errorbar(df.ChTemp[df.ChLoc==0], df.ChPed[df.ChLoc==0], df.ChPedErr[df.ChLoc==0], ls='none', marker="o", label="Even",  c=fColor1)
    axs3[0][2].errorbar(df.ChTemp[df.ChLoc==1], df.ChPed[df.ChLoc==1], df.ChPedErr[df.ChLoc==1], ls='none', marker="o", label="Odd",  c=fColor1)
    axs3[0][2].legend()"""
    axs3[0][2].scatter(df.ChTemp, df.ChPed, marker="o", c=fColor1)
    axs3[0][2].grid()
    axs3[0][2].set_xlabel(r"$\overline{\rm T}$  [$^\circ$ C]")
    axs3[0][2].set_ylabel(r"$\overline{\rm B}$ [ADC]")


    """axs3[1][2].errorbar(df.ChTemp[df.ChLoc==0], df.ChRMS[df.ChLoc==0], df.ChRMSErr[df.ChLoc==0], ls='none', marker="o", label="Even",  c=fColor1)
    axs3[1][2].errorbar(df.ChTemp[df.ChLoc==1], df.ChRMS[df.ChLoc==1], df.ChRMSErr[df.ChLoc==1], ls='none', marker="o", label="Odd",  c=fColor2)
    axs3[1][2].legend()"""
    axs3[1][2].scatter(df.ChTemp, df.ChRMS, marker="o",  c=fColor1)
    axs3[1][2].grid()
    axs3[1][2].set_xlabel(r"$\overline{\rm T}$  [$^\circ$ C]")
    axs3[1][2].set_ylabel(r"$\overline{\sigma_{\rm B}}$ [ADC]")


    outputFilepath = os.path.dirname(os.path.realpath(__file__))+"/../plots/stats/"
    # Create the directory if it doesnt exist
    if not os.path.exists(outputFilepath):
        os.makedirs(outputFilepath)
    os.makedirs(outputFilepath+"png/", exist_ok=True)
    os.makedirs(outputFilepath+"pdf/", exist_ok=True)
    
    # save the data frame
    df.to_csv(outputFilepath+"baseline_stats.csv")

    SaveSubplots(fig1, axs1, outputFilepath, "statistics_baseline1" )
    SaveSubplots(fig2, axs2, outputFilepath, "statistics_baseline2" )
    SaveSubplots(fig3, axs3, outputFilepath, "statistics_temperature" )
##########################################

