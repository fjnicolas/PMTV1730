import numpy as np
import ROOT

fBaseline = 0
fNChannels = 16
fNOpDet = 320
fWfSize = 5000

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
def ReadFromAnaROOT(filename, mode="all", ch_skip=[], fBoardIDList=[], fVerbose=0, fMaxEvents=-1, fRunID=-1, fCh=-1):

    print("Reading from AnaTree")
    file = ROOT.TFile.Open(filename)
    file.Print()
    tree_path="caenv1730dump/nt_wvfm"#"caenv1730ana/nt_wvfm;1"
    tuple =  file.Get(tree_path)
   
    wvMeanChDict = {}
    wvRMSChDict = {}
    chTempDict = {}
    eventIDDict = {}
    runIDDict = {}
    timeStampDict = {}

    # initialize dictionaries
    #for ch in range(fNChannels):
    for ch in range(fNOpDet):
        eventIDDict[ch]=[]
        runIDDict[ch]=[]
        timeStampDict[ch]=[]
        wvMeanChDict[ch]=[]
        wvRMSChDict[ch]=[]
        chTempDict[ch]=[]

    eventCounter = 0
    previousEvent=-1
    for tree_entry in range( tuple.GetEntries() ):
        if(fMaxEvents!=-1 and eventCounter>fMaxEvents): continue        
        tuple.GetEntry(tree_entry)
        
        if(int(tuple.ch) in ch_skip): continue
        if( (tuple.boardId in fBoardIDList)==False and len(fBoardIDList)>0): continue

        if(fRunID!=-1 and tuple.art_run!=fRunID): continue
        if(fCh!=-1 and tuple.ch!=fCh): continue
        eventIDDict[int(tuple.ch)].append(tuple.art_ev)
        runIDDict[int(tuple.ch)].append(tuple.art_run)
        timeStampDict[int(tuple.ch)].append(tuple.stamp_time)

        if(mode=="all"):
            wvMeanChDict[int(tuple.ch)].append(tuple.ped)
            wvRMSChDict[int(tuple.ch)].append(tuple.rms)
        elif(mode=="start"):
            wvMeanChDict[int(tuple.ch)].append(tuple.pedStart)
            wvRMSChDict[int(tuple.ch)].append(tuple.rmsStart)
        elif(mode=="end"):
            wvMeanChDict[int(tuple.ch)].append(tuple.pedEnd)
            wvRMSChDict[int(tuple.ch)].append(tuple.rmsEnd)
        chTempDict[int(tuple.ch)].append(tuple.temp)
        if(tuple.art_ev!=previousEvent):
            eventCounter+=1
        previousEvent=tuple.art_ev

        if(fVerbose>=1): 
            print(int(tuple.ch), tuple.ped, tuple.rms, tuple.temp)
            print(int(tuple.art_run), int(tuple.art_ev), int(tuple.caen_ev), "Ch:",int(tuple.ch), tuple.stamp_time)

    eventIDDictFinal = {}
    runIDDictFinal = {}
    timeStampDictFinal = {}
    wvMeanChDictFinal = {}
    wvRMSChDictFinal = {}
    chTempDictFinal = {}
    for ch in range(fNOpDet):
        if(len(eventIDDict[ch])>0):
            eventIDDictFinal[ch] =  np.array(eventIDDict[ch])
            runIDDictFinal[ch] =  np.array(runIDDict[ch])
            timeStampDictFinal[ch] =  np.array(timeStampDict[ch])
            wvMeanChDictFinal[ch] = np.array(wvMeanChDict[ch])
            wvRMSChDictFinal[ch] =  np.array(wvRMSChDict[ch])
            chTempDictFinal[ch] =  np.array(chTempDict[ch])
        
    return runIDDictFinal, eventIDDictFinal, timeStampDictFinal, wvMeanChDictFinal, wvRMSChDictFinal, chTempDictFinal
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