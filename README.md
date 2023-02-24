# PMTV1730
Soft tools for PMT V1730 analysis



## Running DAQ instructions
* ssh sbnd@sbnd-gateway01.fnal.gov
* ssh sbnd@sbnd-evb02.fnal.gov (open two terminals)
* Setup script (`pmtv1730_launchdaq.sh` in home area)
* In the first terminal open the monitoring. Command is just DAQInterface 
  * Have a look to the file https://github.com/SBNSoftware/sbndaq/blob/develop/sbn-nd/DAQInterface/boot.txt
  * We have to set it up the eve machine we use
  * We are currently using evb02 and partition 1
  * Message people in the Slack channel #sbnd_daq about what server and partition we are using
* Then we modified configs/standard/EventBuilder2.fcl:
   * These are the hicls to run the evnetbuilder. Tells what evbn server we are. 
   * We modified the path where the output goes (somewhere in /scratch/fnicolas/myfancydirectory)
* Then we have a look to these fhicls: pmtx01.fcl and  pmt_standard.fcl  (apparently we didn’t modify anything?)
* Other file we looked at: https://github.com/SBNSoftware/sbndaq/blob/develop/sbn-nd/DAQInterface/known_boardreaders_list
* Always run in the DAQInterfacee (BOTH MONITROING AND ./RUN)

# Files we needed to modify to run PMT DAQ
* Using `sbndaq v01_05_00`
* sbndaq/sbn-nd/DAQInterface/boot.txt
  * ~~PMT host: sbnd-evb04-daq~~
  * PMT host: sbnd-evb02-daq
  * ~EventBuilder host: sbnd-evb04-daq~
  * ~EventBuilder label: EventBuilder4~
  * EventBuilder host: sbnd-evb02-daq
  * EventBuilder label: EventBuilder2
  * ~Dispatcher host: sbnd-evb04-daq~
  * Dispatcher host: sbnd-evb02-daq
  
  * This file contains the info of the event builder machine (evb) we are use to make our run. We need to make sure anyone else is running on the same machine. There is a [spreadsheet](https://docs.google.com/spreadsheets/d/1xJb7Dge_ktMcXOaUkF0sAGDepLXewSUJwcNiCsHDAOA/edit#gid=0) for that. Modify the event builder machine we are using. Eventually we will want to change the partition in the evb machine we are using (`Subsystem id: 1`)

* sbndaq/sbn-nd/DAQInterface/configs/standard/EventBuilder2.fcl
  * ~outputs.normalOutput.fileName: "/scratch/data/data_evb02_run%R_%#_%to.root"~
  * outputs.normalOutput.fileName: "/scratch/pmtv1730_fnicolas/noiseonly/data_evb02_run%R_%#_%to.root"
  
  * There are 4 `EventBuildeer` fhicls, one per evb server. We need to modify the fhicl corresponding to the server we'll use. Change the output path (where the eveents will be recorded). Make sure to create the directory in the server before running.
  
  

*sbndaq/sbn-nd/DAQInterface/configs/standard/pmtx01.fcl
  * Added daq.fragment_receiver.fragment_id: 15
  
  * This is the fhicl that will call the DAQ generator/producer. It's what we call from the `run` script

*sbndaq/sbn-nd/DAQInterface/run
  * ~-setdaqcomps.sh pmtx03 crt01PULL spectdc ptb01PULL~
  * setdaqcomps.sh pmtx01
  
  * Scrit to start the DAQ. We specify the fhicl with our DAQ configuration. In our case is `pmtx01.fcl`.

## Running analyzer instructions

* Purpose: get readable waveforms
* For ana we setup this environment: source ana_launchdaq.sh
* We run this [fhicl-file](https://github.com/SBNSoftware/sbndaq-artdaq/blob/develop/sbndaq-artdaq/ArtModules/Common/dump_CAENV1730.fcl)
* Correr en /sbnd/fnicolas (“HOME” area)

Notes on long data takings:
Cannot disconnect DAQInterface. Recommend running in tmux or screen so it's running in the background. If we close the monitor window, there's no way to get it back (not critical).
 There's an option in the EventBuilder to prescale the number of recorded events.


## Running CAEN wavedump

- Instructions: [Running wavedump](https://github.com/SBNSoftware/SBNSoftware.github.io/blob/master/running_wvdump.md)
- Useeful commands: `caen` will show info about the connected boards

