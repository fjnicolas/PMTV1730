# PMTV1730
Soft tools for PMT V1730 analysis



## Running DAQ instructions
* ssh sbnd@sbnd-gateway01.fnal.gov
* ssh sbnd@sbnd-evb02.fnal.gov (open two terminals)
* Setup script (pmt_launch.sh in home area)
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


## Running analyzer instructions

* Purpose: get readable waveforms
* For ana we setup this environment: source ana_launchdaq.sh
* We run this fhicl-file: https://github.com/SBNSoftware/sbndaq-artdaq/blob/develop/sbndaq-artdaq/ArtModules/Common/dump_CAENV1730.fcl
* Correr en /sbnd/fnicolas (“HOME” area)
