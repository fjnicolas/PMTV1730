#include "plotone.h"
#include "plotone.C"

void MacroPlotV1730(int nentries=10, std::string file_name = "processed_data_result.root", std::string dir_name = "caenv1730dump"){

  plotone *fPlotOne = new plotone(file_name, dir_name);

  for(unsigned int i=0; i<nentries; i++){
    fPlotOne->Loop(i, 1);
  }

  return 0; 
}