#! /bin/bash

root -l << EOF
.L makeLUXEFastSimFullSimDumpPlotsFromText.C++
makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_NoECutNtrn_Processed_Sorted.txt", 1.0, 32,false)
makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_NoECutNtrn_Processed_Sorted.txt", 1.0, 31, false)
makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId32.txt", 1.0, 32, true)
makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId31.txt", 1.0, 31, true)
EOF