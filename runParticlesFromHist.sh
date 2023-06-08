#! /bin/bash

varName=${1:-"Part1"}
echo "-----------------------------------"
echo "---working on part "${varName}
echo "-----------------------------------"

root -l -b << EOF
.L makeDumpParticlesFromHistogramLUXE_C.so
//# makeDumpParticlesFromHistogramLUXE("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33_NoECutNtrn.root", 1, 135310, 549773, "${varName}")
makeDumpParticlesFromHistogramLUXE("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33_NoECutNtrn.root", 1, 56555, 210410, "${varName}SmallStat")
EOF