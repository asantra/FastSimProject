#! /bin/bash

# root -l << EOF
# .L makeLUXEFastSimFullSimDumpPlotsFromText.C++
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_NoECutNtrn_Processed_Sorted.txt", 1.0, 32,false)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_NoECutNtrn_Processed_Sorted.txt", 1.0, 31, false)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId32.txt", 1.0, 32, true)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId31.txt", 1.0, 31, true)
# EOF

for weight in 1.0 1.1 1.5 2.0 2.5; do
    #weight=${1:-"1.0"}
echo "working on weight: "${weight}
root -l << EOF
    .L makeDumpPlotsFromText.C++
    makeDumpPlotsFromText("${weight}")
EOF

root -l << EOF
    .L makeDumpParticlesFromHistogramLUXE.C++
    makeDumpParticlesFromHistogramLUXE("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33_NoECutNtrn_CoarseBinning_${weight}timesBackwardInThetaMore3AndRLess300.root", 549773, 135310)
EOF
python plotFastSamplingFullSim.py -times ${weight}
done
