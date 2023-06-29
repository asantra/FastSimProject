#! /bin/bash

# root -l << EOF
# .L makeLUXEFastSimFullSimDumpPlotsFromText.C++
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_NoECutNtrn_Processed_Sorted.txt", 1.0, 32,false)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_NoECutNtrn_Processed_Sorted.txt", 1.0, 31, false)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId32.txt", 1.0, 32, true)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId31.txt", 1.0, 31, true)
# EOF

# for weight in 1.0 1.1 1.5 2.0 2.5; do
for weight in 3.0; do
    #weight=${1:-"1.0"}
weightNt=1.3
weightPh=2.5
echo "working on weight: "${weight}
echo "Neutron weight: "${weightNt}
echo "Photon weight: "${weightPh}
# root -l << EOF
#     .L makeDumpPlotsFromText.C++
#     makeDumpPlotsFromText("${weightNt}", "${weightPh}")
# EOF

root -l << EOF
    .L makeDumpParticlesFromHistogramLUXE.C++
    makeDumpParticlesFromHistogramLUXE("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33_NoECutNtrn_CoarseBinning_${weightNt}timesNeutron_${weightPh}timesPhoton_BackwardInThetaMore3AndRLess300.root", 549773, 135310)

EOF
python plotFastSamplingFullSim.py -times 1.0 -timesNt ${weightNt} -timesPh ${weightPh}
done
