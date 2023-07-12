#! /bin/bash



# for weight in 1.0 1.1 1.5 2.0 2.5; do
# for weight in 3.0; do
    #weight=${1:-"1.0"}
weightNt=1.005
weightPh=4.0
version="12"
echo "working on weight: "${weight}
echo "Neutron weight: "${weightNt}
echo "Photon weight: "${weightPh}
root -l << EOF
    .L makeDumpPlotsFromText.C++
    makeDumpPlotsFromText("${weightNt}", "${weightPh}","${version}")
EOF

root -l << EOF
    .L makeDumpParticlesFromHistogramLUXE.C++
    makeDumpParticlesFromHistogramLUXE("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId33_NoECutNtrn_CoarseBinning_${weightNt}timesNeutron_${weightPh}timesPhoton_v${version}.root", 1099546, 270620)
EOF
python plotFastSamplingFullSim.py -times 1.0 -timesNt ${weightNt} -timesPh ${weightPh} -version ${version}
# done
