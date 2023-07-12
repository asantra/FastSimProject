#! /bin/bash

weightNt=1.005
weightPh=4.0
version="12"
echo "working on weight: "${weight}
echo "Neutron weight: "${weightNt}
echo "Photon weight: "${weightPh}
# root -l -b << EOF
# .L makeLUXEFastSimFullSimDumpPlotsFromText.C++
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_${weightNt}timesNeutron_${weightPh}timesPhoton_v${version}.txt", 1.0, 32,false)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FastSim_0p06BX_${weightNt}timesNeutron_${weightPh}timesPhoton_v${version}.txt", 1.0, 31, false)
# EOF
python plotFastSimFullSim.py -timesNt ${weightNt} -timesPh ${weightPh} -version ${version} -det 32
python plotFastSimFullSim.py -timesNt ${weightNt} -timesPh ${weightPh} -version ${version} -det 31


# root -l -b << EOF
# .L makeLUXEFastSimFullSimDumpPlotsFromText.C++
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId32.txt", 1.0, 32, true)
# makeLUXEFastSimFullSimDumpPlotsFromText("/Volumes/OS/LUXEBkgOutputFile/Geant4Files/OutputFile/LUXEDumpFiles_FullSim_0p06BX_DetId31.txt", 1.0, 31, true)
# EOF