# Using GAN simulation to prepare FastSim samples
1. First, make a Geant4 like tree from GAN output (CSV file):
- Run makeDumpParticlesFromGANLUXE.C
- Go inside root and run the following commands:
```
root -l
.L makeDumpParticlesFromGANLUXE.C++
makeDumpParticlesFromGANLUXE(<give the location and name of the CSV file generated by GAN>)
```

2. This will give a root tree which will be compatible with Geant4
