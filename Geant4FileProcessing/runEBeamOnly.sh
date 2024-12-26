#! /bin/bash

# For Tracks tree:
echo "FastSim/FullSim samples"

### put the text file containing the root file from Geant4 processing
filename="LUXEFastSim_PreparedFromPart2C.txt"

#### this block is to read one line of a list at a time
counter=0
increment=1

while IFS= read -r line
do
  counter=$(( $counter + $increment ))
  echo "$line  $counter"
  
  root -l -b -q process_track_tree_draw_v8.C\(\"$line\",\"$counter\",\"$filename\"\)
done < "$filename"


duration=$SECONDS
echo "Total time taken for this process ---- $(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
