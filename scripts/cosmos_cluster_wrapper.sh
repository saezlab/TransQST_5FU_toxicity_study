#!/bin/sh


# number of elements in the RDS file containing the COSMOS inputs
nElements=7
nReps=10

# for each index
for j in $(seq 1 $nReps); do

  for i in $(seq 1 $nElements); do

    # Avoiding submission of more jobs than permitted
    #while [ $(qstat -u $USER | wc -l) -ge 35 ]; do
    #    sleep 60
    #done

    # print message
    echo "Submitting job $i / $nElements, rep: $j / $nReps"

    # submit job using index
    qsub -q short scripts/cosmos_cluster.sh -v index=${i},rep=${j}

    sleep 5
  done
done
