#!/bin/sh


# number of elements in the RDS file containing the CARNIVAL inputs
nElements= ccc

# for each index
for i in $(seq 1 $nElements); do

  # Avoiding submission of more jobs than permitted
  while [ $(qstat -u $USER | wc -l) -ge 35 ]; do
        sleep 60
  done
  
  # print message
  echo "Submitting job $i / $nElements"
  
  # submit job using index
  qsub scripts/carnival_cluster.sh -v index=${i}

done