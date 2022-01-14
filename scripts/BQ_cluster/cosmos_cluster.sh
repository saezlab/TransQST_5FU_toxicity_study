#!/bin/sh

#PBS -e /net/data.isilon/ag-saez/bq_agabor/WP8_data_analysis/data/carnival_logs/carnival_${index}.err
#PBS -o /net/data.isilon/ag-saez/bq_agabor/WP8_data_analysis/data/carnival_logs/carnival_${index}.out
#PBS -l nodes=2,mem=30gb

# start conda
source /net/data.isilon/ag-saez/bq_agabor/SOFTWARE/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate r_env

# set wd on project
cd /net/data.isilon/ag-saez/bq_agabor/WP8_data_analysis

# create log directory
mkdir -p /net/data.isilon/ag-saez/bq_agabor/WP8_data_analysis/data/carnival_logs/

# launch script, printing output on default output file
/net/data.isilon/ag-saez/bq_mgarrido/SOFTWARE/miniconda3/envs/r_env/bin/Rscript scripts/cosmos_cluster.R ${index}

