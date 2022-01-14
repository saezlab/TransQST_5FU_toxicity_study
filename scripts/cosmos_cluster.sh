#!/bin/sh

#PBS -e /beegfs/work/hd_az242/WP8_data_analysis/WP8_data_analysis/data/carnival_logs/carnival_${index}_${rep}.err
#PBS -o /beegfs/work/hd_az242/WP8_data_analysis/WP8_data_analysis/data/carnival_logs/carnival_${index}_${rep}.out
#PBS -l nodes=1:ppn=1,mem=120gb,walltime=48:00:00
#PBS -N run_cosmos_n1_128gb${index}${rep}

# start conda
source /beegfs/work/hd_az242/SOFTWARE/miniconda3/etc/profile.d/conda.sh

# activate conda env
conda activate r_env

# set wd on project
cd /beegfs/work/hd_az242/WP8_data_analysis/WP8_data_analysis

# create log directory
mkdir -p /beegfs/work/hd_az242/WP8_data_analysis/WP8_data_analysis/data/carnival_logs/

# launch script, printing output on default output file
/beegfs/work/hd_az242/SOFTWARE/miniconda3/envs/r_env/bin/Rscript scripts/cosmos_cluster.R ${index} ${rep}

