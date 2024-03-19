#!/bin/bash
#SBATCH -J cont-single
#SBATCH -p icelake
#SBATCH --time=5:0:0
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk
#SBATCH --array=1-27

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.3.1-icelake

Rscript simu_continuous_ipw_cali_seq_single.R
