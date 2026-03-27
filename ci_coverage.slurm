#!/bin/bash
#SBATCH -J CI_simu_weak_conf
#SBATCH -p icelake
#SBATCH --time=10:0:0
#SBATCH -A mrc-bsu2-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk
#SBATCH --array=300,500,1000,2500

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.3.1-icelake

Rscript mle_ci_comp_time_weak_conf.R
