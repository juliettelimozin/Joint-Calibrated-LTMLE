#!/bin/bash
#SBATCH -J singletrial-med-sce2
#SBATCH -p icelake
#SBATCH --time=36:0:0
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk
#SBATCH --array=1-9

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.3.1-icelake

Rscript HPC_code_singletrial_scenario2.R 
