#!/bin/bash
#SBATCH -J trueval-sce1
#SBATCH -p icelake
#SBATCH --time=36:0:0
#SBATCH -A mrc-bsu-sl2-cpu
#SBATCH --cpus-per-task=67
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=jml219@cam.ac.uk

. /etc/profile.d/modules.sh   
module purge                    
module load rhel8/default-icl    

module load R/4.3.1-icelake

Rscript true_HR_MRD_gen_singletrial.R 

