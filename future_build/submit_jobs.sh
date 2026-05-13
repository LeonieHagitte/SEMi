#!/bin/bash
#SBATCH --job-name=mnlfa_trees_sim
#SBATCH --array=1-100
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

mkdir logs
Rscript simulation_rt.R ${SLURM_ARRAY_TASK_ID} 100