#!/bin/bash
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J MultiQC
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 30:00
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL

module load MultiQC/1.28-foss-2024a

source ../../.env

cd "$DIR/analyses/01_STAR_mapping"

multiqc out/*/Log.final.out