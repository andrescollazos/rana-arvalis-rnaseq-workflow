#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -p pelle
#SBATCH -J gff_to_gtf
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 30:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL

module load gffread/0.12.7-GCCcore-13.3.0

source ../../.env

gffread "$ref_genome_gff3" -T -F -o "$DIR/data/metadata/moor_frog_annotation.gtf" 