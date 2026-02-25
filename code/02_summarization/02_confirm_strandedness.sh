#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -p pelle
#SBATCH -J determine_strand
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 30:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL

module load RSeQC/5.0.4 BEDOPS/2.4.41

source ../../.env

gtf2bed --attribute-key=transcript_id < "$ref_genome_gtf" > "$DIR/data/metadata/moor_frog_annotation.bed"

infer_experiment.py \
    -r "$DIR/data/metadata/moor_frog_annotation.bed" \
    -i "$DIR/analyses/01_STAR_mapping/out/P32262_272/Aligned.sortedByCoord.out.bam"