#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J Samtools_Sorting
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH --array=0-38
#SBATCH -t 3:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL

module load samtools/1.20
source ../../.env

manifest="$meta_data/01_mapping/reads_manifest_sort.tsv"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$manifest")
sample_id=$(echo "$line" | cut -f3)
sample_folder="$DIR/analyses/01_STAR_mapping/out/${sample_id}"

set -e

samtools sort -@ $SLURM_CPUS_PER_TASK -m 4G \
  "$sample_folder/Aligned.out.bam" \
  -o "$sample_folder/Aligned.sortedByCoord.out.bam"

samtools index "$sample_folder/Aligned.sortedByCoord.out.bam"

rm "$sample_folder/Aligned.out.bam"
