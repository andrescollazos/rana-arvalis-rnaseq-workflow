#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J STAR_Mapping
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH --array=0-68
#SBATCH -t 3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

module load STAR/2.7.11b-GCC-13.3.0
source ../../.env

manifest="$meta_data/01_mapping/reads_manifest.tsv"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" "$manifest")
read1=$(echo "$line" | cut -f1)
read2=$(echo "$line" | cut -f2)
sample_id=$(echo "$line" | cut -f3)
sample_folder="$DIR/analyses/01_STAR_mapping/out/${sample_id}"
mkdir -p "$sample_folder"

STAR \
  --runThreadN $SLURM_CPUS_PER_TASK \
  --genomeDir $DIR/analyses/01_STAR_mapping/star_index \
  --readFilesIn "$read1" "$read2" \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outFileNamePrefix "$sample_folder/" \
  --outSAMtype BAM SortedByCoordinate