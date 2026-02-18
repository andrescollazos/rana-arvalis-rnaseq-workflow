#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -M pelle
#SBATCH -J STAR_Index
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

module load STAR/2.7.11b-GCC-13.3.0
source ../../.env

mkdir -p $DIR/analyses/01_STAR_mapping/star_index
lfs setstripe -c 8 -S 16M $DIR/analyses/01_STAR_mapping/star_index

STAR \
  --runThreadN $SLURM_CPUS_PER_TASK \
  --runMode genomeGenerate \
  --genomeDir $DIR/analyses/01_STAR_mapping/star_index \
  --genomeFastaFiles $ref_genome_fasta \
  --sjdbGTFfile $ref_genome_gff3 \
  --sjdbGTFtagExonParentTranscript Parent \
  --sjdbOverhang 150