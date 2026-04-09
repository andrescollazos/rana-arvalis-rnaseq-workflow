#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -p pelle
#SBATCH -J BLAST
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

module load BLAST+/2.17.0-gompi-2024a

source ../../.env

cd $DIR/code/05_GO_annotation

# 1. Build database if needed
if [ ! -f xtrop_proteins.fasta.pin ]; then
    makeblastdb -in xtrop_proteins.fasta -dbtype prot
fi

# 2. Run BLAST
blastp \
  -query proteins.fa \
  -db xtrop_proteins.fasta \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads $SLURM_CPUS_PER_TASK \
  -out blast_results_${SLURM_JOB_ID}.tsv

sort -k1,1 -k6,6gr -k5,5g blast_results_${SLURM_JOB_ID}.tsv | \
awk '!seen[$1]++' > blast_best_${SLURM_JOB_ID}.tsv