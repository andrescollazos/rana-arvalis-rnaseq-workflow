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

# 1. Build Ensembl X. tropicalis protein BLAST database
if [ ! -f Xenopus_tropicalis.UCB_Xtro_10.0.pep.all.fa.pin ]; then
    makeblastdb \
      -in Xenopus_tropicalis.UCB_Xtro_10.0.pep.all.fa \
      -dbtype prot
fi

# 2. Run BLAST against Ensembl proteome
blastp \
  -query proteins.fa \
  -db Xenopus_tropicalis.UCB_Xtro_10.0.pep.all.fa \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads "${SLURM_CPUS_PER_TASK:-1}" \
  -out "blast_results_${SLURM_JOB_ID}.tsv"

# 3. Keep best hit per query
sort -k1,1 -k6,6gr -k5,5g "blast_results_${SLURM_JOB_ID}.tsv" | \
awk '!seen[$1]++' > "blast_best_${SLURM_JOB_ID}.tsv"