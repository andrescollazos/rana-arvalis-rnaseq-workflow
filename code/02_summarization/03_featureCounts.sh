#!/bin/bash -l
#SBATCH -A uppmax2025-2-482
#SBATCH -p pelle
#SBATCH -J featureCounts
#SBATCH -e %x.%A_%a.er
#SBATCH -o %x.%A_%a.out
#SBATCH -t 30:00
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

source ../../.env

manifest="$meta_data/01_mapping/reads_manifest.tsv"

bam_files=$(awk -v dir="$DIR/analyses/01_STAR_mapping/out" 'NR>0 {print dir"/"$3"/Aligned.sortedByCoord.out.bam"}' "$manifest")

featureCounts \
	-p --countReadPairs \
    -B -C \
	-F GTF \
	-t exon \
	-g gene_id \
	-a "$ref_genome_gtf" \
	-G "$ref_genome_fasta" \
	-T $SLURM_CPUS_PER_TASK \
	-s 1 \
	-o "$DIR/analyses/02_summarization/featureCounts.txt" \
	$bam_files