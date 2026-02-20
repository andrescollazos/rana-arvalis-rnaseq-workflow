echo -e "Sample\tInputReads\tUnique%\tMulti%\tTooManyLoci%\tUnmappedTooShort%\tUnmappedOther%\tMismatchRate%\tAvgMappedLength" > summary.tsv

for f in *.Log.final.out; do
  sample=${f%.Log.final.out}

  input=$(grep "Number of input reads" "$f" | awk '{print $NF}')
  unique=$(grep "Uniquely mapped reads %" "$f" | awk '{print $NF}')
  multi=$(grep "% of reads mapped to multiple loci" "$f" | awk '{print $NF}')
  toomany=$(grep "% of reads mapped to too many loci" "$f" | awk '{print $NF}')
  unmapped_short=$(grep "% of reads unmapped: too short" "$f" | awk '{print $NF}')
  unmapped_other=$(grep "% of reads unmapped: other" "$f" | awk '{print $NF}')
  mismatch=$(grep "Mismatch rate per base" "$f" | awk '{print $NF}')
  avg_len=$(grep "Average mapped length" "$f" | awk '{print $NF}')

  echo -e "$sample\t$input\t$unique\t$multi\t$toomany\t$unmapped_short\t$unmapped_other\t$mismatch\t$avg_len" >> summary.tsv
done