## Join pdfs:
`pdfjam synthesis.pdf paper1.pdf paper2.pdf paper3.pdf paper4.pdf --outfile papers.pdf`

## Convert pdf to booklet
`pdfbook2 output.pdf`

## Add blank pages to a pdf
`gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER \
  -sOutputFile=output.pdf \
  toprint.pdf \
  <(printf "<</PageSize [595 842]>> setpagedevice %.0s showpage" $(seq 1 $NUMBER_PAGES))`
