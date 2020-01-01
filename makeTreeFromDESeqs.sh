mothur "#filter.seqs(fasta=DESeq16S.silva.align, vertical=T)"

iqtree \
    -s DESeq16S.silva.filter.fasta \
    -alrt 1000 \
    -bb 1000 \
    -nt AUTO \
    -wt

