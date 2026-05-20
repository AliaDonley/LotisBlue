#!/usr/bin/perl

## this is to figure out which SNPs are variable in the fasta

foreach $i (1..23){
    system "grep -v \"^>\" sub_max_chromV2ad1_beastfiltered_fff_o_lycpool_chrom$i.fasta | perl -pe 'tr/ACGTN/12345/' | sed 's/./& /g' > text_chrom$i.fasta\n";
}
