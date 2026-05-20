#!/usr/bin/perl
#

for $i (1..23){
        system "../wins/MASTERWORK/bin/slidingwindow SubABMxSINxTBY_sub_max_chrom$i.fasta MappingABMxSINxTBY.txt > winout$i.tsv\n";
}
