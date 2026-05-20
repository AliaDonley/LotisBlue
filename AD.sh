#!/usr/bin/bash
#
# extract allele depth AD from biallelic SNPs that passed filtering
#

for f in fff*vcf
do
        echo "Processing $f"
        out="$(echo $f | sed -e 's/vcf/txt/')"
        echo "Output is ad1_$out"
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\1/g' > ad1_$out
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\2/g' > ad2_$out
done
