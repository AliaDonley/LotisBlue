#!/usr/bin/perl
#

foreach $i (1..23){

        system "/uufs/chpc.utah.edu/common/home/u6047808/bin/ASTER-Linux/bin/caster-site -i sub_CASTchromad1_fff_o_lycpool_chrom$i.fasta -o CASTcout_max_$i --root Lyc-MEN12 --thread 24\n";
#       system "ASTER-Linux/bin/caster-site -i sub_CASTchromad1_fff_o_lycpool_chrom$i.fasta -o CASTcout_max_$i --root MEN --thread 24\n";
}




#For afterGATK filtering
#!/usr/bin/perl
#
foreach $i (1..23){
        system "/uufs/chpc.utah.edu/common/home/u6047808/bin/ASTER-Linux/bin/caster-site -i CAST_chrom_noBAT49_filteredV2_$i.fasta -o cout_noBAT49_filteredV2_$i --root MEN12 --thread 24\n";
}
