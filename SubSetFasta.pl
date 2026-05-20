#!/usr/bin/perl
#
# this subsets and concatenates a set of SNPs from fasta

foreach $i (1..23){
        open(IN,"keepSNPs_max_chrom$i") or die "failed to open snps file $i\n";
        #open(IN,"keepSNPs_chrom$i") or die "failed to open snps file $i\n";
        $j = 0;
        while(<IN>){
                chomp;
                push (@{$snps[$i]},$_);
        }
        close(IN);
}

open(OUT, ">lyc_genomemaxV2.fasta") or die "failed to write\n";

%seq;
foreach $i (1..23){
        open(IN,"sub_max_chrom$i.fasta") or die "failed to open snps file $i\n";
        while(<IN>){
                chomp;
                if(m/^>(\S+)/){
                        $id = $1;
                        if($i == 1){
                                @{$seq{$id}} = ();
                        }
                } else {
                        foreach $snp (@{$snps[$i]}){
                                $c = substr $_,$snp-1, 1;
                                unless(length($c)==1){
                                        print "$c\n";
                                }
                                push(@{$seq{$id}}, $c);
                        }
                }
        }
        close(IN);
}

foreach $pop (sort keys %seq){
        $str = join("",@{$seq{$pop}});
        unless($pop =~ m/rep/){
                $pop =~ s/Lyc-//;
                ##$pop =~ s/\d+//;
                print OUT ">$pop\n";
                print OUT "$str\n";
           }
}
close(OUT);
