#!/usr/bin/perl
#
# keep only a subset of taxa
#

open(IN, "SubTaxa.txt") or die;
while(<IN>){
        chomp;
        $keep{$_} = 1;
}
close(IN);

foreach $fa (@ARGV){
        open(IN, $fa) or die "failed to read $fa\n";
        open(OUT, "> sub_$fa") or die "failed to write for $fa\n";
        while(<IN>){
                chomp;
                if(m/^>(\S+)/){
                        $id = $1;
                        if(defined $keep{$id}){
                                print OUT "$_\n";
                                $a = <IN>;
                                print OUT $a;
                        }
                }
        }
        close(IN);
        close(OUT);
}
