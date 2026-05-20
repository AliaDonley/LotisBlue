#!/usr/bin/perl
#
# keep only a subset of taxa
#

## subfile contains taxa to keep
$subfile = shift(@ARGV);
open(IN, $subfile) or die;
while(<IN>){
        chomp;
        $keep{$_} = 1;
}
close(IN);

$prefix = $subfile;
$prefix =~ s/\.txt// or die "faield at $prefix\n";

foreach $fa (@ARGV){
        open(IN, $fa) or die "failed to read $fa\n";
        open(OUT, "> $prefix"."_$fa") or die "failed to write for $fa\n";
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
