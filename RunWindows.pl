#!/usr/bin/perl
# Zach's RunWindows.pl with quartet tag, binary and prefix as arguments
($q, $bin, $prefix) = @ARGV;
die "usage: RunWindows.pl QUARTET_TAG WINDOW_BINARY ALN_PREFIX\n" unless defined $prefix;
die "window binary not found: $bin\n" unless -x $bin;
for $i (1..23){
$aln = "Sub${q}_${prefix}$i.fasta";
die "missing $aln\n" unless -e $aln;
system("$bin $aln Mapping$q.txt > winout$i.tsv") == 0
or die "window run failed on chromosome $i\n";
}
