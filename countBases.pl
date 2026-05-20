#!/usr/bin/perl

open(IN, "lyc_genomemax.fasta") or die "failed to read\n";
#open(IN, "Lmel_dovetailPacBio_genome.fasta") or die "failed to read\n";
while(<IN>){
	chomp;
	if(m/^>(\S+)/){
		$scaf = $1;
		$cnts{$scaf};
	} else{
		$A = $_ =~ tr/Aa/Aa/;
		$C = $_ =~ tr/Cc/Cc/;
		$G = $_ =~ tr/Gg/Gg/;
		$T = $_ =~ tr/Tt/Tt/;
		$cnts{$scaf}{'a'} += $A;
		$cnts{$scaf}{'c'} += $C;
		$cnts{$scaf}{'g'} += $G;
		$cnts{$scaf}{'t'} += $T;
	}
}

foreach $scaf (sort keys %cnts){
	print "$scaf";
	foreach $base (sort keys %{$cnts{$scaf}}){
		print " $cnts{$scaf}{$base}";
	}
	print "\n";
}
