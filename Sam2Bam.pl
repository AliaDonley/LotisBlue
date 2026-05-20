#!/usr/bin/perl
#
# alignment with bwa mem 
#
use Parallel::ForkManager;
my $max = 30;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $bam (@ARGV){
        $pm->start and next FILES; ## fork
        $out = $bam;
        $out =~ s/bam/sam/ or die "failed sub here: $out\n";
        system "samtools view -h -o $out $bam\n";
        $pm->finish;
}

$pm->wait_all_children;
