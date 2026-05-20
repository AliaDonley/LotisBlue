#!/usr/bin/perl
#
# run mapdamage
#


use Parallel::ForkManager;
my $max = 30;
my $pm = Parallel::ForkManager->new($max);

foreach $file (@ARGV){
        $pm->start and next; ## fork
        $file =~ m/_([a-zA-Z0-9]+)\.bam/;
        $id = $1;
        $out = "output_$id"."_filteredMapDamage";
        print "mapDamage -i $file -r /uufs/chpc.utah.edu/common/home/gompert-group3/data/butterflygenomes/Lmelissa/Lmel_dovetailPacBio_genome.fasta.masked -d $out\n";
        system "mapDamage -i $file -r /uufs/chpc.utah.edu/common/home/gompert-group3/data/butterflygenomes/Lmelissa/Lmel_dovetailPacBio_genome.fasta.masked -d $out\n";
        $pm->finish;
}
$pm->wait_all_children;
