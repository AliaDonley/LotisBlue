#!/usr/bin/perl
#
# filter vcf with GATK and tabix 
#

use Parallel::ForkManager;
my $max = 26;
my $pm = Parallel::ForkManager->new($max);



my $ref = "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";
my $gatk = "/uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar"; 


foreach $vcf (@ARGV){
        $pm->start and next; ## fork
        $in = $vcf; ##don't need to gunzip here..?
        $in =~ s/\.gz//; ## drops .gz from the filename
        $o = "fff_$vcf"; ## outfile with name fff_xxxx.vcf.gz ##if not zipped, switch variable to $in      

        system "~/bin/tabix $vcf\n"; ##tabix needs bgzipped (should be .gz)
        system "bgzip -d $vcf\n"; ## unzipe file- correct
        
        system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar IndexFeatureFile -I $in\n";
        system "java -jar /uufs/chpc.utah.edu/sys/installdir/gatk/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar VariantFiltration -R $ref -V $in -O $o --filter-name \"bqbz\" --filter-expression \"BQBZ > 3.0 || BQBZ < -3.0\" --filter-name \"mqbz\" --filter-expression \"MQBZ > 3.0 || MQBZ < -3.0\" --filter-name \"rpbz\" --filter-expression \"RPBZ > 3.0 || RPBZ < -3.0\" --filter-name \"depth\" --filter-expression \"DP < 1350\" --filter-name \"mapping\" --filter-expression \"MQ < 30\" --verbosity ERROR\n";
        system "bgzip -d $o\n";

        $pm->finish;

}

$pm->wait_all_children;
