#!/usr/bin/perl
#
# filter bam files to remove short sequences 
#
use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);



foreach $bam (@ARGV){ #takes bam
        $pm->start and next; ## fork
        $o = "filt_$bam";

        my $cmd = qq{ samtools view -h "$bam" | awk 'length(\$10) >= 50 || \$1 ~ /^@/' | samtools view -Sb - > "$o"};
        system ('bash','-lc',$cmd);
        $pm->finish;
}
$pm-> wait_all_children;
