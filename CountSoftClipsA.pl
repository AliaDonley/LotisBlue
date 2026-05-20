#!/usr/bin/perl
#
#softclipping
#
use Parallel::ForkManager;
my $max = 30;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $sam (@ARGV){
        $pm->start and next FILES; ## fork
        $out = $sam;
        $out =~ s/sam/txt/ or die "failed sub here: $out\n";
        system "perl CountSoftClip.pl $sam > $out\n";
        $pm->finish;
}

$pm->wait_all_children;
