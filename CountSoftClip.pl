#!/usr/bin/perl
#
# counts the total number of soft-clipped bases in a sam file
# also reports the number of reads with at least one soft-clipped base
# and the total number of reads
#
# usage: CountSoftClip.pl myfile.sam
# or redirect the output to a file like this
# usage: CountSoftClip.pl myfile.sam > myoutput.txt

## set up counters
$tcnt = 0; ## total number of clipped bases
$scnt = 0; ## total number of sequences with clipped bases
$xcnt = 0; ## total number of sequences

## get the file and open it
$sam = shift(@ARGV); 
open(IN, $sam) or die "failed to read the sam file: $sam\n";
while(<IN>){
        chomp;
        unless(m/^\@/){ # skip header lines
                @line = split(/\s+/,$_); # split line
                ## CIGAR string is in column 6
                @match = ($line[5] =~ m/(\d+)S/g);
                $cnt = 0; # reset counter
                foreach $s (@match){## count soft-clipped bases
                        $cnt += $s;
                }
                $tcnt += $cnt;
                if($cnt > 0){
                        $scnt++; ## increment sequence counter
                }
                $xcnt++; ## increment total sequence counter
        }
}
print "Done processing $sam\n";
print "Found $tcnt soft-clipped bases\n";
print "$scnt out of $xcnt sequences contained at least one soft-clipped base\n";
close(IN);
