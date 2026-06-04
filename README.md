Divergence time estimation for the extinct Lotis blue butterfly. Using Zach Gompert's LycAdmixMosaic method for genomic analysis and population structure analysis. 

# Data

Raw pool-seq data currently in: /uufs/chpc.utah.edu/common/home/gompert-group5/data/lycaeides_poolseq/Alignment. The set of samples I used here are form a single round of sequencing

population    Samples(N)    Nominal taxon
Insert table here

| Population | Sample (n) | Species |
|------------|------------|---------|
| ABM        | ABM20 (48) | L. melissa |
| BAT		 | BAT20 (	  |         |
| BCR        | BCR17 (48) | JH (admixed) |
| BHP        | BHP19(48)  | L. melissa |
| BKM        | BKM19 (33) | Warners (admixed) |
| BTB	  	 | BTB17 (48) | JH (admixed) |
| CLH        | CLH19 (36) | White Mt. (admixed) |
| CP         | CP19  (48) | Sierra (admixed) |
| EP         | EP19  (48) | Warners (admixed) |
| GNP		 | GNP17 (56) | L. iads  |
| HJ         | HJ20 (48)  | L. melissa |
| HNV        | HNV17 (48) | JH (admixed) |
| LOTIS      | LOTIS (14)   | L. anna lotis |
| LS		 | LS19 (48	  | L. anna |
| MEN        | MEN12 (10) | L. argyrognomon (France) |
| MR         | MR20 (48)  | Sierra (admixed) |
| MTU        | MTU20 (48) | L. melissa |
| SBW		 | SBW18 (20) | L. iads (Alaska) |
| SHC        | SHC11 (46) | L. anna ricei |
| SIN        | SIN10 (48) | L. melissa |
| SUV        | SUV20 (51) | L. melissa|
| TBY		 | TBY51 ()	  | L. idas sublivens |
| TBY        | TBY11 (24) | L. idas sublivens|
| TIC        | TIC19 (48) | Sierra (admixed) |
| VE         | VE20 (48)  | L. melissa |
| YG		 | YG20 (48)  | L. anna |


The outgroup, MEN, is Plebejus argyrognomon. The data were generated and cleaned up by BGI (with soapnuke). Here is the report from BGI: BGI_F22FTSUSAT0310-01_LYCgpswR_report_en.pdf. The sequence data (in bam format) are now in the NCBI SRA (PRJNA1375794). 

#check report numbers here


# DNA Sequence Alignment

We aligned the DNA sequence data to the updated PacBio _L. melissa_ genome. I used bwa-mem2 version #### (https://github.com/bwa-mem2/bwa-mem2)
I had access to the already indexed reference genome done by Zach Gompert which can be found in: /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome. 
## index genome with bwa-mem2
```sh
/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta
```
We set up the alignment using the submission script
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bwa-mem2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignments

perl BwaMemFork.pl ../F22FTSUSAT0310-01_LYCgpswR/soapnuke/clean/*/*1.fq.gz 

```


Which runs
```perl
#!/usr/bin/perl
#
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";

FILES:
foreach $fq1 (@ARGV){
	$pm->start and next FILES; ## fork
	$fq2 = $fq1;
	$fq2 =~ s/_1\.fq\.gz/_2.fq.gz/ or die "failed substitution for $fq1\n";
        $fq1 =~ m/clean\/([A-Za-z0-9]+)/ or die "failed to match id $fq1\n";
	$ind = $1;
	$fq1 =~ m/([A-Za-z_\-0-9]+)_1\.fq\.gz$/ or die "failed match for file $fq1\n";
	$file = $1;
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:Lyc-"."$ind\\tLB:Lyc-"."$ind\\tSM:Lyc-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$file.bam - && samtools index -@ 2 $ind"."_$file.bam\n";

	$pm->finish;
}
```
## Zach aligned ^^^
All downstream analysis can be found in /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

We used samtools (version ##) to sort and index the alignments. Following this, the .bam files for each population were merged using samtools (version ##). This was submitted using the shell script:
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=merge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Alignment

perl ../Scripts/MergeFork.pl
``` 
Which runs:
```pl
#!/usr/bin/perl
#
# merge alignments for each population sample with samtools version XX 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

open(IDS,"pids.txt");
while(<IDS>){
	chomp;
	push(@IDs,$_);
}
close(IDS);

FILES:
foreach $id (@IDs){
	$pm->start and next FILES; ## fork
        system "samtools merge -c -p -o Merged/$id.bam $id"."_*.bam\n";
	system "samtools index -@ 2 Merged/$id.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```
# Removing PCR Duplicates
We used samtools (version ##) to remove PCR duplicates following the standard protocol.  I am using the default option (same as -m t) to measure positions based on template start/end. And I am using -r to not just mark but remove duplicates. The submission script is:
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=dedup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)


cd /scratch/general/nfs1/dedup

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts/RemoveDupsFork.pl *bam
```
Which runs
```pl
#!/usr/bin/perl
#
# PCR duplicate removal with samtools
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $bam (@ARGV){
	$pm->start and next FILES; ## fork
	$bam =~ m/^([A-Za-z0-9]+)/ or die "failed to match $bam\n";
	$base = $1;
	system "samtools collate -o co_$base.bam $bam /scratch/general/nfs1/dedup/t$bam\n";
	system "samtools fixmate -m co_$base.bam fix_$base.bam\n";
	system "samtools sort -o sort_$base.bam fix_$base.bam\n";
	## using default definition of dups
	## measure positions based on template start/end (default). = -m t
	system "markdup -T /scratch/general/nfs1/dedup -r sort_$base.bam dedup_$base.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```
# Filter ends 
Filter (ct and ga) before mapdamage. Tried doing it in perl, gave up and just ran interactively with 
FilterFork.sh
```sh
#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert #SBATCH --partition=kingspeak
#SBATCH --job-name=filter
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load samtools
cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis
for bam in *.bam; do base=$(basename "$bam" .bam) out="${base}.filtered${minlen}.bam"
samtools view -h "$bam"
| awk -v m="$minlen" '($1 ~ /^@/) || (length($10) >= 50)'
| samtools view -b -o "$out" -

samtools index "$out" done
```
Output: 







# MapDamage
We ran all aligned and indexed samples through mapdamage2 to assess if the aDNA was contaminated. The data needed to be filtered before being run through mapdamage to eliminate noise. We did one run with all bases and one run where we filtered out the softclipped bases. To do this we used SubMapDamFork.sh


Which runs MapDamFork.pl



# Variant Calling

Variants were called with bcftools (version##). We did not perform INDEL realignment--WHY?
The following submission script was submitted:
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bcf_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /scratch/general/nfs1/dedup

perl /uufs/chpc.utah.edu/common/home/gompert-group2/data/Lycaeides_poolSeq/Scripts/BcfForkLg.pl chrom*list 
```
FIX TO MY PATH

which runs
```pl
#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#
use Parallel::ForkManager;
my $max = 26;
my $pm = Parallel::ForkManager->new($max);

my $genome ="/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";

foreach $chrom (@ARGV){
	$pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
	$out = "o_lycpool_chrom$1";
	system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o $out"."vcf\n";
	$pm->finish;

}

$pm->wait_all_children;
```
Each chromosome scaffold is being processes seperately as a chrom*list

The variant data can be found in: 

We filtered the .vcf file with GATK (version ##), keeping only those (bases or snps)? with mapping quality>30, depth>1350, and bias scores less than +/-3. This was done with [VarFiltFork2.pl](VarFiltFork2.pl)
```pl
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
```

Allele depths were extracted from the filtered .vcf files. INDELS and multialleleic data were dropped here too. 
[AD.sh](AD.sh)
```sh
#!/usr/bin/bash
#
# extract allele depth AD from biallelic SNPs that passed filtering
#

for f in fff*vcf
do
        echo "Processing $f"
        out="$(echo $f | sed -e 's/vcf/txt/')"
        echo "Output is ad1_$out"
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\1/g' > ad1_$out
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\2/g' > ad2_$out
done
```
output:

We also grapped the SNP information with [SNP.sh](SNP.sh)
```sh
#!/usr/bin/bash
#
# extract alleles from biallelic SNPs that passed filtering
#

for f in fff*vcf
do
	echo "Processing $f"
	out="$(echo $f | sed -e 's/vcf/txt/')"
	echo "Output is snps_$out"
	grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | cut -f 4,5 > beastfilteredsnps_$out 
done
```
# Population Genetic Structure













# This is down to where I thoughtfully went through it all. Everything below is just formatted. 







# Counting Softclipping for real dataset
```sh
for f in max_chromV2ad1_beastfiltered_fff_o_lycpool_chrom*.fasta; do
    chrom=$(echo $f | grep -oP '\d+(?=\.fasta)')
    nchar=$(grep -v ">" $f | head -1 | tr -d '\n' | wc -c)
    size=$(du -sh $f | cut -f1)
    echo -e "chrom${chrom}\t${nchar}\t${size}"
done | sort -t'm' -k2 -n
```

# Mapdamage on Real Dataset 
```pl
##MapDamFork.pl
#!/usr/bin/perl
#run mapdamage

use Parallel::ForkManager;
my $max = 30;
my $pm = Parallel::ForkManager->new($max);

foreach $file (@ARGV){
        $pm->start and next; ## fork
        $file =~ m/_([a-zA-Z0-9]+)\.bam/;
        $id = $1;
        $out = "output_$id"."_filteredMapDamage";
        print "mapDamage -i $file -r ~/../gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta -d $out\n";
        system "mapDamage -i $file -r ~/../gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta -d $out\n";
        $pm->finish;
}
$pm->wait_all_children;
```

Run by SubMapDamFork.sh
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=xercesmapdamagefiltered
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load perl
module load samtools
module load bwa

cd /scratch/general/nfs1/u6000989/LycLotis/AllBams

perl MapDamFork.pl *bam
```

RAN THROUGH ALL XERCES ANALYSIS FOR SOFTCLIPPING TO INFORM THIS PROJECT. XERCES DATA TURNED OUT TO BE SHIT, BUT THROUGH MAP DAMAGE WE FOUND THAT THE LOTIS AND OTHER ADNA DATA FROM THIS SET LOOKED GOOD. ZACH HAD ALREADY ALIGNED, VARIANT CALLED, ETC (EVERYTHING UP TO VC FOR THIS ACTUAL DATA SET). I'M TAKING OVER FOR EVERYTHING THROUGH VARIANT CALLING. BASICALLY MAKING THE TREE AND PUTTING THE ANCIENT SAMPLES INTO AN ALRADY EXISTING TREE IN AN UPCOMING PAPER. 

USING THE FRAMEWORK LAYED OUT IN ZACH'S LYC-ADMIXTURE MOZAIC GITHUB
# VARIANT CALLING
USING THE FRAMEWORK LAYED OUT IN ZACH'S LYC-ADMIXTURE MOZAIC GITHUB
DIDN'T HAVE THE INDEXED BAI FILES, NEEDED TO INDEX USING THIS SCRIPT: 

Index.sh
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=Indexing
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /uufs/chpc.utah.edu/common/home/u6047808/ZLycLotis/AllBams

for bam in *.bam
do 
echo "Indexing $bam..."
samtools index "$bam"
done
```
# Variant calling
Ran with VariantCall.sh
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=VC_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /uufs/chpc.utah.edu/common/home/u6047808/ZLycLotis/AllBams
perl /uufs/chpc.utah.edu/common/home/u6047808/ZLycLotis/AllBams/VariantCall.pl chrom*list

```

Which ran: VariantCall.pl
```pl
#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#

use Parallel::ForkManager;
my $max = 26;
my $pm = Parallel::ForkManager->new($max);

my $genome ="/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";

foreach $chrom (@ARGV){
        $pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
        $out = "o_lycpool_chrom$1";
        system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o $out"."vcf\n";
        $pm->finish;

}

$pm->wait_all_children;
```
output: o_lycpool_chrom#.vcf

# Filtering with GATK
Ran: VarFiltFork2.sh
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=VarFilt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /uufs/chpc.utah.edu/common/home/u6047808/ZLycLotis/AllBams

perl /uufs/chpc.utah.edu/common/home/u6047808/ZLycLotis/AllBams/VarFiltFork2.pl
```

Which Runs: VarFiltFork2.pl
```pl
#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#

use Parallel::ForkManager;
my $max = 26;
my $pm = Parallel::ForkManager->new($max);

my $genome ="/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacB$

foreach $chrom (@ARGV){
        $pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
        $out = "o_lycpool_chrom$1";
        system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD $
        $pm->finish;

}
$pm->wait_all_children;
```
# Allele depth files
Next, extract allele depths from filtered vcf files. Also drops indels and multiallelic data. This gives ad1 and ad2. Ran with: AD.sh

```sh
#!/usr/bin/bash
#
# extract allele depth AD from biallelic SNPs that passed filtering 
#

for f in fff*vcf
do
        echo "Processing $f"
        out="$(echo $f | sed -e 's/vcf/txt/')"
        echo "Output is ad1_$out"
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\1/g' > ad1_$out   
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\2/g' > ad2_$out
done
```










# SECOND RUN AFTER SCRATCH SPACE WAS WIPED

##All bam and bai index files are in /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis
# Indexing
Indexed bam files using:
Index.sh
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=Indexingbai
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu


module load samtools
## version 1.16
module load bcftools
## version 1.16

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

for bam in *.bam
do 
echo "Indexing $bam..."
samtools index "$bam"
done
```
# Filter ends for mapdamage
Needed to Filter ends for mapdamage. Ran using:
FilterFork.sh
```sh
#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=kingspeak
#SBATCH --job-name=filter
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load samtools

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis 

for bam in *.bam; do
  base=$(basename "$bam" .bam)
  out="${base}.filtered${minlen}.bam"


  samtools view -h "$bam" \
    | awk -v m="$minlen" '($1 ~ /^@/) || (length($10) >= 50)' \
    | samtools view -b -o "$out" -

  samtools index "$out"
done
```
tried using perl script but gave up. 

# Ran MapDamage on all populations using 
MapDamFork.pl
```pl
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
        print "mapDamage -i $file -r ~/../gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta -d $out\n";
        system "mapDamage -i $file -r ~/../gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta -d $out\n";
        $pm->finish;
}
$pm->wait_all_children;
```

Run by SubMapDamFork.sh
```sh
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=xercesmapdamagefiltered
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load perl
module load samtools
module load bwa

cd /scratch/general/nfs1/u6000989/LycLotis/AllBams

perl MapDamFork.pl *bam
```

# Variant Calling
Make a folder containing all bams. Had to make the folder with all the full paths to the bams. Left the index files (.bai) in the main directories 
```sh
ls -1 *.bam > bams
head bams
```
## Run Variant calling on all bams with 
VariantCall.pl
```pl
#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#

use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta";

foreach $chrom (@ARGV){
        $pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
        $out = "o_lycpool_chrom$1";
        system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o $out"."vcf\n";
        $pm->finish;

}

$pm->wait_all_children;
```
output: o_lycpool_chrom*.vcf
Submitted with:SubVariantCall.sh
Note: chromlist is a series of empty files that will serve as the scaffold names for all the chromosomes for all populations. For example, chrom10.list file says Scaffold_1639;HRSCAF_2219 inside 
chrom1.list-chrom23.list

```sh
#!/bin/sh
#SBATCH --time=13:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=VC_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load samtools
module load bcftools
## version 1.16

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

perl /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/VariantCall.pl chrom*list
```


# Filtered the vcf file with GATK 
version (4.1.4.1), keeping only those with mapping quality > 30, depth > 1350 and bias scores less than +- 3. Using:
VarFiltFork2.pl
```pl
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
```
output: fff_o_lycpool_chrom10.vcf and fff_o_lycpool_chrom10.vcf.gz.tbi

The above was submitted by SubVarFiltFork2.sh
```sh
#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=VarFilt
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu


module load samtools
module load bcftools
## version 1.16

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

perl /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/VarFiltFork2.pl *.vcf.gz
```

# Allele depth and SNP
Used: AD.sh
```sh
#!/usr/bin/bash
#
# extract allele depth AD from biallelic SNPs that passed filtering
#

for f in fff*vcf
do
        echo "Processing $f"
        out="$(echo $f | sed -e 's/vcf/txt/')"
        echo "Output is ad1_$out"
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\1/g' > ad1_$out
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\2/g' > ad2_$out
done
```
output: ad1_fff_o_lycpool_chrom*.txt and ad2_fff_o_lycpool_chrom*.txt for each chromosome 
This creates allele depth files for each allele (ad1* and ad2*) and chromosome, which I can use for downstream analyses.

I also grapped the SNP information (alleles) using:
SNP.sh
```sh
#!/usr/bin/bash
#
# extract alleles from biallelic SNPs that passed filtering
#

for f in fff*vcf
do
        echo "Processing $f"
        out="$(echo $f | sed -e 's/vcf/txt/')"
        echo "Output is snps_$out"
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | cut -f 4,5 > snps_$out &
done
```
output: snps_fff_o_lycpool_chrom*.txt

# POPULATION GENETIC STRUCTURE FINALLY BITCH
Need to get this working fully: 
PCA and FST in R on chpc

getwd()
setwd()

library(data.table)

a1f<-list.files(pattern="ad1_fff")
a1f<-a1f[1:23]
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_fff_lycSpecPool_chrom","",a1f)
chrom<-gsub(".txt","",temp)
reps<-grep(pattern="rep",ids[,1])



pdf("WG_LG_PCAsV2.pdf",width=7,height=10.5)
par(mfrow=c(3,2))
par(mar=c(4.5,5.5,2.5,1.5))
cl<-1.3;ca<-1.1;cm<-1.4
for(i in 1:N){
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n<-a1+a2
	p<-a2/(a1+a2) ## non-ref

	p[is.na(p)]<-0.001
	## pca
	pc<-prcomp(t(p[,-c(reps,17)]),center=TRUE,scale=FALSE) ## drop replicates and MEN12
	o<-summary(pc)
	pct<-round(o$importance[2,1:3] * 100,1)

	plot(pc$x[,1],pc$x[,2],pch=rep(15:20,each=7),col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC2 (",pct[2],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("Chromosome ",chrom[i],sep=""),cex.main=cm)
	text(pc$x[,1],pc$x[,2],ids[-c(reps,17),1],cex=.7)
	xa<-min(pc$x[,1]) * .9
	ya<-pc$x[22,2] * .9
	#if(ya > 0){
	#	legend(xa,ya,ids[,1],pch=rep(15:20,each=7),col=ids[,2],ncol=3,cex=.6)
	#}
	plot(pc$x[,1],pc$x[,3],pch=rep(15:20,each=7),col=ids[-c(reps,17),2],xlab=paste("PC1 (",pct[1],")"),ylab=paste("PC3 (",pct[3],")"),cex.lab=cl,cex.axis=ca)
	title(main=paste("Chromosome ",chrom[i],sep=""),cex.main=cm)
	text(pc$x[,1],pc$x[,3],ids[-c(reps,17),1],cex=.7)

}
dev.off()
##view file in directory in chpc with xdg-open

# FST
Not fully working yet
```r
P<-vector("list",23)
n<-vector("list",23)
H<-vector("list",23)
for(i in 1:N){
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n[[i]]<-a1+a2
	P[[i]]<-a2/(a1+a2) ## non-ref
	H[[i]] <- 2 * P[[i]] * (1 - P[[i]])

}
```
## mn P
```
reps<-grep(pattern="rep",ids[,1])
mnP<-vector("list",23)
Nas<-vector("list",23)
for(i in 1:N){
	mnP[[i]]<-apply(P[[i]][,-reps],1,mean,na.rm=TRUE)
	Nas[[i]]<-apply(is.na(P[[i]][,-reps]),1,sum)

}
#retain 0.02 percent with the appropriate conditions
keepSNPs<-vector("list",23)
prop<-0.0002
for(i in 1:N){
	xx<-which(mnP[[i]] > 0.01 & mnP[[i]] < 0.99 & Nas[[i]]==0)
	keepSNPs[[i]]<-sort(sample(xx,floor(length(mnP[[i]])*prop),replace=FALSE))
}
###this gave me an error code so I tried this to check the number of eligible snps:
for(i in 1:N){
  xx <- which(mnP[[i]] > 0.01 & mnP[[i]] < 0.99 & Nas[[i]] == 0)
  cat("chrom", i, "eligible SNPs:", length(xx), "\n")
}


keepSNPs <- vector("list", 23)
prop <- 0.0002
for(i in 1:N){
    xx <- which(mnP[[i]] > 0.01 & mnP[[i]] < 0.99 & Nas[[i]] == 0)
    if(length(xx) > 0){ n_keep <- floor(length(xx) * prop)
        if(n_keep > 0){keepSNPs[[i]] <- sort(sample(xx, n_keep, replace = FALSE))
        } else {keepSNPs[[i]] <- integer(0)}
 } else {keepSNPs[[i]] <- integer(0) }
}
```
##sampled from filtered snps, not empy vectors, scale proportion based on length
```
for(i in 1:N){
	out<-paste("keepSNPs_chrom",chrom[i],sep="")
	write.table(keepSNPs[[i]],file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

save(list=ls(),file="fst.rdatV2.1")

##compute all of the pairwise Fst by pair and chromosome (not for each SNP)
Npop<-27
Nx<-(Npop*(Npop-1))/2
fstGw<-matrix(NA,nrow=Nx,ncol=23)
fst90<-matrix(NA,nrow=Nx,ncol=23)
x<-1
for(i in 1:(Npop-1)){
	for(j in (i+1):Npop){
		for(k in 1:23){
			keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
			pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
			Ht<-2*pbar*(1-pbar)
			Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
			
			fstGw[x,k]<-mean(Ht-Hs)/mean(Ht)
			fst<-(Ht-Hs)/Ht
			fst90[x,k]<-quantile(fst,.9,na.rm=TRUE)
		}
		x<-x+1
	}
}

##plot of mean and 90 quantile per chromsome for each pair
mord<-c(1,12,18:23,2:11,13:17)
pdf("FstLycSpcV2.1.pdf",width=8,height=10)
par(mfrow=c(5,3))
par(mar=c(4,5,2.5,.5))
x<-1
for(i in 1:(Npop-1)){
	for(j in (i+1):Npop){
		plot(fstGw[x,mord],pch=19,ylim=c(0,1),xlab="Chromosome",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=.9)
		segments(x0=1:24,y0=fstGw[x,mord],x1=1:24,y1=fst90[x,mord])
		title(main=paste(ids[i,1]," x ",ids[j,1],sep=""),cex.main=1.2)
		mn<-round(mean(fstGw[x,mord]),2)
		text(5,.9,mn)
		x<-x+1
	}
}
dev.off()

##simple plot of Fst continuum, this includes replicates, just sorted median (across chroms) Fst
plot(sort(apply(fstGw,1,median)),pch=19)
##pretty sure the uptick at the far right is Alaskc and France
##same thing just Z
plot(sort(fstGw[,24]),pch=19) ## very similar and against each other
plot(apply(fstGw,1,median),fstGw[,24],pch=19,xlab="Genome median",ylab="Z chrom")
abline(a=0,b=1)
mostly correlated, Z higher

### window examples
populations are from the ID.txt file and designated by i<-1 and j<-27

##ABM20 and YG20

mord<-c(1,12,18:23,2:11,13:17)
pdf("FstWinsABM20_YG20V2.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-1;j<-27
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

##MEN and SBW18
pdf("FstWinsMEN12_SBW18.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-16;j<-19
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

```
Lotis and SHC11
```
pdf("FstWinsLOTIS_SHC11V2.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-14;j<-20
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

########### Corrected code and the errors: Ht and Hs are already computed and have 
###length =length(keep). But Ht[keep] and Hs[keep] was indexing a short vector and 
##producing NAs. Here is the corrected code (I think)

pdf("FstWinsLOTIS_SHC11V2.pdf", width=8, height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i <- 14; j <- 20
for(k in mord){
  keep <- which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
  Nw <- floor(length(keep)/200)
  if(Nw < 1) next
  win <- rep(1:Nw, each=200)
  pbar <- (P[[k]][keep,i] + P[[k]][keep,j]) / 2
  Ht   <- 2*pbar*(1-pbar)
  Hs   <- P[[k]][keep,i]*(1-P[[k]][keep,i]) + P[[k]][keep,j]*(1-P[[k]][keep,j])
  Num <- tapply(Ht[1:(Nw*200)] - Hs[1:(Nw*200)], win, mean, na.rm=TRUE)
  Den <- tapply(Ht[1:(Nw*200)], win, mean, na.rm=TRUE)
  fst <- Num / Den
  plot(fst, xlab="SNP window", ylab=expression(F[ST]), cex.lab=1.3, cex.axis=1, ylim=c(0,1), type="l")
  title(main=paste("LG ", chrom[k], sep=""), cex.main=1.3)
}
dev.off()
Doing the same for my broken ones (will be labelled V2)
Lotis and YG20


pdf("FstWinsLOTIS_YG20V2.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-14;j<-27
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()
```
TBY51 and TBY11
```r
pdf("FstWinsTBY51_TBY11.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-24;j<-23
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()
```
# Version 2)- forgot to label as V2 in chpc
```
pdf("FstWinsTBY51_TBY11.pdf", width=8, height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i <- 24; j <- 23
for(k in mord){
  keep <- which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
  Nw <- floor(length(keep)/200)
  if(Nw < 1) next
  win <- rep(1:Nw, each=200)
  pbar <- (P[[k]][keep,i] + P[[k]][keep,j]) / 2
  Ht   <- 2*pbar*(1-pbar)
  Hs   <- P[[k]][keep,i]*(1-P[[k]][keep,i]) + P[[k]][keep,j]*(1-P[[k]][keep,j])
  Num <- tapply(Ht[1:(Nw*200)] - Hs[1:(Nw*200)], win, mean, na.rm=TRUE)
  Den <- tapply(Ht[1:(Nw*200)], win, mean, na.rm=TRUE)
  fst <- Num / Den
  plot(fst, xlab="SNP window", ylab=expression(F[ST]), cex.lab=1.3, cex.axis=1, ylim=c(0,1), type="l")
  title(main=paste("LG ", chrom[k], sep=""), cex.main=1.3)
}
dev.off()
BAT49 and BAT20

pdf("FstWinsBAT49_BAT20.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-3;j<-2
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

## V2
pdf("FstWinsBAT49_BAT20V2.pdf", width=8, height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i <- 2; j <- 2
for(k in mord){
  keep <- which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
  Nw <- floor(length(keep)/200)
  if(Nw < 1) next
  win <- rep(1:Nw, each=200)
  pbar <- (P[[k]][keep,i] + P[[k]][keep,j]) / 2
  Ht   <- 2*pbar*(1-pbar)
  Hs   <- P[[k]][keep,i]*(1-P[[k]][keep,i]) + P[[k]][keep,j]*(1-P[[k]][keep,j])
  Num <- tapply(Ht[1:(Nw*200)] - Hs[1:(Nw*200)], win, mean, na.rm=TRUE)
  Den <- tapply(Ht[1:(Nw*200)], win, mean, na.rm=TRUE)
  fst <- Num / Den
  plot(fst, xlab="SNP window", ylab=expression(F[ST]), cex.lab=1.3, cex.axis=1, ylim=c(0,1), type="l")
  title(main=paste("LG ", chrom[k], sep=""), cex.main=1.3)
}
dev.off()

###GNP17 and HNV17
pdf("FstWinsGNP17_HNV17.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-11;j<-13
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()

####GNP17 and Lotis


pdf("FstWinsGNP17_LOTIS.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-11;j<-14
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()






pdf("FstWinsMR20_CP19.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-17;j<-9
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()




pdf("FstWinsBCR17_BTB17.pdf",width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(4.5,5,2.5,.5))
i<-4;j<-7
for(k in mord){
	keep<-which(n[[k]][,i] > 10 & n[[k]][,j] > 10)
	Nw<-floor(length(keep)/200)
	win<-rep(1:Nw,each=200)
	pbar<-(P[[k]][keep,i]+P[[k]][keep,j])/2
	Ht<-2*pbar*(1-pbar)
	Hs<-P[[k]][keep,i] * (1-P[[k]][keep,i]) + P[[k]][keep,j] * (1-P[[k]][keep,j])
	Num<-tapply(X=Ht[keep][1:(Nw*200)]-Hs[keep][1:(Nw*200)],INDEX=win,mean)
	Den<-tapply(X=Ht[keep][1:(Nw*200)],INDEX=win,mean)
	plot(Num/Den,xlab="SNP window",ylab=expression(F[ST]),cex.lab=1.3,cex.axis=1,ylim=c(0,1),type='l')
	title(main=paste("LG ",chrom[k],sep=""),cex.main=1.3)
}
dev.off()
```



## Compute the invariant sites
## NOT WORKING YET

ComputeInvariant.R
determines the number of invariant bases to add to the beast XML file. Reads in base counts (A, C, G, T) by scaffold from the genome 
get this with: perl countBases.pl > baseCounts.txt, from /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/
cnts<-read.table("baseCounts.txt",header=FALSE)
```sh
## get big scaffolds, top 23
totals<-apply(cnts[,-1],1,sum)
rev(sort(totals))[1:23]
# [1] 31424484 26391115 25798047 25285067 25079126 25067872 22431864 21864934
# [9] 21544644 21226148 20022703 19681008 18787320 18503979 18445203 17279337
#[17] 17086116 16627315 16626411 16602730 15714686 13063727  9211676

chr<-which(totals >= 9211676)

## get total A, C, G, T
bcnt<-apply(cnts[chr,-1],2,sum)

## multiply by 0.00035 to match subsetting for SNP data
prop<-0.00035 ##Change this depending on what was computed
sbcnt<-floor(bcnt*prop)

## get SNP bases
## perl countBases.pl > snpCounts.txt
snps<-read.table("snpCounts.txt",header=FALSE)
## average across Lycaeides
snpCnts<-floor(apply(snps[-13,-1],2,mean))
invar<-sbcnt-snpCnts
invar
#   A    C    G    T 
#50190 28322 28284 50114
These are zachs numbers

```






# Input files for BEAST 
Generated with:
mkBeastDat.R
```r
library(data.table)

a1f<-list.files(pattern="ad1_fff")
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
asnp<-a1f
asnp<-gsub("ad1","snps",asnp)
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_fff_lycSpecPool_chrom","",a1f)
chrom<-gsub(".txt","",temp)

SSeq<-vector("list",27)

##this is where it breaks becasue snp chrom 9 doesnt exist?

for(i in 1:N){
        SSeq[[i]]<-vector("list",length(a1f))
}

for(i in 1:N){
        cat(i,"\n")
        out<-paste("max_chrom",chrom[i],".fasta",sep="")
        a1<-as.matrix(fread(a1f[i],header=F))
        a2<-as.matrix(fread(a2f[i],header=F))
        n<-a1+a2
        p<-a2/(a1+a2) ## non-ref
        p[n < 5]<-NA
        J<-dim(p)[2]
        L<-dim(p)[1]
        snps<-as.data.frame(fread(asnp[i],header=FALSE))
        for(j in 1:J){
                nx<-as.numeric(p[,j] > .5) + 1
                ss<-rep("N",L)
                jx<-which(is.na(p[,j])==FALSE)
                for(l in jx){
                        ss[l]<-snps[l,nx[l]]
                }
                SS1<-paste(ss,collapse="")
                #SSeq[[j]][[i]]<-paste(ss,collapse="")
                cat(">",ids[j,1],"\n",file=out,append=TRUE,sep="")
                cat(SS1,"\n",file=out,append=TRUE,sep="")
        }
        
}
```
output max_chromad1_fff_o_lycpool_chrom*.fasta
for round 2 output is max_chromad1V2_fff_o_lycpool_chrom*.fasta  . the output here should be ATGC by population
Did not need to use SubAlign.pl to drop rep pop samples

## Used mkNumericFasta.pl to convert alignments to numeric 
mkNumericFasta.pl
```pl
#!/usr/bin/perl

## this is to figure out which SNPs are variable in the fasta

foreach $i (1..23){
        system "grep -v \"^>\" sub_max_chromad1_fff_o_lycpool_chrom$i.fasta | perl -p -i -e 'tr/ACGTN/12345/' | sed 's/./& /g' > text_max_chrom$i.fasta\n";
        #system "grep -v \"^>\" sub_chrom$i.fasta | perl -p -i -e 'tr/ACGT/1234/' | sed 's/./& /g' > text_chrom$i.fasta\n";


###updated used: max_chromV2ad1_fff_o_lycpool_chrom$i.fasta
foreach $i (1..23){ 

	system "grep -v \"^>\" max_chromV2ad1_beastfiltered_fff_o_lycpool_chrom$i.fasta | perl -pe 'tr/ACGTN/12345/' | sed 's/./& /g' > sub_max_chromV2ad1_beastfiltered_fff_o_lycpool_chrom$i.fasta\n"; 
 
```
output: sub_max_chromV2ad1_fff_o_lycpool_chrom*.fasta
output is now numeric 
(plain sub_max_chrom is AGTC format, not numeric)

# Filter snps
Used GetSNPSubstMax.R to generate a filtered and reduced set of SNPs. This needed to be found in Zachs directory /uufs/chpc.utah.edu/common/home/gompert-group5/projects/LycAdmix/Beast
GetSNPSubstMax.R
```r
##identify SNPs from the alignments that are variable in the alignment
##subsample these

library(data.table)

miss<-vector("list",23)
for(i in 1:23){
        ifile<-paste("sub_max_chromV2ad1_fff_o_lycpool_chrom",i,".fasta",sep="")
        dat<-fread(ifile,header=FALSE)
        miss[[i]]<-apply(dat==5,2,mean)

}
## retain 0.035 percent with the appropriate conditions
keepSNPs<-vector("list",23)
prop<-0.00035
for(i in 1:23){
        xx<-which(miss[[i]]==0) ## no missing data
        keepSNPs[[i]]<-sort(sample(xx,floor(length(miss[[i]])*prop),replace=FALSE))
}
for(i in 1:23){
        out<-paste("keepSNPs_max_chrom",i,sep="")
        write.table(keepSNPs[[i]],file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
save(list=ls(),file="snps_max.rdat")
```
output: snps_max.rdat AND keepSNPS_max_chrom*
Then created combined across chromosomes fasta alignment filent file with just subset of SNPs with SubSetFasta.pl and converted to .nex file for beauti

SubSetFasta.pl
```pl

#!/usr/bin/perl
#
# this subsets and concatenates a set of SNPs from fasta
foreach $i (1..23){
	open(IN,"keepSNPs_max_chrom$i") or die "failed to open snps file $i\n";
	#open(IN,"keepSNPs_chrom$i") or die "failed to open snps file $i\n";
	$j = 0;
	while(<IN>){
		chomp;
		push (@{$snps[$i]},$_);
	}
	close(IN);
}
open(OUT, "> lyc_genomemax.fasta") or die "failed to write\n";
%seq;
foreach $i (1..23){
	open(IN,"sub_max_chrom$i.fasta") or die "failed to open snps file $i\n";
	while(<IN>){
		chomp;
		if(m/^>(\S+)/){
			$id = $1;
			if($i == 1){
				@{$seq{$id}} = ();
			}
		} else {
			foreach $snp (@{$snps[$i]}){
				$c = substr $_,$snp-1, 1;
				unless(length($c)==1){
					print "$c\n";
				}
				push(@{$seq{$id}}, $c);
			}
		}
	}
	close(IN);
}
foreach $pop (sort keys %seq){
	$str = join("",@{$seq{$pop}});
	unless($pop =~ m/rep/){
		$pop =~ s/Lyc-//;
		$pop =~ s/\d+//;
		print OUT ">$pop\n";
		print OUT "$str\n";
	}
}
close(OUT);
```
In terminal: seqmagick convert --output-format nexus --alphabet dna lyc_genomemax.fasta lyc_genomemax.nex
output: lyc_genomemax.fasta and .nexx

# Beauti model parameters
```
ml beast
beauti
```
used lyc_genomemax.nex as input for beauti (import alignment from partitions)- run inclding BAT49 and TBY51

didn't use tip dates
Site model- 
	gamma. Gamma category count: 4. Shape= 1
	Proportion invariant- 0
	Subst Model- GTR
	All Rates (AC, AG, GT, etc.) = 1 and are estimates. 
	Rate CT is not an estimate and = 1
	frequencies = estimated
Clock Model
	Random local clock
	check scaling....?
	mean clock rate = .0029
	
Priors:
For these all estimate boxes are left empty (in code they will = false)
	Tree- coalescent bayseian skyline (BSP)
	RRate changes- Poisson 
	clockrates- gamma w/ 1-> 1e-9, dimensions -44, apha= .05, beta=10, shapescale, offset=0
	freqparameters- dirichlet 4, 4, 4, 4, 1, dimensions 0 -> 1, alpha = 4,4,4,4, sum=1, offset=0
	gammaShapes- exponential, dim= 1 -> inf, mean= 1
	rateAC- gamma,  dim= 1 -> inf, alpha= .05, beta=10, offset= 0, mode=shapescale?
	rateAG- gamma,  dim= 1 -> inf, alpha= .05, beta=20, offset= 0, mode=shapescale?
	rateAT- gamma,  dim= 1 -> inf, alpha= .05, beta=10, offset= 0, mode=shapescale?
	rateCG- gamma,  dim= 1 -> inf, alpha= .05, beta=10, offset= 0, mode=shapescale?
	rateGT- gamma,  dim= 1 -> inf, alpha= .05, beta=20, offset= 0, mode=shapescale?
	TMRCAALL.prior-all taxon, mean = 2.125, sigma=.53, offset= 0, monophyletic- true, normal distribution 
	TMRCALyc.prior- all taxon but outgrou MEN, mean=1.29, sigma= .53, offset= 0, monophyletic= true, normal distribution 
MCMC
Coupled MCMC
	Chains:5
resample every 1000
delta temp: .025
optimise= true
optimise delay: 100?
target: wtf does this mean .234?
chainlength: 500000000
Store every 25000
Pre burnin 400
Number initialization attempts: 800

last run saved as lyc_wgs_max_ranlc4.xml- with both BAT49 and TBY51
Note: if i were to do this run again, need to update the means and sigmas for the priors on lyc and lycall. 
Also, in the .xml file, need to put the invariant site counts. For now I'm using Zachs until I figure out my own:  
```
<data id='lyc_genomemax' spec='FilteredAlignment' filter='-' data='@lyc_genomemaxOrig' constantSiteWeights='50190 28322 28284 50114'>
    </data>
```
Be sure to change file name too

# CASTER 
For cASTER to work, run: ml gcc/13.3.0
Caster and Beast use same input file for mkBeastDat.R, so these already exist, but to be sure and to match with zachs, ran with sub_max instead of max

output sub_max_chromad1_fff_o_lycpool_chrom*.fasta
Zach also used SubAlign.pl here, but I do not need to because I have no replicate pops

From here I am running RunCasterSite.pl 
```pl
#!/usr/bin/perl
#

foreach $i (1..23){

        system "/uufs/chpc.utah.edu/common/home/u6047808/bin/ASTER-Linux/bin/caster-site -i sub_max_chromad1_fff_o_lycpool_chrom$i.fasta -o CASTcout_max_$i --root MEN --thread 24\n";
#       system "ASTER-Linux/bin/caster-site -i sub_CASTchromad1_fff_o_lycpool_chrom$i.fasta -o CASTcout_max_$i --root MEN --thread 24\n";
}

##output: count_max
### Then using ape (version ?) to plot the 23 trees ## while rotating around nodes to maximixe visula similarity for comparison..... apparently 

###using plotTrees.R
#### the input from this is coming from step above

library(ape)
## ape version 5.8

pdf("casterTrees.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<-paste("CASTcout_max_",i,sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE)
	title(main=paste("Chrom.",i),cex.main=1.3)
}

dev.off()

pdf("casterTreesMax.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<-paste("CASTcout_max_",i,sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE,type="cladogram")
	title(main=paste("Chromosome",i),cex.main=1.3)
}

dev.off()
```

## trying to clean things up
```r
trees<-vector("list",23)
for(i in 1:23){
	inf<-paste("CASTcout_max_",i,sep="")
	trees[[i]]<-read.tree(inf)
}

ref_tree<-trees[[1]]
ref_order<-ref_tree$tip.label
## rotateConstr is part of ape
alntrees<-lapply(trees, function(tr) rotateConstr(tr, ref_order))
pdf("casterTreesMax.pdf",width=9,height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	plot.phylo(alntrees[[i]],cex=.7,use.edge.length=FALSE,type="cladogram")
	title(main=paste("Chromosome",i),cex.main=1.3)
}

dev.off()
```
## formated version of figure
```
pdf("fig_caster.pdf",width=8,height=10)
par(mfrow=c(4,6))
par(mar=c(1,1,1,1))
for(i in 1:23){
	plot.phylo(alntrees[[i]],cex=.7,use.edge.length=FALSE,type="cladogram")
	mtext(paste("Chr.",i),cex=1.1,line=-2,side=3,adj=.1)
}

dev.off()
```
should have output pdf's here

# Window based analyses with cASTER

im afraid hahahahahaaaaa

Now I have no idea what's going on I need somebody to hold my hand

## running anumber of sliding window analyses with CASTER. 
Zach says: For each set of 4 taxa (A, B, C and outgroup) I compute scores in 10 kb windows, then average over sets of 5 windows to plot 
normalized (sum to 1) scores across the genome. I have one sub-directory for each set of taxa. The basic commands look like this (executed from 
within the sub-directory with a mapping file, Sub* file and symbolic links to the sub_max* alignments)######

What did I do? Ran WinSubAlign.pl
```pl
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
```
AND
RunWindows.pl:
```pl
#!/usr/bin/perl
#
for $i (1..23){
        system "../wins/MASTERWORK/bin/slidingwindow SubABMxSINxTBY_sub_max_chrom$i.fasta MappingABMxSINxTBY.txt > winout$i.tsv\n";
}
```
The results were summarized with 
Win.Test.R

```R
## 10 kb windows HJxVExSIN
# A = HJ, B = VE, C = SIN
# intraspecific only? 

scAB<-vector("list",23)
scBC<-vector("list",23)
scAC<-vector("list",23)
for(ch in 1:23){
	inf<-paste("winout",ch,".tsv",sep="")
	dat<-read.table(inf,header=TRUE)
	L<-dim(dat)[1]
	NL<-L-5
	ab<-rep(NA,NL)
	ac<-rep(NA,NL)
	bc<-rep(NA,NL)
	for(i in 1:NL){
		ab[i]<-mean(dat$A.B[i:(i+4)])
		bc[i]<-mean(dat$B.C[i:(i+4)])
		ac[i]<-mean(dat$A.C[i:(i+4)])
		sc<-c(ab[i],bc[i],ac[i])
		ssc<-sum(sc)
		ab[i]<-ab[i]/ssc
		bc[i]<-bc[i]/ssc
		ac[i]<-ac[i]/ssc
	}
	scAB[[ch]]<-ab
	scAC[[ch]]<-ac
	scBC[[ch]]<-bc
}

mat<-rbind(unlist(scAB),unlist(scAC),unlist(scBC))

apply(mat,1,mean)
#[1] 0.6354200 0.2473369 0.1172431
# looks like HJ x SIN more than VE x SIN? could be real?


szs<-unlist(lapply(scAB,length))
ch<-rep(c(1:23),szs)

pdf("winHJxVExSIN.pdf",width=9,height=4)
par(mar=c(5,5,1,1))
xx<-barplot(mat,border=NA,axes=FALSE,xlab="Chromosome",ylab="Score")
mids<-tapply(X=xx,INDEX=ch,mean)
bnds<-tapply(X=xx,INDEX=ch,max)[-23]
abline(v=bnds)
axis(2)
axis(1,mids,c(1:22,"Z"))
dev.off()
```

TO LOOK AT THE BEAST TREE: 
java -jar ~/../gompert-group5/projects/LycAdmix/Beast/FigTree_v1.4.4/lib/figtree.jar Combined_consensusranlc4.trees OR TESTZ or Combinedtreesfinal.trees
# Checking coverage for TBY, BAT, and LOTIS with ad1 and ad2 files 
```r
R
be in correct working directory
library(data.table)
a1 <- as.matrix(fread("ad1_fff_o_lycpool_chrom1.txt", header=FALSE))
a2 <- as.matrix(fread("ad2_fff_o_lycpool_chrom1.txt", header=FALSE))
cov <- a1 + a2

##rows= snps
##colums= pops


nms <- read.table("bams", header=FALSE)
clean_nms <- gsub("dedup_|.bam", "", nms[,1])
colnames(cov) <- clean_nms

lotis_cov<- cov[, "/uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filt_LOTIS"]
bat20_cov<- cov[, "/uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filt_BAT20"]
bat49_cov<- cov[, "/uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filt_BAT49"]
tby11_cov<- cov[, "/uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filt_TBY11"]
tby51_cov<- cov[, "/uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filt_TBY51"]
shc11_cov<- cov[, "/uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filt_SHC11"]

summary(lotis_cov)
summary(bat49_cov)
summary(tby51_cov)
summary(bat20_cov)

```
finding a sweet spot between the two extremes from the ancient populations. Somewhere between 8-20x coverage and 500x coverage. Will then read in each ad file for each chromosome, index it, and filter for each population 
based on what I decided
```r
sum(lotis_cov > 20)
sum(lotis_cov > 20 & bat49_cov > 20)
sum(lotis_cov > 20   & lotis_cov < 500)

library(data.table)
a1 <- as.matrix(fread("ad1_fff_o_lycpool_chrom23.txt", header=FALSE))
a2 <- as.matrix(fread("ad2_fff_o_lycpool_chrom23.txt", header=FALSE))
snps <- as.matrix(fread("beastfilteredsnps_fff_o_lycpool_chrom23.txt", header=FALSE))
cov <- a1 + a2
min_pops <- 27
keep <- apply(cov, 1, function(x) sum(x >= 20 & x <= 500) >= min_pops)


##nrow(keep)		# total SNPs before filtering
length(keep)
sum(keep)       # SNPs kept after filtering
sum(!keep)		#SNPs removed

a1_filtered <- a1[keep,]
a2_filtered <- a2[keep,]
snps_filtered<- snps[keep,]

out1 <- "ad1_beastfilteredV2_fff_o_lycpool_chrom23.txt"
out2 <- "ad2_beastfilteredV2_fff_o_lycpool_chrom23.txt"
out3 <- "beastfilteredsnpsV2_fff_o_lycpool_chrom23.txt"


fwrite(as.data.table(a1_filtered), file=out1, sep="\t", col.names=FALSE)	
fwrite(as.data.table(a2_filtered), file=out2, sep="\t", col.names=FALSE)
fwrite(as.data.table(snps_filtered), file=out3, sep="\t", col.names=FALSE)

a1_original <- as.matrix(fread("ad1_fff_o_lycpool_chrom14.txt", header=FALSE))
a1_filtered_check <- as.matrix(fread("ad1_beastfilteredV2_fff_o_lycpool_chrom14.txt", header=FALSE))

snps_origional <- as.matrix(fread("beastfilteredsnps_fff_o_lycpool_chrom10.txt", header=FALSE))
snps_filtered_check <- as.matrix(fread("beastfilteredsnpsV2_fff_o_lycpool_chrom14.txt", header=FALSE))

nrow(a1_original)      
nrow(a1_filtered_check)

nrow(snps_origional)
nrow(snps_filtered_check)


head(fread("beastfilteredad1_fff_o_lycpool_chrom.txt", header=FALSE))
```
check how many SNPS total

wc -l beastfiltered*
24560


Needed to regenerate the SNP information from SNPs.sh. Still used my origional filtered set from round 1
where filtered the vcf file with GATK version (4.1.4.1), keeping only those with mapping quality > 30, depth > 1350 and bias scores less than +- 3

Using VarFiltFork2.pl. 
output was beastfilteredsnps_fff_o_lycpool_chrom_.txt
allele depth files are taken care of from code up above. 

NOTE: output has been switched to ad1 and 2_ beastfiltered_fff_o_lycpool_chrom


Need to also filter the beastfilteredsnp files. Running these through the same code as ad1 and ad2 files in R. Code 
will be below. first I copied the origional SNP files snps_fff_o_lycpool_chrom and names the set to be filtered beastfilteredsnps_fff_o_lycpool_chrom the new files will have that new corrected name beastfilteredsnps_fff_o_lycpool_chrom. 
```r
library(data.table)
a1 <- as.matrix(fread("fff_o_lycpool_chrom10.vcf", header=FALSE))
cov <- a1 
min_pops <- 27
keep <- apply(cov, 1, function(x) sum(x >= 20 & x <= 500) >= min_pops)


nrow(a1)      # total SNPs before filtering
sum(keep)       # SNPs kept after filtering
sum(!keep)		#SNPs removed

snps_filtered <- a1[keep,]

out1 <- "beastfilteredsnpsTEST_fff_o_lycpool_chrom10.txt"

fwrite(snps_filtered, file=out1, sep="\t", col.names=FALSE)	

a1_original <- as.matrix(fread("beastfilteredsnps_fff_o_lycpool_chrom10.txt", header=FALSE))
a1_filtered_check <- as.matrix(fread("beastfilteredsnpsTEST_fff_o_lycpool_chrom10.txt", header=FALSE))

nrow(a1_original)      
nrow(a1_filtered_check)

head(fread("beastfilteredad1_fff_o_lycpool_chrom2.txt", header=FALSE))
```

Input files for BEAST

```r
mkBeastDat.R

library(data.table)

a1f<-list.files(pattern="ad1_beastfiltered")
a1f<-a1f[1:23]
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)

asnp<-a1f
asnp<-list.files(pattern="beastfilteredsnpsV2")
asnp<-asnp[1:23]
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_fff_lycSpecPool_chrom","",a1f)
chrom<-gsub(".txt","",temp)

SSeq<-vector("list",27)
for(i in 1:N){
	SSeq[[i]]<-vector("list",length(a1f))
}

for(i in 1:N){
	cat(i,"\n")
	out<-paste("CAST_chrom_V2",chrom[i],".fasta",sep="")
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n<-a1+a2
	p<-a2/(a1+a2) ## non-ref
	p[n < 5]<-NA
	J<-dim(p)[2]
	L<-dim(p)[1]
	snps<-as.data.frame(fread(asnp[i],header=FALSE))
	for(j in 1:J){
		nx<-as.numeric(p[,j] > .5) + 1
		ss<-rep("N",L)
		jx<-which(is.na(p[,j])==FALSE)
		for(l in jx){
			ss[l]<-snps[l,nx[l]]
		}
		SS1<-paste(ss,collapse="")
		#SSeq[[j]][[i]]<-paste(ss,collapse="")
		cat(">",ids[j,1],"\n",file=out,append=TRUE,sep="")
		cat(SS1,"\n",file=out,append=TRUE,sep="")
	}
	
}
```
output: max_chromad1_beastfilteredV2_fff_o_lycpool_chrom   
had output: max_chromFinalad1_beastfilteredV2_fff_o_lycpool_chrom - but pretty sure it was updated to line above.   


# Generate Caster input files
using: mkCaster.R
```r
library(data.table)

a1f<-list.files(pattern="ad1_beastfiltered")
a1f<-a1f[1:23]
a2f<-a1f
a2f<-gsub("ad1","ad2",a2f)
asnp<-a1f
asnp<-list.files(pattern="beastfilteredsnpsV2")
asnp<-asnp[1:23]
N<-length(a1f)
ids<-read.table("IDs.txt",header=FALSE)
temp<-gsub("ad1_fff_lycSpecPool_chrom","",a1f)
chrom<-gsub(".txt","",temp)

SSeq<-vector("list",27)
for(i in 1:N){
	SSeq[[i]]<-vector("list",length(a1f))
}

for(i in 1:N){
	cat(i,"\n")
	out<-paste("CAST_chrom_V2",chrom[i],".fasta",sep="")
	a1<-as.matrix(fread(a1f[i],header=F))
	a2<-as.matrix(fread(a2f[i],header=F))
	n<-a1+a2
	p<-a2/(a1+a2) ## non-ref
	p[is.na(p)]<-0.001
	J<-dim(p)[2]
	L<-dim(p)[1]
	snps<-as.data.frame(fread(asnp[i],header=FALSE))
for(j in 1:J){
		nx<-rbinom(n=L,size=1,prob=p[,j]) + 1
		ss<-rep(NA,L)
		for(l in 1:L){
			ss[l]<-snps[l,nx[l]]
		}
		SS1<-paste(ss,collapse="")
		#SSeq[[j]][[i]]<-paste(ss,collapse="")
		cat(">",ids[j,1],"\n",file=out,append=TRUE,sep="")
		cat(SS1,"\n",file=out,append=TRUE,sep="")
	}
	
}

```
## check if this still true
-for this I just made my updated max_chrom_finalad1_beastfiltered files my sub_max files to avoid naming confusion and mislabelling later on. 
mv max_chrom_Finalad1_beastfilteredV2_fff_o_lycpool_chrom9.fasta sub_max_chrom9.fasta
didn't need to run SubAlign.pl 
Output is now  sub_max_chromV2ad1_beastfiltered_fff_o_lycpool_chrom*


## snp subsetting
```pl
mkNumericFasta.pl
#!/usr/bin/perl

## this is to figure out which SNPs are variable in the fasta

foreach $i (1..23){
	system "grep -v \"^>\" max_chrom$i.fasta | perl -p -i -e 'tr/ACGTN/12345/' | sed 's/./& /g' > text_max_chrom$i.fasta\n";
	#system "grep -v \"^>\" sub_chrom$i.fasta | perl -p -i -e 'tr/ACGT/1234/' | sed 's/./& /g' > text_chrom$i.fasta\n";

}

foreach $i (1..23){
    system "grep -v \"^>\" sub_max_chromV2ad1_beastfiltered_fff_o_lycpool_chrom$i.fasta | perl -pe 'tr/ACGTN/12345/' | sed 's/./& /g' > text_chrom$i.fasta\n";
```
output: text_chrom*.fasta

Next filtered the SNPs for Beast, doing .5 (half) instead of .0035 using GetSNPSubstMax.R
Identify SNPs from the alignments that are variable in the alignments and subsample these
```r
library(data.table)

miss<-vector("list",23)
for(i in 1:23){
        ifile<-paste("text_chrom",i,".fasta",sep="")
        dat<-fread(ifile,header=FALSE)
        miss[[i]]<-apply(dat==5,2,mean)

}
## retain 0.5 percent with the appropriate conditions
keepSNPs<-vector("list",23)
prop<-0.5
for(i in 1:23){
        xx<-which(miss[[i]]==0) ## no missing data
        keepSNPs[[i]]<-sort(sample(xx,floor(length(miss[[i]])*prop),replace=FALSE))
}
for(i in 1:23){
        out<-paste("keepSNPs_max_chrom",i,sep="")
        write.table(keepSNPs[[i]],file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
save(list=ls(),file="snps_maxV2.rdat")
```
output: snps_maxV2.rdat




### Combined fasta alignments with subset of snps 
Using SubSetFasta.pl:
```pl

#!/usr/bin/perl
#
# this subsets and concatenates a set of SNPs from fasta

foreach $i (1..23){
	open(IN,"keepSNPs_max_chrom$i") or die "failed to open snps file $i\n";
	#open(IN,"keepSNPs_chrom$i") or die "failed to open snps file $i\n";
	$j = 0;
	while(<IN>){
		chomp;
		push (@{$snps[$i]},$_);
	}
	close(IN);
}
open(OUT, ">lyc_genomemaxV2.fasta") or die "failed to write\n";
%seq;
foreach $i (1..23){
	open(IN,"sub_max_chrom$i.fasta") or die "failed to open snps file $i\n";
	while(<IN>){
		chomp;
		if(m/^>(\S+)/){
			$id = $1;
			if($i == 1){
				@{$seq{$id}} = ();
			}
		} else {
			foreach $snp (@{$snps[$i]}){
				$c = substr $_,$snp-1, 1;
				unless(length($c)==1){
					print "$c\n";
				}
				push(@{$seq{$id}}, $c);
			}
		}
	}
	close(IN);
}
foreach $pop (sort keys %seq){
	$str = join("",@{$seq{$pop}});
	unless($pop =~ m/rep/){
		$pop =~ s/Lyc-//;
		##$pop =~ s/\d+//;
		print OUT ">$pop\n";
		print OUT "$str\n";
	}
}
close(OUT);
```

seqmagick convert --output-format nexus --alphabet dna lyc_genomemaxV2.fasta lyc_genomemaxV2.nex

when I put the .nex file into beauti, I have 6135 bases after this filtering round
put into beauti and save as output: lyc_wgs_ranlc.....
put log files into logcombiner with output: combinedlogsFinal.log




## Counting the proportion of SNPs per chromosome before and after filtering to determine whether or not to cut BAT49 and TBY51
still need to cut out BAT49

##Pre-filtered allele depth files are under: ad1_fff_o_lycpool_chrom*.txt 
	pre filtered includes all pops, all chroms, filtered with GATK as a vcf. Does not including the strict constraints of 20-500x coverage
	pre filtered snp files under: snps_fff_o_lycpool_chrom*.txt
	(made 4/23/26, about 50-90M)

##Post-filtered ad files under: ad1_beastfilteredV2_fff_o_lycpool_chrom
	post filtered includes all pops, all chroms, with a strict filter of 20-500x coverage across all 27 pops. I have not filtered for lotis, bat49, or tby51 seperately yet. 
	post filtered snp files under: beastfilteredsnps_fff_o_lycpool..... - I fear these may have been overwritten. 
		In /noBAT49 I also have beastfilteredV2_fff_... that I believe is the pre filtered snp files, but run through more stringent filters. 


# checking coverage and proportion per each chrom
Get row counts for each chromosome file
```
wc -l ad1_beastfilteredV2_fff_o_lycpool_chrom*.txt
```
Total= 12280

Confirm a file has no header (first row should be all numbers)
```
head -1 ad1_beastfilteredV2_fff_o_lycpool_chrom1.txt
```
this gave me the number of SNPs per chromosome. 12,280 total snps over 23 chromosomes and across 27 pops. Lowest coverage was population 3, BAT49

SNP counts and proportions per chromosome 
All this code can be found in SNPsProportionCounts.sh
we know the total SNPs is 12280
```sh
total=12280

for f in ad1_beastfilteredV2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportions2.txt

cat snp_counts_proportions.txt
```
output: snp_counts_proportions.txt and second run was snp_counts_proportions2.txt

PROPORTION OF SNPS ON EACH CHROMOSOME

next going to caluclate the mean depth per population (column)
```r
for f in ad1_beastfilteredV2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    awk -v chrom="$chrom" '
    {
        for (i=1; i<=NF; i++) {
            sum[i] += $i
            count[i]++
        }
    }
    END {
        for (i=1; i<=NF; i++)
            print chrom, i, sum[i]/count[i]
    }' $f
done > mean_coverage_per_chrom_pop2.txt
```
Note for code: 
chrom=$(echo $f | grep -oP 'chrom\d+') - extracts chromosome label followed by one or more digits and prints only the match. P- perl style regex sorts ouput. n sorts numerically, k2 sorts second field (alphebetically and then numerically). t'm' sets M in chrom as the end or delimiter (cutoff)

Check
```
head mean_coverage_per_chrom_pop.txt
```
Output format: chrom  pop_index  mean_depth
output file: mean_coverage_per_chrom_pop*.txt - run 1 and 2 same so far

Just added population names to values (mean depths and snp counts) using my IDs.txt file. Created sample_names.txt to avoid overwriting my IDs.txt file again
```
awk '{print NR"\t"$1}' IDs.txt > sample_names.txt
##Check that it worked
cat sample_names.txt
##Now join the files 

awk 'NR==FNR{name[$1]=$2; next} {print $1, name[$2], $3}' \
  sample_names.txt mean_coverage_per_chrom_postfilter.txt \
  > mean_coverage_named.txt

mean_coverage_total_postfilter.txt
```
Output format: chrom  sample_name  mean_depth 
head mean_coverage_named.txt

Showing us BAT49 has the worst depth coverage (~30) conpared to lotis (~40-50) and TBY51~(40-50). 
To check coverage for each population: grep BAT49 mean_coverage_named.txt

## Chromosome by population matrix (##rows= chrom, columns= pops)
### Now on the post filtered
```
total=12280
for f in beastfilteredsnpsV2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportionspost.txt

cat snp_counts_proportionspost.txt
```
output: snp_counts_proportionspost.txt and second run was snp_counts_proportions2.txt

### PROPORTION OF SNPS ON EACH CHROMOSOME
next going to caluclate the mean depth per population (column)
```
for f in ad1_beastfilteredV2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    awk -v chrom="$chrom" '
    {
        for (i=1; i<=NF; i++) {
            sum[i] += $i
            count[i]++
        }
    }
    END {
        for (i=1; i<=NF; i++)
            print chrom, i, sum[i]/count[i]
    }' $f
done > mean_coverage_per_chrom_pop2.txt
```
head mean_coverage_per_chrom_pop.txt
Output format: chrom  pop_index  mean_depth
output file: mean_coverage_per_chrom_pop*.txt - run 1 and 2 same so far

Just added population names to values (mean depths and snp counts) using my IDs.txt file. Created sample_names.txt to avoid overwriting my IDs.txt file again
```
awk '{print NR"\t"$1}' IDs.txt > sample_names.txt
##Check that it worked
cat sample_names.txt
##Now join the files 

awk 'NR==FNR{name[$1]=$2; next} {print $1, name[$2], $3}' \
  sample_names.txt mean_coverage_per_chrom_pop.txt \
  > mean_coverage_named.txt
```
Output: chrom  sample_name  mean_depth 
head mean_coverage_named.txt
Showing us BAT49 has the worst depth coverage (~30) conpared to lotis (~40-50) and TBY51~(40-50). to check coverage for each population: 
```
grep BAT49 mean_coverage_named.txt
```
##Chromosome by population matrix (##rows= chrom, columns= pops)
###Fixed this by combining ad 1 and ad2
```
# Confirm both sets exist
ls ad1_beastfilteredV2_fff_o_lycpool_chrom*.txt | wc -l
ls ad2_beastfilteredV2_fff_o_lycpool_chrom*.txt | wc -l

# Confirm row counts match between ad1 and ad2 for a few chromosomes
wc -l ad1_beastfilteredV2_fff_o_lycpool_chrom1.txt \
       ad2_beastfilteredV2_fff_o_lycpool_chrom1.txt \
       ad1_beastfilteredV2_fff_o_lycpool_chrom10.txt \
       ad2_beastfilteredV2_fff_o_lycpool_chrom10.txt

##these should have the same number as their counterparts

###Mean total coverage per chrom per sample
for f in ad1_beastfilteredV2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    f2=$(echo $f | sed 's/ad1/ad2/')
    paste $f $f2 | awk -v chrom="$chrom" -v ncol=27 '{
        for (i=1; i<=ncol; i++) {
            sum[i] += $i + $(i+ncol)
            count[i]++
        }
    }
    END {
        for (i=1; i<=ncol; i++)
            print chrom, i, sum[i]/count[i]
    }'
done > mean_coverage_total_postfilter.txt

head mean_coverage_total_postfilter.txt
```

Adding sample names 
```
awk '{print NR"\t"$1}' IDs.txt > sample_names.txt
##Check that it worked
cat sample_names.txt
##Now join the files 

awk 'NR==FNR{name[$1]=$2; next} {print $1, name[$2], $3}' \
sample_names.txt mean_coverage_total_postfilter.txt \
mean_coverage_total_postfilter.txt

head mean_coverage_total_postfilter.txt
```
### Build coverage matrix 

```
echo -e "chrom\t$(awk '{print $1}' IDs.txt | tr '\n' '\t' | sed 's/\t$//')" > coverage_matrix_postfilter.txt

# Matrix
awk '
{
    val[$1][$2] = $3
    pops[$2] = 1
}
END {
    n = asorti(pops, poplist)
    for (c = 1; c <= 23; c++) {
        chrom = "chrom" c
        printf chrom
        for (i = 1; i <= n; i++)
            printf "\t" val[chrom][poplist[i]]
        printf "\n"
    }
}' mean_coverage_total_postfilter.txt >> coverage_matrix_postfilter.txt

cat coverage_matrix_postfilter.txt
```

###prefiltered
```
total_pre=$(wc -l ad1_fff_o_lycpool_chrom*.txt | grep total | awk '{print $1}')
echo $total_pre

for f in ad1_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    f2=$(echo $f | sed 's/ad1/ad2/')
    paste $f $f2 | awk -v chrom="$chrom" -v ncol=27 '{
        for (i=1; i<=ncol; i++) {
            sum[i] += $i + $(i+ncol)
            count[i]++
        }
    }
    END {
        for (i=1; i<=ncol; i++)
            print chrom, i, sum[i]/count[i]
    }'
done > mean_coverage_total_prefilter.txt

##adding pop names
awk '{print NR"\t"$1}' IDs.txt > sample_names.txt
##Check that it worked
cat sample_names.txt
##Now join the files 

awk 'NR==FNR{name[$1]=$2; next} {print $1, name[$2], $3}' \
sample_names.txt mean_coverage_total_prefilter.txt \
mean_coverage_total_prefilter.txt

head mean_coverage_total_prefilter.txt

##matrix
awk '
{
    val[$1][$2] = $3
    pops[$2] = 1
}
END {
    n = asorti(pops, poplist)
    for (c = 1; c <= 23; c++) {
        chrom = "chrom" c
        printf chrom
        for (i = 1; i <= n; i++)
            printf "\t" val[chrom][poplist[i]]
        printf "\n"
    }
}' mean_coverage_total_prefilter.txt >> coverage_matrix_prefilter.txt

cat coverage_matrix_prefilter.txt
```
## Post and pre filtering SnP counts and proportions 
```
###find the totals: 
For pre filtered: wc -l beastfilteredsnps_fff_o_lycpool_chrom*.txt
15387914

total_pre=15387914

for f in beastfilteredsnps_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total_pre" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportions_prefilter.txt

cat snp_counts_proportions_prefilter.txt
```

## Post filtereed: wc -l beastfilteredsnpsV2_fff_o_lycpool_chrom*.txt
12280
Find props
```
total_post=12280

for f in beastfilteredsnpsV2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total_post" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportions_postfilter.txt

cat snp_counts_proportions_postfilter.txt



get rid of BAT49
set a min for lotis and tomboy (10x) and an AVERAGE min for the rest- mean of 20x and below 500
don't include lotis and tomboy in average- maybe
check for outside of the bell (multicopy)- knock out top few percent 
still need to do thinning of snps beast
```

##checking contents of origional files 

#Confirm column order in original ad1 file
```
head -1 ad1_fff_o_lycpool_chrom1.txt

# Check original SNP file
head -1 beastfilteredsnps_fff_o_lycpool_chrom1.txt

# Confirm your VCF header column order
grep "^#CHROM" fff_o_lycpool_chrom1.vcf | tr '\t' '\n' | nl | head -15

## made new directory for no BAT49
# Remove column 3 from all original ad1, ad2, and SNP files

for f in ../ad1_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    cut -f1-2,4- $f > ad1_fff_o_lycpool_${chrom}_noBAT49.txt
done

for f in ../ad2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    cut -f1-2,4- $f > ad2_fff_o_lycpool_${chrom}_noBAT49.txt
done

for f in ../beastfilteredsnps_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    cut -f1-2,4- $f > beastfilteredsnps_fff_o_lycpool_${chrom}_noBAT49.txt
done
```
### Verify that it worked
#Check columns - should be 26
```
awk '{print NF}' ad1_fff_o_lycpool_chrom1_noBAT49.txt | sort -u
## Check SNP file is still 2 columns
awk '{print NF}' beastfilteredsnps_fff_o_lycpool_chrom1_noBAT49.txt | sort -u
## Rows could match originals 
wc -l ad1_fff_o_lycpool_chrom1_noBAT49.txt ../ad1_fff_o_lycpool_chrom1.txt

# Check original column order from VCF
grep "^#CHROM" ../fff_o_lycpool_chrom1.vcf | tr '\t' '\n' | nl | grep -E "LOTIS|TBY51|TBY11"

# The VCF and .txt files have different headers, so need to convert .vcf to .txt. Subtract metacolumns up to 9 (format)
# to get to the populations. 

## To check that the positions of the populations, run:

grep "^#CHROM" ../fff_o_lycpool_chrom1.vcf | tr '\t' '\n' | \
    tail -n +10 | \
    grep -v "BAT49" | \
    nl
```
## R code to filter heavily on aDNA, after getting rid of BAT49
Lotis and TBY51, min coverage of 10
```
library(data.table)

lotis_col <- 13
tby51_col <- 23
other_cols <- setdiff(1:26, c(lotis_col, tby51_col))

for (chrom in 1:23) {
  cat("Processing chrom", chrom, "\n")
  
  a1 <- as.matrix(fread(paste0("ad1_fff_o_lycpool_chrom", chrom, "_noBAT49.txt"), header=FALSE))
  a2 <- as.matrix(fread(paste0("ad2_fff_o_lycpool_chrom", chrom, "_noBAT49.txt"), header=FALSE))
  snps <- as.matrix(fread(paste0("beastfilteredsnps_fff_o_lycpool_chrom", chrom, "_noBAT49.txt"), header=FALSE))
  
  cov <- a1 + a2
  
  keep <- apply(cov, 1, function(x) {
    lotis_ok <- x[lotis_col] >= 10 & x[lotis_col] <= 500
    tby51_ok <- x[tby51_col] >= 10 & x[tby51_col] <= 500
    others_ok <- mean(x[other_cols]) >= 20 & all(x[other_cols] <= 500)
    lotis_ok & tby51_ok & others_ok
  })
  
  cat("Total SNPs:", length(keep), "\n")
  cat("SNPs kept:", sum(keep), "\n")
  cat("SNPs removed:", sum(!keep), "\n\n")
  
  a1_filtered <- a1[keep,]
  a2_filtered <- a2[keep,]
  snps_filtered <- snps[keep,]
  
  fwrite(as.data.table(a1_filtered), 
         file=paste0("ad1_beastfilteredV2_fff_o_lycpool_chrom", chrom, "_noBAT49.txt"), 
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(a2_filtered), 
         file=paste0("ad2_beastfilteredV2_fff_o_lycpool_chrom", chrom, "_noBAT49.txt"), 
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(snps_filtered), 
         file=paste0("beastfilteredsnpsV2_fff_o_lycpool_chrom", chrom, "_noBAT49.txt"), 
         sep="\t", col.names=FALSE)
}
```
check values:
##Total pre-filtered snps: 15,387,914
##Total post-filter snps kept: 632,430
##Total removed: 14,784,484 (95% ish...)

#which files are being counted
```
wc -l beastfilteredsnpsV2_fff_o_lycpool_chrom*.txt
```
# count only the noBAT49 files
```
wc -l beastfilteredsnpsV2_fff_o_lycpool_chrom*_noBAT49.txt | tail -1
```
# non-noBAT49 files separately
```
wc -l beastfilteredsnpsV2_fff_o_lycpool_chrom*.txt | grep -v "noBAT49"
```

## Snp counts and proportion
# SNP counts and proportions
```
total=632430

for f in beastfilteredsnpsV2_fff_o_lycpool_chrom*_noBAT49.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportions_noBAT49.txt

cat snp_counts_proportions_noBAT49.txt

# Mean total coverage (ad1 + ad2) per chromosome per sample
for f in ad1_beastfilteredV2_fff_o_lycpool_chrom*_noBAT49.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    f2=$(echo $f | sed 's/ad1/ad2/')
    paste $f $f2 | awk -v chrom="$chrom" -v ncol=26 '{
        for (i=1; i<=ncol; i++) {
            sum[i] += $i + $(i+ncol)
            count[i]++
        }
    }
    END {
        for (i=1; i<=ncol; i++)
            print chrom, i, sum[i]/count[i]
    }'
done > mean_coverage_total_noBAT49.txt

head mean_coverage_total_noBAT49.txt
```
# Regenerate sample names without BAT49
```
grep "^#CHROM" ../fff_o_lycpool_chrom1.vcf | tr '\t' '\n' | \
    tail -n +10 | \
    grep -v "BAT49" | \
    nl -nrz -w1 -v1 > sample_names_noBAT49.txt

cat sample_names_noBAT49.txt


awk 'NR==FNR{name[$1]=$2; next} {print $1, name[$2], $3}' \
    sample_names_noBAT49.txt mean_coverage_total_noBAT49.txt \
    > mean_coverage_total_noBAT49_named.txt

head mean_coverage_total_noBAT49_named.txt

# Header
echo -e "chrom\t$(grep -v "BAT49" ../fff_o_lycpool_chrom1.vcf | grep "^#CHROM" | \
    tr '\t' '\n' | tail -n +10 | grep -v "BAT49" | tr '\n' '\t' | sed 's/\t$//')" \
    > coverage_matrix_noBAT49.txt

# Matrix
awk '
{
    val[$1][$2] = $3
    pops[$2] = 1
}
END {
    n = asorti(pops, poplist)
    for (c = 1; c <= 23; c++) {
        chrom = "chrom" c
        printf chrom
        for (i = 1; i <= n; i++)
            printf "\t" val[chrom][poplist[i]]
        printf "\n"
    }
}' mean_coverage_total_noBAT49_named.txt >> coverage_matrix_noBAT49.txt

cat coverage_matrix_noBAT49.txt
```

# mkBeast.R- no BAT49


# check what files get picked up with new pattern
list.files(pattern="ad1_beastfilteredV2_fff_o_lycpool_chrom.*noBAT49")

# Check IDs file
```r
read.table("IDs.txt", header=FALSE)
library(data.table)

a1f <- list.files(pattern="ad1_beastfilteredV2_fff_o_lycpool_chrom.*noBAT49")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
asnp <- list.files(pattern="beastfilteredsnpsV2_fff_o_lycpool_chrom.*noBAT49")
asnp <- asnp[1:23]

N <- length(a1f)
ids <- read.table("sample_names_noBAT49.txt", header=FALSE)

# Extract chromosome numbers from filenames
temp <- gsub("ad1_beastfilteredV2_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub("_noBAT49.txt", "", temp)

J_total <- nrow(ids)  # 26 samples

for(i in 1:N){
    cat(i, "\n")
    out <- paste("BEAST_chrom_noBAT49_", chrom[i], ".fasta", sep="")
    
    a1 <- as.matrix(fread(a1f[i], header=F))
    a2 <- as.matrix(fread(a2f[i], header=F))
    n <- a1 + a2
    p <- a2/(a1 + a2)  ## non-ref allele frequency
    p[n < 5] <- NA
    
    J <- dim(p)[2]  # should be 26
    L <- dim(p)[1]  # number of SNPs for this chromosome
    
    snps <- as.data.frame(fread(asnp[i], header=FALSE))
    
    for(j in 1:J){
        nx <- as.numeric(p[,j] > .5) + 1
        ss <- rep("N", L)
        jx <- which(is.na(p[,j]) == FALSE)
        for(l in jx){
            ss[l] <- snps[l, nx[l]]
        }
        SS1 <- paste(ss, collapse="")
        cat(">", ids[j,2], "\n", file=out, append=TRUE, sep="")
        cat(SS1, "\n", file=out, append=TRUE, sep="")
    }
}

```

### For thinning for BEAST: need to thin somewhere between .12 and .35 given chrom size. Trying two options:

## Number 1- fixed number of snps per chrom
```r
library(data.table)
miss <- vector("list", 23)
for(i in 1:23){
    ifile <- paste("text_max_chrom", i, ".fasta", sep="")
    dat <- fread(ifile, header=FALSE)
    miss[[i]] <- apply(dat==5, 2, mean)
}

## retain exactly 5000 SNPs with no missing data
keepSNPs <- vector("list", 23)
target <- 5000

for(i in 1:23){
    xx <- which(miss[[i]] == 0)  ## no missing data
    cat("Chrom", i, "- SNPs with no missing data:", length(xx), "\n")
    n_keep <- min(target, length(xx))  ## in case fewer than 5000 available
    keepSNPs[[i]] <- sort(sample(xx, n_keep, replace=FALSE))
}

for(i in 1:23){
    out <- paste("keepSNPs_max_chrom", i, sep="")
    write.table(keepSNPs[[i]], file=out, row.names=FALSE, col.names=FALSE, quote=FALSE)
}
save(list=ls(), file="snps_maxFixed.rdat")
```
This shit didn't work, go to option 2

In R - how many SNPs have NO missing data per chromosome?


## Number 2- proportion for all chromosomes (4-6k)

```r
## had to find proportion first:
in R:
counts <- c(35697,32239,29561,31106,37354,28824,28523,30358,23719,
            28649,27816,26005,29827,22933,27622,25948,27168,24320,
            29748,22991,21464,17048,23510)

total_target <- 6000
props <- counts / sum(counts)
per_chrom <- floor(props * total_target)
per_chrom
sum(per_chrom)

6000/632430
Prop= .00949

GetSNPSubsetMax.R

library(data.table)
miss<-vector("list",23)
for(i in 1:23){
        ifile<-paste("text_chrom",i,".fasta",sep="")
        dat<-fread(ifile,header=FALSE)
        miss[[i]]<-apply(dat==5,2,mean)

}
## retain 6k snps each  with the appropriate conditions
keepSNPs<-vector("list",23)
prop<-0.00949

for(i in 1:23){
        xx<-which(miss[[i]]==0) ## no missing data
        keepSNPs[[i]]<-sort(sample(xx,floor(length(miss[[i]])*prop),replace=FALSE))
}

for(i in 1:23){
        out<-paste("keepSNPs_maxProp_chrom",i,sep="")
        write.table(keepSNPs[[i]],file=out,row.names=FALSE,col.names=FALSE,quote=FALSE)
}
save(list=ls(),file="snps_maxProp.rdat")

```
SubSetFasta.pl
```pl

foreach $i (1..23){
	open(IN,"keepSNPs_maxProp_chrom$i") or die "failed to open snps file $i\n";
	#open(IN,"keepSNPs_chrom$i") or die "failed to open snps file $i\n";
	$j = 0;
	while(<IN>){
		chomp;
		push (@{$snps[$i]},$_);
	}
	close(IN);
}

open(OUT, ">lyc_genomemax_noBAT49.fasta") or die "failed to write\n";

%seq;
foreach $i (1..23){
	open(IN,"BEAST_chrom_noBAT49_$i.fasta") or die "failed to open snps file $i\n";
	while(<IN>){
		chomp;
		if(m/^>(\S+)/){
			$id = $1;
			if($i == 1){
				@{$seq{$id}} = ();
			}
		} else {
			foreach $snp (@{$snps[$i]}){
				$c = substr $_,$snp-1, 1;
				unless(length($c)==1){
					print "$c\n";
				}
				push(@{$seq{$id}}, $c);
			}
		}
	}
	close(IN);
}

foreach $pop (sort keys %seq){
	$str = join("",@{$seq{$pop}});
	unless($pop =~ m/rep/){
		$pop =~ s/Lyc-//;
		##$pop =~ s/\d+//;
		print OUT ">$pop\n";
		print OUT "$str\n";
	}
}
close(OUT);
```


# Removing BAT49 and TBY51
this can all be found in /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/noBAT49_noTBY51
## Removing both populations
### Remove columns 3 and 24 (BAT49 and TBY51) from all original ad1, ad2, and SNP files
```sh
for f in ../ad1_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    cut -f1-2,4-23,25- $f > ad1_fff_o_lycpool_${chrom}_noBAT49_noTBY51.txt
done

for f in ../ad2_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    cut -f1-2,4-23,25- $f > ad2_fff_o_lycpool_${chrom}_noBAT49_noTBY51.txt
done

for f in ../beastfilteredsnps_fff_o_lycpool_chrom*.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    cut -f1-2,4-23,25- $f > beastfilteredsnps_fff_o_lycpool_${chrom}_noBAT49_noTBY51.txt
done
```
Verify
```sh
#should be 25 columns
awk '{print NF}' ad1_fff_o_lycpool_chrom1_noBAT49_noTBY51.txt | sort -u
#snp file should be 2 columns
awk '{print NF}' beastfilteredsnps_fff_o_lycpool_chrom1_noBAT49_noTBY51.txt | sort -u
# row counts should match origionals
wc -l ad1_fff_o_lycpool_chrom1_noBAT49_noTBY51.txt ../ad1_fff_o_lycpool_chrom1.txt
#check what column Lotis is
grep "^#CHROM" ../fff_o_lycpool_chrom1.vcf | tr '\t' '\n' | \
>     tail -n +10 | \
>     grep -v "BAT49" | \
>     grep -v "TBY51" | \
>     nl | grep "LOTIS"
##column 13
```
Filter in R
```R
library(data.table)

lotis_col <- 13
other_cols <- setdiff(1:25, lotis_col)  # all 24 other pops

for (chrom in 1:23) {
  cat("Processing chrom", chrom, "\n")
  
  a1 <- as.matrix(fread(paste0("ad1_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51.txt"), header=FALSE))
  a2 <- as.matrix(fread(paste0("ad2_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51.txt"), header=FALSE))
  snps <- as.matrix(fread(paste0("beastfilteredsnps_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51.txt"), header=FALSE))
  
  cov <- a1 + a2
  
  keep <- apply(cov, 1, function(x) {
    lotis_ok <- x[lotis_col] >= 10 & x[lotis_col] <= 500
    others_ok <- mean(x[other_cols]) >= 20 & all(x[other_cols] <= 500)
    lotis_ok & others_ok
  })
  
  cat("Total SNPs:", length(keep), "\n")
  cat("SNPs kept:", sum(keep), "\n")
  cat("SNPs removed:", sum(!keep), "\n\n")
  
  a1_filtered <- a1[keep,]
  a2_filtered <- a2[keep,]
  snps_filtered <- snps[keep,]
  
  fwrite(as.data.table(a1_filtered),
         file=paste0("ad1_beastfilteredV2_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51.txt"),
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(a2_filtered),
         file=paste0("ad2_beastfilteredV2_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51.txt"),
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(snps_filtered),
         file=paste0("beastfilteredsnpsV2_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51.txt"),
         sep="\t", col.names=FALSE)
}
```
Check total in bash to make sure it makes sense
```sh
wc -l beastfilteredsnpsV2_fff_o_lycpool_chrom*_noBAT49_noTBY51.txt | tail -1
```
Total: 2562462 total

Get snp counts and proportions, and mean coverage
```sh
# SNP counts and proportions
total=2562462

for f in beastfilteredsnpsV2_fff_o_lycpool_chrom*_noBAT49_noTBY51.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportions_noBAT49_noTBY51.txt

cat snp_counts_proportions_noBAT49_noTBY51.txt
```
Now the mean
```sh
for f in ad1_beastfilteredV2_fff_o_lycpool_chrom*_noBAT49_noTBY51.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    f2=$(echo $f | sed 's/ad1/ad2/')
    paste $f $f2 | awk -v chrom="$chrom" -v ncol=25 '{
        for (i=1; i<=ncol; i++) {
            sum[i] += $i + $(i+ncol)
            count[i]++
        }
    }
    END {
        for (i=1; i<=ncol; i++)
            print chrom, i, sum[i]/count[i]
    }'
done > mean_coverage_total_noBAT49_noTBY51.txt

head mean_coverage_total_noBAT49_noTBY51.txt
```
Add the pop names and build coverage matrix
```sh
# Generate sample names without BAT49 and TBY51
grep "^#CHROM" ../fff_o_lycpool_chrom1.vcf | tr '\t' '\n' | \
    tail -n +10 | \
    grep -v "BAT49" | \
    grep -v "TBY51" | \
    nl -nrz -w1 -v1 > sample_names_noBAT49_noTBY51.txt

cat sample_names_noBAT49_noTBY51.txt

awk 'NR==FNR{name[$1]=$2; next} {print $1, name[$2], $3}' \
    sample_names_noBAT49_noTBY51.txt mean_coverage_total_noBAT49_noTBY51.txt \
    > mean_coverage_total_noBAT49_noTBY51_named.txt

head mean_coverage_total_noBAT49_noTBY51_named.txt

# Header
echo -e "chrom\t$(grep "^#CHROM" ../fff_o_lycpool_chrom1.vcf | tr '\t' '\n' | \
    tail -n +10 | grep -v "BAT49" | grep -v "TBY51" | \
    tr '\n' '\t' | sed 's/\t$//')" > coverage_matrix_noBAT49_noTBY51.txt

# Matrix
awk '
{
    val[$1][$2] = $3
    pops[$2] = 1
}
END {
    n = asorti(pops, poplist)
    for (c = 1; c <= 23; c++) {
        chrom = "chrom" c
        printf chrom
        for (i = 1; i <= n; i++)
            printf "\t" val[chrom][poplist[i]]
        printf "\n"
    }
}' mean_coverage_total_noBAT49_noTBY51_named.txt >> coverage_matrix_noBAT49_noTBY51.txt

cat coverage_matrix_noBAT49_noTBY51.txt

#check that it all looks good
# all rows (chroms) should have exactly 26 fields (chrom + 25 pops)
awk '{print NF, $1}' coverage_matrix_noBAT49_noTBY51.txt | sort -u
```
# Compute Invariants for BEAST one liner
Need to generate the basecounts.txt file using countBases_noTBY51.pl
```pl
#!/usr/bin/perl

open(IN, "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta") or die "failed to read\n";
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
```
Output: baseCounts.txt
Next need snpCounts.txt, generated by countBases_noTBY51.pl
```pl
#!/usr/bin/perl
open(IN, "lyc_genomemax_noBAT49_noTBY51.fasta") or die "failed to read\n";
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
perl countSNPs_noTBY51.pl > snpCounts.txt
```
Check output
```sh
wc -l snpCounts.txt
head -5 snpCounts.txt
```
Use ComputeInvariant.R for the invariant sites across pops
```R
cnts <- read.table("baseCounts.txt", header=FALSE)

## get big scaffolds, top 23
totals <- apply(cnts[,-1], 1, sum)
chr <- which(totals >= 9211676)

## get total A, C, G, T
bcnt <- apply(cnts[chr,-1], 2, sum)

## proportion to match noBAT49_noTBY51 subsetting (6000/2562462)
prop <- 0.002341
sbcnt <- floor(bcnt * prop)

## get SNP bases
snps <- read.table("snpCounts.txt", header=FALSE)

## average across all 25 populations
snpCnts <- floor(apply(snps[,-1], 2, mean))
invar <- sbcnt - snpCnts
invar
```
Values for noBAT.49_noTBY51
	A: 343080
	C: 197056
	G: 196907
	T: 342642

Did the same for just noBAT49
	A: 1,395,334
	C: 803,333
	G: 802,810
	T: 1,393,667

# BEAST and CASTER input files
```sh
# Copy scripts from noBAT49 dir
cp ../noBAT49/filter_noBAT49.R .
cp ../noBAT49/mkNumericFasta_noBAT49.pl .
cp ../noBAT49/SubSetFasta.pl .
cp ../noBAT49/GetSNPSubsetMax.R .
```
Run BEAST input in R
```R
library(data.table)

a1f <- list.files(pattern="ad1_beastfilteredV2_fff_o_lycpool_chrom.*noBAT49_noTBY51")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
asnp <- list.files(pattern="beastfilteredsnpsV2_fff_o_lycpool_chrom.*noBAT49_noTBY51")
asnp <- asnp[1:23]

N <- length(a1f)
ids <- read.table("sample_names_noBAT49_noTBY51.txt", header=FALSE)

temp <- gsub("ad1_beastfilteredV2_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub("_noBAT49_noTBY51.txt", "", temp)

for(i in 1:N){
    cat(i, "\n")
    out <- paste("BEAST_chrom_noBAT49_noTBY51_", chrom[i], ".fasta", sep="")
    
    a1 <- as.matrix(fread(a1f[i], header=F))
    a2 <- as.matrix(fread(a2f[i], header=F))
    n <- a1 + a2
    p <- a2/(a1 + a2)
    p[n < 5] <- NA
    
    J <- dim(p)[2]  # 25 samples
    L <- dim(p)[1]
    
    snps <- as.data.frame(fread(asnp[i], header=FALSE))
    
    for(j in 1:J){
        nx <- as.numeric(p[,j] > .5) + 1
        ss <- rep("N", L)
        jx <- which(is.na(p[,j]) == FALSE)
        for(l in jx){
            ss[l] <- snps[l, nx[l]]
        }
        SS1 <- paste(ss, collapse="")
        cat(">", ids[j,2], "\n", file=out, append=TRUE, sep="")
        cat(SS1, "\n", file=out, append=TRUE, sep="")
    }
}
```
Ran with mkBeastDat_noBAT49_noTBY51.R
output: BEAST_chrom_nobat or tby
Check sequence length and sample size
```sh
# Check 25 samples per file
grep -c ">" BEAST_chrom_noBAT49_noTBY51_1.fasta

# Check sequence length matches SNP count for chrom1 (159,408)
grep -v ">" BEAST_chrom_noBAT49_noTBY51_1.fasta | head -1 | wc -c
```
Make CASTER input files
```sh
library(data.table)

a1f <- list.files(pattern="ad1_beastfilteredV2_fff_o_lycpool_chrom.*noBAT49_noTBY51")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
asnp <- list.files(pattern="beastfilteredsnpsV2_fff_o_lycpool_chrom.*noBAT49_noTBY51")
asnp <- asnp[1:23]

N <- length(a1f)
ids <- read.table("sample_names_noBAT49_noTBY51.txt", header=FALSE)

temp <- gsub("ad1_beastfilteredV2_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub("_noBAT49_noTBY51.txt", "", temp)

for(i in 1:N){
    cat(i, "\n")
    out <- paste("CAST_chrom_noBAT49_noTBY51_", chrom[i], ".fasta", sep="")
    
    a1 <- as.matrix(fread(a1f[i], header=F))
    a2 <- as.matrix(fread(a2f[i], header=F))
    n <- a1 + a2
    p <- a2/(a1 + a2)
    p[is.na(p)] <- 0.001
    
    J <- dim(p)[2]  # 25 samples
    L <- dim(p)[1]
    
    snps <- as.data.frame(fread(asnp[i], header=FALSE))
    
    for(j in 1:J){
        nx <- rbinom(n=L, size=1, prob=p[,j]) + 1
        ss <- rep(NA, L)
        for(l in 1:L){
            ss[l] <- snps[l, nx[l]]
        }
        SS1 <- paste(ss, collapse="")
        cat(">", ids[j,2], "\n", file=out, append=TRUE, sep="")
        cat(SS1, "\n", file=out, append=TRUE, sep="")
    }
}
```
Run with RscriptmkCaster_noBAT49_noTBY51.R
output: CAST_chrom_no bat or tby
Check file content
```sh
# Check 25 samples per file
grep -c ">" CAST_chrom_noBAT49_noTBY51_1.fasta

# Check sequence length matches SNP count for chrom1 (159,408)
grep -v ">" CAST_chrom_noBAT49_noTBY51_1.fasta | head -1 | wc -c
```
SNP subsetting for BEAST using mkNumericFasta.pl
```sh
#!/usr/bin/perl
## this is to figure out which SNPs are variable in the fasta
foreach $i (1..23){
    system "grep -v \"^>\" BEAST_chrom_noBAT49_noTBY51_$i.fasta | perl -pe 'tr/ACGTN/12345/' | sed 's/./& /g' > text_chrom_noTBY51_$i.fasta\n";
}
```
Run with perl mkNumericFasta.pl
 Check:
 ```sh
wc -l text_chrom_noTBY51_1.fasta  # should be 25
head -1 text_chrom_noTBY51_1.fasta | wc -w  # should be 159408
```
Calculate the new proportion for ~6k snps
in R
```R
counts <- c(159408, 139747, 123916, 130047, 151460, 125298, 121332,
            130806, 98258, 115154, 117189, 106853, 112221, 99622,
            104960, 100377, 103636, 90251, 108431, 90845, 77426,
            57636, 97589)

total_target <- 6000
prop <- total_target / sum(counts)
cat("prop =", prop, "\n")

props <- counts / sum(counts)
per_chrom <- floor(props * total_target)
per_chrom
sum(per_chrom)
##sum= 5990
```
SNP subsetting (random sampling to get snps down to about 6K)
In R
```R
library(data.table)

miss <- vector("list", 23)
for(i in 1:23){
    ifile <- paste("text_chrom_noTBY51_", i, ".fasta", sep="")
    dat <- fread(ifile, header=FALSE)
    miss[[i]] <- apply(dat==5, 2, mean)
}

prop <- 6000 / 2562462  # = 0.002341

keepSNPs <- vector("list", 23)
for(i in 1:23){
    xx <- which(miss[[i]] == 0)
    n_keep <- floor(length(miss[[i]]) * prop)
    cat("Chrom", i, "- keeping:", n_keep, "from", length(xx), "complete SNPs\n")
    keepSNPs[[i]] <- sort(sample(xx, n_keep, replace=FALSE))
}

for(i in 1:23){
    out <- paste("keepSNPs_maxProp_noTBY51_chrom", i, sep="")
    write.table(keepSNPs[[i]], file=out, row.names=FALSE, col.names=FALSE, quote=FALSE)
}
save(list=ls(), file="snps_maxProp_noTBY51.rdat")
```
Can run in R or as Rscript GetSNPSubsetMax_noTBY51.R
output- keepSNPs_maxProp_noTBY51_chrom*
Check with: wc -l keepSNPs_maxProp_noTBY51_chrom* or ls -lh

Run SubSetFasta.pl
```pl
#!/usr/bin/perl
#
# this subsets and concatenates a set of SNPs from fasta

foreach $i (1..23){
	open(IN,"keepSNPs_maxProp_noTBY51_chrom$i") or die "failed to open snps file $i\n";
	#open(IN,"keepSNPs_chrom$i") or die "failed to open snps file $i\n";
	$j = 0;
	while(<IN>){
		chomp;
		push (@{$snps[$i]},$_);
	}
	close(IN);
}

open(OUT, ">lyc_genomemax_noBAT49_noTBY51.fasta") or die "failed to write\n";

%seq;
foreach $i (1..23){
	open(IN,"BEAST_chrom_noBAT49_noTBY51_$i.fasta") or die "failed to open snps file $i\n";
	while(<IN>){
		chomp;
		if(m/^>(\S+)/){
			$id = $1;
			if($i == 1){
				@{$seq{$id}} = ();
			}
		} else {
			foreach $snp (@{$snps[$i]}){
				$c = substr $_,$snp-1, 1;
				unless(length($c)==1){
					print "$c\n";
				}
				push(@{$seq{$id}}, $c);
			}
		}
	}
	close(IN);
}

foreach $pop (sort keys %seq){
	$str = join("",@{$seq{$pop}});
	unless($pop =~ m/rep/){
		$pop =~ s/Lyc-//;
		##$pop =~ s/\d+//;
		print OUT ">$pop\n";
		print OUT "$str\n";
	}
}
close(OUT);
```
Verify
```sh
#pops
grep ">" lyc_genomemax_noBAT49_noTBY51.fasta
###snp total
grep -v ">" lyc_genomemax_noBAT49_noTBY51.fasta | head -1 | wc -c
```
Convert to nexus in bash
```sh
seqmagick convert --output-format nexus --alphabet dna \
    lyc_genomemax_noBAT49_noTBY51.fasta \
    lyc_genomemax_noBAT49_noTBY51.nex
#check
head -10 lyc_genomemax_noBAT49_noTBY51.nex
```
# Run BEAST
Make the 
cp ../noBAT49/SubBeast.sh .
Ran using SubBeast.sh, changed the input and output to fit this run
Latest tree: lyc_wgs_ranlc_noBAT49_V2.trees
Calibration correct here, need to fix invariant sites


# Trying to compute the invariants.....again
Done in: /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/noBAT49
Just doing noBAT49 first. Run with countBases.pl to make baseCounts_noBAT49_v2.txt

countBases.pl
```pl
#!/usr/bin/perl

open(IN, "/uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta") or die "failed to read\n";
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
```
perl countBases.pl > baseCounts_noBAT49_V2.txt
Output: baseCounts_noBAT49_V2.txt
Total filtered SNPs -632430
```
wc -l beastfilteredsnpsV2_fff_o_lycpool_chrom*_noBAT49.txt | tail -1
```
Subsetted SNPs- 5991
```
grep -v ">" lyc_genomemax_noBAT49.fasta | head -1 | wc -c
```
5991 (+1 for added row)

##find prop of bases that are acutally invariant (compared between total genome and subset genome for _BAT49)
In R
```r
cnts <- read.table("baseCounts_noBAT49_V2.txt", header=FALSE)
totals <- apply(cnts[,-1], 1, sum)
chr <- which(totals >= 9211676)
bcnt <- apply(cnts[chr,-1], 2, sum)

prop <- 5990 / 632430
cat("Proportion:", prop, "\n")

sbcnt <- floor(bcnt * prop)
cat("Scaled genome counts (A, C, G, T):", sbcnt, "\n")
```
prop= .0095
A: 1393954 C: 803332 G: 802890 T: 1392333

Count snps from subsetted genome using countSNPs_noBAT49.pl
```
#!/usr/bin/perl
open(IN, "lyc_genomemax_noBAT49.fasta") or die "failed to read\n";
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
```
Run with: perl countSNPs_noBAT49.pl > snpCounts_noBAT49.txt
Now calculate invariant sites in R
```r
#read SNP base counts and subtract from scaled genome counts
snps <- read.table("snpCounts_noBAT49.txt", header=FALSE)
dim(snps)  # should be 26 x 5
snpCnts <- floor(apply(snps[,-1], 2, mean))
cat("Mean SNP counts (A, C, G, T):", snpCnts, "\n")

#calculate invariant counts
invar <- sbcnt - snpCnts
cat("Invariant counts (A, C, G, T):", invar, "\n")
```
A: 1392597 C:801756 G: 801234 T:1390934
Pllug these into the data line in: lyc_wgs_max_ranlc_noBAT49.xml
Run beast with SubBeast.sh

Now do the same for tby51
```
wc -l beastfilteredsnpsV2_fff_o_lycpool_chrom*_noBAT49_noTBY51.txt | tail -1
2562462 total
grep -v ">" lyc_genomemax_noBAT49_noTBY51.fasta | head -1 | wc -c
5991
```
run: 
perl countSNPs_noTBY51.pl > snpCounts_noTBY51.txt
```
wc -l snpCounts_noTBY51.txt
head -3 snpCounts_noTBY51.txt
```
In R
```r
# Step 1: read genome base counts
cnts <- read.table("baseCounts.txt", header=FALSE)
dim(cnts)  # should be 1651 x 5

# Step 2: get 23 big chromosomes and sum base counts
totals <- apply(cnts[,-1], 1, sum)
chr <- which(totals >= 9211676)
cat("Number of chromosomes:", length(chr), "\n")
bcnt <- apply(cnts[chr,-1], 2, sum)
cat("Total genome base counts (A, C, G, T):", bcnt, "\n")

# Step 3: scale by proportion
prop <- 5990 / 2562462
cat("Proportion:", prop, "\n")
sbcnt <- floor(bcnt * prop)
cat("Scaled genome counts (A, C, G, T):", sbcnt, "\n")

# Step 4: read SNP counts
snps <- read.table("snpCounts_noTBY51.txt", header=FALSE)
dim(snps)  # should be 25 x 5
snpCnts <- floor(apply(snps[,-1], 2, mean))
cat("Mean SNP counts (A, C, G, T):", snpCnts, "\n")

# Step 5: calculate invariant counts
invar <- sbcnt - snpCnts
cat("Invariant counts (A, C, G, T):", invar, "\n")
```
342578 196768 196619 342140
put into .xml file and run


# An attempt at demographic inference

Need to generate an allele frequency file from my allele depth files. 
Do I want to try this with dadi and then feed that file into moments?
Can I use dadi language in moments?




# Filtering GACT out
need to redo everything after variant calling. Pipelne goes: alignment- remove pcr duplicates - map damage- variant calling- hard filtering 

Hard filtered with HardFilterGACT.sh
```sh
#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=LycLotis_hardfilter_GACT
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load bcftools


cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filtered_GACT

for vcf in fff_o_lycpool_chrom*.vcf; do
    chrom=$(basename "$vcf" .vcf)
    echo "Filtering $chrom..."

    bcftools filter \
        -e '(REF="C" & ALT="T") || (REF="G" & ALT="A")' \
        "$vcf" \
        -o "${chrom}.filtered.vcf"

    before=$(bcftools view -H "$vcf" | wc -l)
    after=$(grep -v "^#" "${chrom}.filtered.vcf" | wc -l)
    echo "$chrom: $before -> $after variants (removed $((before - after)))"
done
```
output: fff_o_lycpool_chrom*.filtered.vcf

Ran AD.sh to get allele depth counts
```sh
#!/usr/bin/bash
#
# extract allele depth AD from biallelic SNPs that passed filtering
#

for f in fff*filtered.vcf
do
        echo "Processing $f"
        out="$(echo $f | sed -e 's/vcf/txt/')"
        echo "Output is ad1_$out"
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\1/g' > ad1_$out
        grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | perl -p -i -e 's/^.+AD\s+//' | perl -p -i -e 's/\S+:(\d+),(\d+)/\2/g' > ad2_$out
done
```
output: ad1 and ad2 _fff_o_lycpool_chrom*.filtered.txt

Ran SNP.sh for snp info
```sh
#!/usr/bin/bash
#
# extract alleles from biallelic SNPs that passed filtering
#

for f in fff*filtered.vcf
do
	echo "Processing $f"
	out="$(echo $f | sed -e 's/vcf/txt/')"
	echo "Output is snps_$out"
	grep ^Sc $f | grep PASS | grep -v [ATCG],[ATCG] | cut -f 4,5 > beastfilteredsnps_$out
done
```
output: beastfilteredsnps_fff_o_lycpool_chrom9.filtered.txt

## skipping pop gen structure for now, will come back to FST and PCA's

# Time calibrated phylogenetic tree
###going to first check the coverage across chromosomes
```sh
# Check number of SNPs per chromosome
wc -l snpinfo_fff_o_lycpool_chrom*.txt

# Check number of columns (should be 27 for 27 samples)
awk '{print NF; exit}' ad1_fff_o_lycpool_chrom1.txt

# Peek at first few lines of both ref and alt
head -3 ad1_fff_o_lycpool_chrom*.filtered.txt
head -3 ad2_fff_o_lycpool_chrom*.filtered.txt

# Extract chrom, pos, ref, alt for PASS biallelic SNPs
for vcf in fff_o_lycpool_chrom*.filtered.vcf; do
    chrom=$(basename "$vcf" .filtered.vcf)
    grep "^Sc" "$vcf" | grep PASS | grep -v "[ATCG],[ATCG]" | \
        awk '{print $1"\t"$2"\t"$4"\t"$5}' > snpinfo_${chrom}.txt
done

# Total SNPs across all chromosomes
wc -l snpinfo_fff_o_lycpool_chrom*.txt
##total: 

# Check ad1 and ad2 line counts match snpinfo
wc -l ad1_fff_o_lycpool_chrom11.txt
wc -l snpinfo_fff_o_lycpool_chrom11.txt
Total: 

# Check columns in ad files
awk '{print NF; exit}' ad1_fff_o_lycpool_chrom11.txt
```
# Checking coverage 
```R
library(data.table)
# Get list of all ad1 files and combine all ad1 files
ad1_files <- list.files(pattern="ad1_fff_o_lycpool_chrom.*\\.filtered\\.txt")
a1 <- as.matrix(rbindlist(lapply(ad1_files, fread, header=FALSE)))
# Do the same for ad2
ad2_files <- list.files(pattern="ad2_fff_o_lycpool_chrom.*\\.filtered\\.txt")
a2 <- as.matrix(rbindlist(lapply(ad2_files, fread, header=FALSE)))
# Check dimensions - should be ~10.9 million rows x 27 columns
dim(a1)
dim(a2)
# check actual coverage values
# Mean coverage per sample
round(colMeans(cov), 1)
# Min and max coverage per sample
round(apply(cov, 2, min), 1)
round(apply(cov, 2, max), 1)
# How many sites have 0 coverage per sample
colSums(cov == 0)
# Distribution of coverage
summary(cov)
```
with this filtering, we still see BAT49 as the worst pop, and Lotis coming in ahead of TBY51 by about the same proporions we saw without GACT filteirng. Going to do a first run without BAT49
```R
library(data.table)

bat49_col <- 3
lotis_col <- 14
tby51_col <- 24
other_cols <- setdiff(1:27, c(bat49_col, lotis_col, tby51_col))

for (chrom in 1:23) {
  cat("Processing chrom", chrom, "\n")
  
  a1 <- as.matrix(fread(paste0("ad1_fff_o_lycpool_chrom", chrom, ".filtered.txt"), header=FALSE))
  a2 <- as.matrix(fread(paste0("ad2_fff_o_lycpool_chrom", chrom, ".filtered.txt"), header=FALSE))
  snps <- as.matrix(fread(paste0("snpinfo_fff_o_lycpool_chrom", chrom, ".txt"), header=FALSE))
  
  # Remove BAT49 column
  a1 <- a1[, -bat49_col]
  a2 <- a2[, -bat49_col]
  
  # Recalculate column indices after removing BAT49
  # LOTIS and TBY51 shift down by 1 since BAT49 (col3) is removed
  lotis_col_new <- lotis_col - 1  # now col 13
  tby51_col_new <- tby51_col - 1  # now col 23
  other_cols_new <- setdiff(1:26, c(lotis_col_new, tby51_col_new))
  
  cov <- a1 + a2
  
  keep <- apply(cov, 1, function(x) {
    lotis_ok <- x[lotis_col_new] >= 10 & x[lotis_col_new] <= 500
    tby51_ok <- x[tby51_col_new] >= 10 & x[tby51_col_new] <= 500
    others_ok <- all(x[other_cols_new] >= 20) & all(x[other_cols_new] <= 500)
    lotis_ok & tby51_ok & others_ok
  })
  
  cat("Total SNPs:", length(keep), "\n")
  cat("SNPs kept:", sum(keep), "\n")
  cat("SNPs removed:", sum(!keep), "\n\n")
  
  a1_filtered <- a1[keep,]
  a2_filtered <- a2[keep,]
  snps_filtered <- snps[keep,]
  
  fwrite(as.data.table(a1_filtered),
         file=paste0("ad1_fff_o_lycpool_chrom", chrom, "_noBAT49.filtered.txt"),
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(a2_filtered),
         file=paste0("ad2_fff_o_lycpool_chrom", chrom, "_noBAT49.filtered.txt"),
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(snps_filtered),
         file=paste0("snpinfo_fff_o_lycpool_chrom", chrom, "_noBAT49.txt"),
         sep="\t", col.names=FALSE)
}
```
To check snps retained per chrom and that BAT49 actually dropped
```r
for(chrom in 1:23){
    f <- paste0("snpinfo_fff_o_lycpool_chrom", chrom, "_noBAT49.txt")
    cat("chrom", chrom, ":", nrow(fread(f, header=FALSE)), "SNPs\n")
}
ncol(fread("ad1_fff_o_lycpool_chrom1_noBAT49.filtered.txt", header=FALSE))

Total snps
sum(c(19486,18442,16663,17919,21176,16947,16293,16824,13054,16072,15734,15014,17331,13143,16067,15001,16144,13913,17358,13060,12481,9854,12627))
```
total: 360603

# SNP counts and proportions 
(in bash)
```sh
total=360603
for f in snpinfo_fff_o_lycpool_chrom*_noBAT49.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportions_noBAT49.txt
cat snp_counts_proportions_noBAT49.txt

# Mean total coverage per chromosome per sample
for f in ad1_fff_o_lycpool_chrom*_noBAT49.filtered.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    f2=$(echo $f | sed 's/ad1/ad2/')
    paste $f $f2 | awk -v chrom="$chrom" -v ncol=26 '{
        for (i=1; i<=ncol; i++) {
            sum[i] += $i + $(i+ncol)
            count[i]++
        }
    }
    END {
        for (i=1; i<=ncol; i++)
            print chrom, i, sum[i]/count[i]
    }'
done > mean_coverage_total_noBAT49.txt
head mean_coverage_total_noBAT49.txt
```

Making BEAST and CAster input files with noBAT49
```sh
library(data.table)

a1f <- list.files(pattern="ad1_fff_o_lycpool_chrom.*_noBAT49\\.filtered\\.txt")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
asnp <- list.files(pattern="snpinfo_fff_o_lycpool_chrom.*_noBAT49\\.txt")
asnp <- asnp[1:23]

N <- length(a1f)
ids <- read.table("sample_names_noBAT49.txt", header=FALSE)

temp <- gsub("ad1_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub("_noBAT49\\.filtered\\.txt", "", temp)

for(i in 1:N){
    cat(i, "\n")
    out <- paste0("BEAST_chrom_noBAT49_", chrom[i], ".fasta")
    
    a1 <- as.matrix(fread(a1f[i], header=F))
    a2 <- as.matrix(fread(a2f[i], header=F))
    n <- a1 + a2
    p <- a2/(a1 + a2)
    p[n < 5] <- NA
    
    J <- dim(p)[2]  # 26 samples
    L <- dim(p)[1]
    
    snps <- as.data.frame(fread(asnp[i], header=FALSE))
    
    for(j in 1:J){
        nx <- as.numeric(p[,j] > .5) + 1
        ss <- rep("N", L)
        jx <- which(is.na(p[,j]) == FALSE)
        for(l in jx){
            ss[l] <- snps[l, nx[l]]
        }
        SS1 <- paste(ss, collapse="")
        cat(">", ids[j,2], "\n", file=out, append=TRUE, sep="")
        cat(SS1, "\n", file=out, append=TRUE, sep="")
    }
}
```
output: BEAST_chrom_noBAT49_#_.fasta
And for CASTER:
```sh
library(data.table)

a1f <- list.files(pattern="ad1_fff_o_lycpool_chrom.*_noBAT49\\.filtered\\.txt")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
asnp <- list.files(pattern="snpinfo_fff_o_lycpool_chrom.*_noBAT49\\.txt")
asnp <- asnp[1:23]

N <- length(a1f)
ids <- read.table("sample_names_noBAT49.txt", header=FALSE)

temp <- gsub("ad1_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub("_noBAT49\\.filtered\\.txt", "", temp)

for(i in 1:N){
    cat(i, "\n")
    out <- paste0("CAST_chrom_noBAT49_", chrom[i], ".fasta")
    
    a1 <- as.matrix(fread(a1f[i], header=F))
    a2 <- as.matrix(fread(a2f[i], header=F))
    n <- a1 + a2
    p <- a2/(a1 + a2)
    p[is.na(p)] <- 0.001
    
    J <- dim(p)[2]  # 26 samples
    L <- dim(p)[1]
    
    snps <- as.data.frame(fread(asnp[i], header=FALSE))
    
    for(j in 1:J){
        nx <- rbinom(n=L, size=1, prob=p[,j]) + 1
        ss <- rep(NA, L)
        for(l in 1:L){
            ss[l] <- snps[l, nx[l] + 2]  # +2 to get ref/alt columns
        }
        SS1 <- paste(ss, collapse="")
        cat(">", ids[j,2], "\n", file=out, append=TRUE, sep="")
        cat(SS1, "\n", file=out, append=TRUE, sep="")
    }
}
```
output: CAST_chrom_noBAT49_*.fasta
Check
```sh
for chrom in $(seq 1 23); do
    seq_len=$(grep -v ">" CAST_chrom_noBAT49_${chrom}.fasta | head -1 | tr -cd 'ATCGN' | wc -c)
    snp_count=$(wc -l < snpinfo_fff_o_lycpool_chrom${chrom}_noBAT49.txt)
    if [ "$seq_len" -eq "$snp_count" ]; then
        echo "chrom${chrom}: PASS (${snp_count} SNPs)"
    else
        echo "chrom${chrom}: FAIL - sequence length ${seq_len} != snp count ${snp_count}"
    fi
done
```
Now subset using mkNumericFasta.pl
```perl
#!/usr/bin/perl
## this is to figure out which SNPs are variable in the fasta

foreach $i (1..23){
    system "grep -v \"^>\" BEAST_chrom_noBAT49_V2_$i.fasta | perl -pe 'tr/ACGTN/12345/' | sed 's/./& /g' > text_chrom_noBAT49_$i.fasta\n";
}
```
Calculate new proportion for about 6K snps in R
Check proporions in snp_counts_proportions_noBAT49.txt to get proportion counts
```less snp_counts_proportions.....
```R

total snps_ 5985
```
# Generate invariant counts for beast
run  countbases.pl to get baseCount
perl countBases.pl > baseCounts_filteredV2_noBAT49.txt Output:baseCounts_filteredV2_noBAT49.txt
perl countSNPs_noBAT49.pl > snpCounts_noBAT49_filteredV2.txt
Rscript ComputeInvariant.R

     V2      V3      V4      V5 
2133273 1229811 1229100 2130801 

# Filtered GACT Removing TBY51 (after BAT49 was already done)
## Generating ad1, ad2 and snp files
```r
library(data.table)

lotis_col <- 13
tby51_col <- 23
other_cols <- setdiff(1:25, c(lotis_col))

for (chrom in 1:23) {
  cat("Processing chrom", chrom, "\n")
  
a1 <- as.matrix(fread(paste0("ad1_fff_o_lycpool_chrom", chrom, "_noBAT49_filteredV2.txt"), header=FALSE))
a2 <- as.matrix(fread(paste0("ad2_fff_o_lycpool_chrom", chrom, "_noBAT49_filteredV2.txt"), header=FALSE))
snps <- as.matrix(fread(paste0("snpinfo_fff_o_lycpool_chrom", chrom, "_noBAT49_filteredV2.txt"), header=FALSE))

  # Remove TBY51 first
  a1 <- a1[, -tby51_col]
  a2 <- a2[, -tby51_col]
  cov <- a1 + a2
  
keep <- apply(cov, 1, function(x) {
    lotis_ok <- x[lotis_col] >= 10 & x[lotis_col] <= 500
    others_ok <- mean(x[other_cols]) >= 20 & all(x[other_cols] <= 500)
    lotis_ok & others_ok  # removed tby51_ok
  })
  
  cat("Total SNPs:", length(keep), "\n")
  cat("SNPs kept:", sum(keep), "\n")
  cat("SNPs removed:", sum(!keep), "\n\n")
  
  a1_filtered <- a1[keep,]
  a2_filtered <- a2[keep,]
  snps_filtered <- snps[keep,]
  
  fwrite(as.data.table(a1_filtered),
         file=paste0("ad1_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51_filtered.txt"),
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(a2_filtered),
         file=paste0("ad2_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51_filtered.txt"),
         sep="\t", col.names=FALSE)
  fwrite(as.data.table(snps_filtered),
         file=paste0("snpinfo_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51_filtered.txt"),
         sep="\t", col.names=FALSE)
}
```
output: ad1, ad2, snp info ad1_fff_o_lycpool_chrom10_noBAT49_noTBY51_filteredV2.txtq

check it
```wc -l snpinfo_fff_o_lycpool_chrom*_noBAT49_noTBY51_filtered.txt | tail -1
  402986 total
```
and snp counts
```sh
for(chrom in 1:23){
    f <- paste0("snpinfo_fff_o_lycpool_chrom", chrom, "_noBAT49_noTBY51_filtered.txt")
    cat("chrom", chrom, ":", nrow(fread(f, header=FALSE)), "SNPs\n")
}
ncol(fread("ad1_fff_o_lycpool_chrom1_noBAT49_noTBY51_filtered.txt", header=FALSE))
```


```sh
# SNP counts and proportions
total=402986
for f in snpinfo_fff_o_lycpool_chrom*_noBAT49_noTBY51_filtered.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    count=$(wc -l < $f)
    echo -e "$chrom\t$count\t$(echo "scale=6; $count/$total" | bc)"
done | sort -t'm' -k2 -n > snp_counts_proportions_noBAT49_noTBY51.txt
cat snp_counts_proportions_noBAT49_noTBY51.txt
```
mean
```
for f in ad1_fff_o_lycpool_chrom*_noBAT49_noTBY51_filtered.txt; do
    chrom=$(echo $f | grep -oP 'chrom\d+')
    f2=$(echo $f | sed 's/ad1/ad2/')
    paste $f $f2 | awk -v chrom="$chrom" -v ncol=25 '{
        for (i=1; i<=ncol; i++) {
            sum[i] += $i + $(i+ncol)
            count[i]++
        }
    }
    END {
        for (i=1; i<=ncol; i++)
            print chrom, i, sum[i]/count[i]
    }'

# add sample names

done > mean_coverage_total_noBAT49_noTBY51.txt
head mean_coverage_total_noBAT49_noTBY51.txt
```
coverage matrix
```
# Generate sample names without BAT49 and TBY51
grep "^#CHROM" fff_o_lycpool_chrom1.filtered.vcf | tr '\t' '\n' | \
    tail -n +10 | \
    grep -v "BAT49" | \
    grep -v "TBY51" | \
    nl -nrz -w1 -v1 > sample_names_noBAT49_noTBY51.txt
cat sample_names_noBAT49_noTBY51.txt

awk 'NR==FNR{name[$1]=$2; next} {print $1, name[$2], $3}' \
    sample_names_noBAT49_noTBY51.txt mean_coverage_total_noBAT49_noTBY51.txt \
    > mean_coverage_total_noBAT49_noTBY51_named.txt
head mean_coverage_total_noBAT49_noTBY51_named.txt

# Header
echo -e "chrom\t$(grep "^#CHROM" fff_o_lycpool_chrom1.filtered.vcf | tr '\t' '\n' | \
    tail -n +10 | grep -v "BAT49" | grep -v "TBY51" | \
    tr '\n' '\t' | sed 's/\t$//')" > coverage_matrix_noBAT49_noTBY51.txt

# Matrix
awk '
{
    val[$1][$2] = $3
    pops[$2] = 1
}
END {
    n = asorti(pops, poplist)
    for (c = 1; c <= 23; c++) {
        chrom = "chrom" c
        printf chrom
        for (i = 1; i <= n; i++)
            printf "\t" val[chrom][poplist[i]]
        printf "\n"
    }
}'
mean_coverage_total_noBAT49_noTBY51_named.txt >> coverage_matrix_noBAT49_noTBY51.txt
cat coverage_matrix_noBAT49_noTBY51.txt

# Check all rows have exactly 26 fields (chrom + 25 pops)
awk '{print NF, $1}' coverage_matrix_noBAT49_noTBY51.txt | sort -u
```
# Run BEAST
```r
library(data.table)

a1f <- list.files(pattern="ad1_fff_o_lycpool_chrom.*_noBAT49_noTBY51_filtered\\.txt")
a1f <- a1f[1:23]
a2f <- gsub("ad1", "ad2", a1f)
asnp <- list.files(pattern="snpinfo_fff_o_lycpool_chrom.*_noBAT49_noTBY51_filtered\\.txt")
asnp <- asnp[1:23]

N <- length(a1f)
ids <- read.table("sample_names_noBAT49_noTBY51.txt", header=FALSE)

temp <- gsub("ad1_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub("_noBAT49_noTBY51_filtered\\.txt", "", temp)

for(i in 1:N){
    cat(i, "\n")
    out <- paste0("BEAST_chrom_noBAT49_noTBY51_", chrom[i], ".fasta")
    
    a1 <- as.matrix(fread(a1f[i], header=F))
    a2 <- as.matrix(fread(a2f[i], header=F))
    n <- a1 + a2
    p <- a2/(a1 + a2)
    p[n < 5] <- NA
    
    J <- dim(p)[2]  # 25 samples
    L <- dim(p)[1]
    
    snps <- as.data.frame(fread(asnp[i], header=FALSE))
    
    for(j in 1:J){
        nx <- as.numeric(p[,j] > .5) + 1
        ss <- rep("N", L)
        jx <- which(is.na(p[,j]) == FALSE)
        for(l in jx){
            ss[l] <- snps[l, nx[l] + 2]
        }
        SS1 <- paste(ss, collapse="")
        cat(">", ids[j,2], "\n", file=out, append=TRUE, sep="")
        cat(SS1, "\n", file=out, append=TRUE, sep="")
    }
}
```
check it
```
for chrom in $(seq 1 23); do
    seq_len=$(grep -v ">" BEAST_chrom_noBAT49_noTBY51_${chrom}.fasta | head -1 | tr -cd 'ATCGN' | wc -c)
    snp_count=$(wc -l < snpinfo_fff_o_lycpool_chrom${chrom}_noBAT49_noTBY51_filtered.txt)
    if [ "$seq_len" -eq "$snp_count" ]; then
        echo "chrom${chrom}: PASS (${snp_count} SNPs)"
    else
        echo "chrom${chrom}: FAIL - sequence length ${seq_len} != snp count ${snp_count}"
    fi
done
```
Numeric Fastas
```perl
#!/usr/bin/perl
## convert BEAST fasta to numeric for SNP subsetting

foreach $i (1..23){
    system "grep -v \"^>\" BEAST_chrom_noBAT49_noTBY51_$i.fasta | perl -pe 'tr/ACGTN/12345/' | sed 's/./& /g' > text_chrom_noBAT49_noTBY51_$i.fasta\n";
}
perl mkNumericFasta_noTBY51.pl
```
output: text_chrom_noBAT49_noTBY51.fasta

Subsetting snps
```r
library(data.table)

miss <- vector("list", 23)
for(i in 1:23){
    ifile <- paste0("text_chrom_noBAT49_noTBY51_", i, ".fasta")
    dat <- fread(ifile, header=FALSE)
    miss[[i]] <- apply(dat==5, 2, mean)
}

# Total SNPs from noBAT49_noTBY51
total_snps <- sum(c(22277, 20489, 18858, 19902, 23617, 18437, 18266, 19136,
                    15003, 17950, 17563, 16639, 19164, 14663, 17616, 16497,
                    17742, 15549, 19042, 14712, 13759, 10907, 15198))
prop <- 6000 / total_snps
cat("prop =", prop, "\n")
cat("total_snps =", total_snps, "\n")

keepSNPs <- vector("list", 23)
for(i in 1:23){
    xx <- which(miss[[i]] == 0)
    n_keep <- floor(length(miss[[i]]) * prop)
    cat("Chrom", i, "- keeping:", n_keep, "from", length(xx), "complete SNPs\n")
    keepSNPs[[i]] <- sort(sample(xx, n_keep, replace=FALSE))
}

for(i in 1:23){
    out <- paste0("keepSNPs_maxProp_noBAT49_noTBY51_chrom", i)
    write.table(keepSNPs[[i]], file=out, row.names=FALSE, col.names=FALSE, quote=FALSE)
}

save(list=ls(), file="snps_maxProp_noBAT49_noTBY51.rdat")
```
And check it
```
total_kept <- sum(sapply(keepSNPs, length))
cat("Total SNPs kept:", total_kept, "\n")

# Confirm no missing data
for(i in 1:23){
    ifile <- paste0("text_chrom_noBAT49_noTBY51_", i, ".fasta")
    dat <- as.matrix(fread(ifile, header=FALSE))
    kept <- keepSNPs[[i]]
    n_missing <- sum(dat[, kept] == 5)
    cat("Chrom", i, "- SNPs kept:", length(kept), "- missing data:", n_missing, "\n")
}
```
SubSetFasta.pl
```perl
#!/usr/bin/perl
foreach $i (1..23){
        open(IN,"keepSNPs_maxProp_noBAT49_noTBY51_chrom$i") or die "failed to open snps file $i\n";
        while(<IN>){
                chomp;
                push (@{$snps[$i]},$_);
        }
        close(IN);
}
open(OUT, ">lyc_genomemax_noBAT49_noTBY51.fasta") or die "failed to write\n";
%seq;
foreach $i (1..23){
        open(IN,"BEAST_chrom_noBAT49_noTBY51_$i.fasta") or die "failed to open snps file $i\n";
        while(<IN>){
                chomp;
                if(m/^>(\S+)/){
                        $id = $1;
                        if($i == 1){
                                @{$seq{$id}} = ();
                        }
                } else {
                        foreach $snp (@{$snps[$i]}){
                                $c = substr $_,$snp-1, 1;
                                unless(length($c)==1){
                                        print "$c\n";
                                }
                                push(@{$seq{$id}}, $c);
                        }
                }
        }
        close(IN);
}

foreach $pop (sort keys %seq){
        $str = join("",@{$seq{$pop}});
        unless($pop =~ m/rep/){
                $pop =~ s/Lyc-//;
                print OUT ">$pop\n";
                print OUT "$str\n";
        }
}
close(OUT);
```
file name: SubsetFasta_noTBY51.pl
output: lyc_genomemax_noBAT49_noTBY51.fasta






















compute invariant sites
```r

cnts <- read.table("baseCounts_filteredV2_noBAT49.txt", header=FALSE)
dim(cnts)
#sum base counts from big 23 chroms
totals <- apply(cnts[,-1], 1, sum)
chr <- which(totals >= 9211676)
cat("Number of chromosomes:", length(chr), "\n")
bcnt <- apply(cnts[chr,-1], 2, sum)
cat("Total genome base counts (A, C, G, T):", bcnt, "\n")

##scale by prop
prop <- 5989 / 402986  # noBAT49_noTBY51 total SNPs
cat("Proportion:", prop, "\n")
sbcnt <- floor(bcnt * prop)
cat("Scaled genome counts (A, C, G, T):", sbcnt, "\n")

##read snp coutns
snps <- read.table("snpCounts_noBAT49_noTBY51.txt", header=FALSE)
dim(snps)  # should be 25 x 5
snpCnts <- floor(apply(snps[,-1], 2, mean))
cat("Mean SNP counts (A, C, G, T):", snpCnts, "\n")

# Step 5: calculate invariant counts
invar <- sbcnt - snpCnts
cat("Invariant counts (A, C, G, T):", invar, "\n")
```
2185096 1259643 1258967 2182583


seqmagick convert --output-format nexus --alphabet dna \
    lyc_genomemax_noBAT49_noTBY51.fasta \
    lyc_genomemax_noBAT49_noTBY51.nex



CASTER trees
```r
library(ape)
pdf ("CasterTrees_noBAT49_noTBY51.pdf", width=9, height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<- paste("cout_noBAT49_noTBY51_", i, sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE)
	title(main=paste("Chrom.",i),cex.main=1.3)
 }
 pdf ("CasterTreesMax_noBAT49_noTBY51.pdf", width=9, height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<- paste("cout_noBAT49_noTBY51_", i, sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE,type="cladogram")
	title(main=paste("Chromosome",i),cex.main=1.3)
}
dev.off()
```
And
```
library(ape)
pdf ("CasterTrees_noBAT49_Filtered_final.pdf", width=9, height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<- paste("cout_noBAT49_filtered_Final", i, sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE)
	title(main=paste("Chrom.",i),cex.main=1.3)
 }
pdf ("CasterTreesMax_noBAT49_filtered_final.pdf", width=9, height=9)
par(mfrow=c(3,3))
par(mar=c(1,1,3,1))
for(i in 1:23){
	inf<-paste("cout_noBAT49_filtered_Final", i, sep="")
	tree<-read.tree(inf)
	plot.phylo(tree,cex=.7,use.edge.length=FALSE,type="cladogram")
	title(main=paste("Chromosome",i),cex.main=1.3)
}
dev.off()
```
outputs: CasterTreesMax_noBAT49_noTBY51.pdf  CasterTrees_noBAT49_noTBY51.pdf
to run caster.pl in bash:
perl -e '
foreach $i (1..23){
    system "/uufs/chpc.utah.edu/common/home/u6047808/bin/ASTER-Linux/bin/caster-site -i CAST_chrom_noBAT49_filteredV2_$i.fasta -o cout_noBAT49_filtered_$i --root Lyc-MEN12 --thread 24\n";
}
'
and 
perl -e '
foreach $i (1..23){
    system "/uufs/chpc.utah.edu/common/home/u6047808/bin/ASTER-Linux/bin/caster-site -i CAST_chrom_noBAT49_noTBY51_$i.fasta -o cout_noBAT49_noTBY51_$i --root Lyc-MEN12 --thread 24\n";
}
'




# Using Treemix
make the input files for treemix, using mkTreeMixin.R

```R
## make input for treemix
library(data.table)

a1f <- list.files(pattern="ad1_fff_o_lycpool_chrom[0-9]+_noBAT49_noTBY51_filtered")
a2f <- gsub("ad1", "ad2", a1f)
asnp <- gsub("ad1", "snpinfo", a1f)

N <- length(a1f)

ids <- read.table("sample_names_clean.txt", header=FALSE)
idx <- sub(pattern="^[0-9]+\t", replacement="", x=ids[,1])
idx <- sub(pattern="Lyc-", replacement="", x=idx)

temp <- gsub("ad1_fff_o_lycpool_chrom", "", a1f)
chrom <- gsub("_noBAT49_noTBY51_filtered.txt", "", temp)

reps <- grep(pattern="rep", x=idx)

for(i in 1:N){
    cat(i, "\n")
    a1 <- as.matrix(fread(a1f[i], header=FALSE))
    a2 <- as.matrix(fread(a2f[i], header=FALSE))
    L <- dim(a1)[1]
    J <- dim(a1)[2]
    combPa <- paste(a1, a2, sep=",")
    PaMat <- matrix(combPa, nrow=L, ncol=J, byrow=FALSE)
    colnames(PaMat) <- idx
    if(length(reps) > 0){
        write.table(file=paste("treemix_in_ch", chrom[i], ".txt", sep=""),
                    PaMat[,-reps], quote=FALSE, row.names=FALSE)
    } else {
        write.table(file=paste("treemix_in_ch", chrom[i], ".txt", sep=""),
                    PaMat, quote=FALSE, row.names=FALSE)
    }
}
```
output: treemix_in_ch*.txt
Ran treemix for 0 to xx migration (admixture) events on the noTBYandBAT dataset for all 23 chroms

Zip the files before running: gzip treemix_in_ch*.txt
Used: run_max_treemix.pl
```pl
#!/usr/bin/perl
#
## version of treemix fork script that obtains the ML result across $N runs


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);

$N = 20;

foreach $fi (@ARGV){
	$fi =~ m/(\d+)/ or die "failed here $fi\n";
	$ch = $1;
	foreach $m (0..8){
		$pm->start and next;
		$out = "tro_ch$ch"."_m$m";
		system "treemix -i $fi -o $out -k 100 -m $m -root MEN12\n";
		open(IN, "$out\.llik");
		$a = <IN>;
		$a = <IN>;
		chomp($a);
		$a =~ m/:\s+(\-[0-9\.]+)/ or die "can't find the ll: $a\n";
		$ll = $1;
		close(IN);
		print "starting ll = $ll\n";
		foreach $i (1..$N){
			$tout = "temp_tro_ch$ch"."_m$m";
			system "treemix -i $fi -o $tout -k 100 -m $m -root MEN12\n";
			open(IN, "$tout\.llik");
			$a = <IN>;
			$a = <IN>;
			chomp($a);
			$a =~ m/:\s+(\-[0-9\.]+)/ or die "can't find the ll: $a\n";
			$llt = $1;
			if($llt > $ll){
				$ll = $llt;
				system "mv $tout\.treeout.gz $out\.treeout.gz\n";
				system "mv $tout\.modelcov.gz $out\.modelcov.gz\n";
				system "mv $tout\.vertices.gz $out\.vertices.gz\n";
				system "mv $tout\.cov.gz $out\.cov.gz\n";
				system "mv $tout\.covse.gz $out\.covse.gz\n";
				system "mv $tout\.edges.gz $out\.edges.gz\n";
				system "mv $tout\.llik $out\.llik\n";
			}
		}
		$pm->finish;
	}
}

$pm->wait_all_children;
```

Make 
Submitted with Sub_Treemix.sh

```sh
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=tmix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

cd /uufs/chpc.utah.edu/common/home/gompert-group6/data/LycLotis/ZLycLotis/filtered_GACT

perl run_max_treemix.pl treemix_in_ch*.txt.gz
```



