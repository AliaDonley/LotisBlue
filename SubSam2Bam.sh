#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=Bam2Sam
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load perl
module load samtools
module load bwa

cd /uufs/chpc.utah.edu/common/home/u6047808/Xerces/ADXerces

perl Sam2Bam.pl filename*.bam
