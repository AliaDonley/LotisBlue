#!/bin/sh
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=lotismapdamagefiltered
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load perl
module load samtools
module load bwa

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

perl MapDamFork.pl *bam
