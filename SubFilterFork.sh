#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=filterfork
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

module load samtools

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

perl /scratch/general/nfs1/u6000989/LycLotisBams/FilterFork.pl dedup_*.bam
