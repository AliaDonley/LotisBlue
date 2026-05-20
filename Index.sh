#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=Indexingbai
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu


module load samtools
## version 1.16


cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

for bam in filt_dedup_*.bam
do
echo "Indexing $bam..."
samtools index "$bam"
done
