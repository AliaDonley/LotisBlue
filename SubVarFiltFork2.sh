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
