#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=beast
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis

ml perl

perl mkNumericFasta.pl
