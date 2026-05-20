#!/bin/sh
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=Softclip
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu


module load perl
cd /uufs/chpc.utah.edu/common/home/u6047808/Xerces/ADXerces

perl CountSoftClipsA.pl algn_xer-*.sam
