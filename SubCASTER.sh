#!/bin/bash
#SBATCH --time=80:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=Caster
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filtered_GACT

perl RunCasterSite.pl
