#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=beast
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alia.donley@usu.edu

cd /uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/noBAT49
ml beast
max=7
count=0

for i in ch1 ch2 ch3 ch4 ch5 ch6;
do
    beast -prefix $i lyc_wgs_max_ranlc_noBAT49.xml &
    ((count++))

    if ((count >= max)); then
        wait -n
        ((count--))
    fi
done

wait
