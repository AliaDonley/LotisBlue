#!/bin/bash
#SBATCH --job-name=casterwin
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=casterwin_%j.out
#SBATCH --error=casterwin_%j.err
# Usage: sbatch run_quartet.slurm A B C OUTGROUP
set -euo pipefail

CASTER_DIR=/uufs/chpc.utah.edu/common/home/u6047808/LycLotis/ZLycLotis/filtered_GACT/Caster
WINBIN=/uufs/chpc.utah.edu/common/home/u6047808/bin/wins/MASTERWORK/bin/slidingwindow
ALN_PREFIX=gact_max_chrom

P1=$1; P2=$2; P3=$3; Po=$4
tag() { echo "${1#Lyc-}"; }
Q="$(tag $P1)x$(tag $P2)x$(tag $P3)"

cd "$CASTER_DIR"
for s in WinSubAlign.pl RunWindows.pl; do
  [ -f "$s" ] || { echo "$s not found" >&2; exit 1; }
done
mkdir -p "win${Q}" && cd "win${Q}"

printf "%s\n%s\n%s\n%s\n" "$P1" "$P2" "$P3" "$Po" > "Sub${Q}.txt"
printf "%s\tA\n%s\tB\n%s\tC\n%s\tOutgroup\n" "$P1" "$P2" "$P3" "$Po" > "Mapping${Q}.txt"
for i in $(seq 1 23); do ln -sf "../${ALN_PREFIX}${i}.fasta" .; done

nh=$(grep -c '^>' "${ALN_PREFIX}1.fasta"); nl=$(grep -vc '^>' "${ALN_PREFIX}1.fasta")
[ "$nh" -eq "$nl" ] || { echo "FASTA is multi-line" >&2; exit 1; }

perl ../WinSubAlign.pl "Sub${Q}.txt" $(for i in $(seq 1 23); do echo ${ALN_PREFIX}${i}.fasta; done)
n=$(grep -c '^>' "Sub${Q}_${ALN_PREFIX}1.fasta")
[ "$n" -eq 4 ] || { echo "Only $n of 4 taxa found; check names" >&2; exit 1; }

perl ../RunWindows.pl "$Q" "$WINBIN" "$ALN_PREFIX"
echo "Done: ${CASTER_DIR}/win${Q}"
