#!/bin/bash
#SBATCH -J MAFFT_align
#SBATCH -A uoa00029         # Project Account
#SBATCH --time=1:00:00     # Walltime
#SBATCH --mem-per-cpu=1024  # memory/cpu (in MB)
#SBATCH --cpus-per-task=8
###SBATCH --gres=gpu ## I could not get this part working in Apr-2019, and commented it out
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
# #SBATCH --mail-user=sknu003@aucklanduni.ac.nz
#SBATCH -o stdout_mafft_align.txt
#SBATCH -e stderr_mafft_align.txt

f01="inp02_output01_accn_ntseq_Mnemiopsis_its1its2.txt"
f02="inp02_algn_Mnelei_18s_09.fas" 
nf01=$(echo $f01 | sed 's/\.txt/.fas/g' )
cat $f01 | sed 's/^/>/g' | \
	sed 's/>/>Mnemiopsis_leidyi_/g' | \
	sed 's/\t/\n/g' > $nf01

FILE="algn_Mnelei_18s_10"
cat "$nf01" > "$FILE"
cat "$f02" >> "$FILE"
# before you start make sure you have MAFFT installed on your unix-system, 
# and that you have the version with all extensions. including linsi. 
# You can check if you have linsi in MAFFT available by typing: linsi -help

# start this bash submission script externally by typing:
# sbatch slurm_submit_triplefin_04_0"$i"_MAFFT.sh

#module spider MAFFT
# --------------------------------------
# MAFFT:
# --------------------------------------
#    Versions:
#        MAFFT/7.164-goolf-1.5.14-with-extensions
#        MAFFT/7.273-foss-2015a

##module load MAFFT/7.273-foss-2015a
module load MAFFT/7.310-gimkl-2017a
#module load MAFFT/7.429-gimkl-2018b

# linsi --help

#get the present working directory and make this a variable that you can call
WD=$(pwd)


linsi $FILE > $FILE.aligned.fasta.txt
