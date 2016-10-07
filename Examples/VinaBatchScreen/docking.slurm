#!/bin/bash -l
#SBATCH --partition=medium
#SBATCH --account=salsburyGrp
#SBATCH --nodes=1
#SBATCH --mail-user=melvrl13@wfu.edu
#SBATCH --tasks-per-node=8
#SBATCH --mem=16gb
#SBATCH --time=01:00:00

module load vina/1.1.2-intel-2012

/home/salsbufr/babel/bin/babel -ismi ${l} --gen3d -f ${x} -l ${x} -opdb zinc.pdb

#/home/salsbufr/babel/bin/obminimize zinc${x}.pdb

/home/luy/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l zinc.pdb

vina --receptor ${r} --exhaustiveness 99 --ligand zinc.pdbqt  --center_x 23.7 --center_y 3.6  --center_z -1.2 --size_x 20 --size_y 20 --size_z 20 --out out.pdbqt > out.log

#If you're doing something comparative
#vina --receptor ${r} --exhaustiveness 99 --ligand zinc.pdbqt  --center_x 33.6 --center_y -8.3 --center_z 22.6 --size_x 28 --size_y 28 --size_z 28 --out F10.pdbqt > F10.log

# I wanted a certain range of energies. You should adjust the awk command for what you want and where you want it stored.
awk '{if (NF==4 && $2<10.0) print $2,$1,FILENAME}' out.log > ../energies/out_${o}.dat

#awk '{if (NF==4 && $2<10.0) print $2,$1,FILENAME}' F10.log > ../energies/F10_${o}.dat

module unload vina/1.1.2-intel-2012

rm -r ../${o}
