#!/bin/bash -l
#SBATCH --partition=medium
#SBATCH --account=salsburyGrp
#SBATCH --nodes=1
#SBATCH --mail-user=salsbufr@wfu.edu
#SBATCH --tasks-per-node=1
#SBATCH --mem=16gb
#SBATCH --time=01:00:00

#note mail-user is not turned on, but username needs to be there
# need to define your own cetner and box size
# add exhaustiveness back in if you want to search longer
#   currently using default search parameters
# currently set to find energies <-8.2 which would be micromola
# -12.3 would be nanomolar
# may need to change where energies are saved
# may need to change where structures are saved
module load vina/1.1.2-intel-2012

/home/salsbufr/babel/bin/babel -ismi ${l} --gen3d -f ${x} -l ${x} -opdb zinc.pdb

#/home/salsbufr/babel/bin/obminimize zinc${x}.pdb

/home/luy/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l zinc.pdb

vina --receptor ${r}  --ligand zinc.pdbqt  --center_x 23.7 --center_y 3.6  --center_z -1.2 --size_x 20 --size_y 20 --size_z 20 --cpu 1 --out out.pdbqt > out.log

#If you're doing two different searches 
#vina --receptor ${r}  --ligand zinc.pdbqt  --center_x 33.6 --center_y -8.3 --center_z 22.6 --size_x 28 --size_y 28 --size_z 28 --cpu 1 --out F10.pdbqt > F10.log

# I wanted a certain range of energies. You should adjust the awk command for what you want and where you want it stored.
#awk '{if (NF==4 && $2<-8.2) print $2,$1,FILENAME}' out.log > /home/salsbufr/DOCK/energies/out_${o}.dat
#fix and uncomment below if you are doing two different searches
#awk '{if (NF==4 && $2<-8.2) print $2,$1,FILENAME}' F10.log > /home/salsbufr/DOCK/energies/F10_${o}.dat

awk '{if (NF==4 && $2<-5.8 ){print $2,$1}}' out.log >> hits.dat
if [ -s hits.dat ]
then
sed '/ENDROOT/d' out.pdbqt > tmp
sed '/ENDBRANCH/d' tmp > tmp2
cp tmp2 /home/salsbufr/DOCK/energies/zinc_out_${o}.pdb
cp hits.dat /home/salsbufr/DOCK/energies/out_${o}.dat 
fi

module unload vina/1.1.2-intel-2012

rm -r ../${o}
