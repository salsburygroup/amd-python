#!/bin/bash -l
#SBATCH --partition=small
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=melvrl13@wfu.edu
#SBATCH --account=salsburyGrp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=16gb
#SBATCH --time=0-03:00:00


module load vina/1.1.2-intel-2012

vina --receptor {x} --ligand ${y}  --center_x 33.6 --center_y -8.4  --center_z 22.6 --size_x 13 --size_y 10 --size_z 14 --out F10_site.pdbqt --log F10site.log

module unload vina/1.1.2-intel-2012

exit



