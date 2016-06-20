#!/bin/bash -l
#SBATCH --partition=small
#SBATCH --account=salsburyGrp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=16gb
#SBATCH --time=0-01:00:00

module load vina/1.1.2-intel-2012

vina --receptor ${r} --ligand ${x}  --center_x ${cx} --center_y ${cy}  --center_z ${cz} --size_x ${sx} --size_y ${sy} --size_z ${sz}

module unload vina/1.1.2-intel-2012