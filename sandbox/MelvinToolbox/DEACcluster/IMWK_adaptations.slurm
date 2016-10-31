#!/bin/bash -l
#SBATCH --partition=medium
#SBATCH --mail-user=melvrl13@wfu.edu
#SBATCH --account=salsburygrp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=120gb
#SBATCH --time=120:00:00


python /home/melvrl13/Repos/amd-python/MelvinToolbox/DEACcluster/IMWK_adaptations.py -top -traj -sel -o