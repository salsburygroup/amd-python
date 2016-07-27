#!/bin/bash -l
#SBATCH --partition=medium
#SBATCH --mail-user=melvrl13@wfu.edu
#SBATCH --account=salsburygrp
#SBATCH --nodes=1
#SBATCH --tasks-per-node=6
#SBATCH --mem=120gb
#SBATCH --time=120:00:00


python $s -top $str -traj $t -sel $a -method $m -o $o 



