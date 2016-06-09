#!/bin/bash -l

module load vina/1.1.2-intel-2012

vina --receptor $1 --ligand $2  --center_x $3 --center_y $4  --center_z $5 --size_x $6 --size_y $7 --size_z $8

module unload vina/1.1.2-intel-2012