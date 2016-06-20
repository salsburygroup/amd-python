#!/usr/bin/env python
import subprocess
from shutil import move
import os
import re

zinc_library = '23_t90.smi' #Smiles from FRS
with open(zinc_library, 'r+') as smiles_strings:
    num_drugs = sum(1 for line in smiles_strings)
for drug_num in range(1, num_drugs + 1):
    babel_command = (
            '/home/salsbufr/babel/bin/babel -i smi /home/salsbufr/babel/23_t90.smi -O zinc{0}.pdb -f {0} -l {0} --gend3d'.format(drug_num)
            )
    subprocess.call(babel_command, shell=True)
    pdb = 'zinc{0}.pdb'.format(drug_num)
    with open(pdb,'r') as fobj:
        text = fobj.read()
    compound = re.findall('^COMPND\s+(ZINC\d+)', text)[0]
    cwd = os.getcwd()
    compound_dir = os.path.join(cwd,compound)
    if not os.path.exists(compound_dir):
        os.makedirs(compound_dir)
    pdb_path = os.path.join(cwd, pdb)
    pdb_resting_place = os.path.join(compound_dir, pdb)
    pdbqt_resting_place = pdb_resting_place + 'qt'
    move(pdb_path,pdb_resting_place)
    ligand_prep_command = (
            '/home/luy/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l '
            + pdb_resting_place + ' -o' + pdbqt_resting_place
            )
    subprocess.call(ligand_prep_command, shell=True)

    F10site_helper = os.path.join(cwd, 'F10_site_docking.slurm')
    # Run F10 site docking
    submit_F10site_command = (
            'sbatch --export=x=' + pdbqt_resting_place + ' ' + F10site_helper
            )
    F10site_job = subprocess.Popen(submit_F10site_command, shell=True, cwd=compound_dir, stdout=subprocess.PIPE)
    out, err = F10site_job.communicate()
    jobid = out.split(' ')[3]
    with open('F10_site_jobs.txt', 'a+') as F10_site_log:
        F10_site_log.write(compound + '\t' + jobid)

    Acidsite_helper = os.path.join(cwd, 'Acid_site_docking.slurm')
    # Run acid site docking
    submit_Acidsite_command = (
            'sbatch --export=x=' + pdbqt_resting_place + ' ' + Acidsite_helper
            )
    print(submit_Acidsite_command)
    Acidsite_job = subprocess.Popen(submit_Acidsite_command, shell=True, cwd=compound_dir, stdout=subprocess.PIPE)
    out, err = Acidsite_job.communicate()
    jobid = out.split(' ')[3]
    with open('Acid_site_jobs.txt', 'a+') as Acid_site_log:
        Acid_site_log.write(compound + '\t' + jobid)

