#!/usr/bin/env python
import subprocess
from shutil import move
import os
import re
import time
import string

squeue_command = 'squeue -u melvrl13'
cwd = os.getcwd()
#dir = os.path.dirname(__file__)
docking_helper = os.path.join(cwd, 'docking.slurm')
zinc_library = '23_t90.smi' #Smiles from FRS
with open(zinc_library, 'r+') as smiles_strings:
    num_drugs = sum(1 for line in smiles_strings)
with open(zinc_library, 'r+') as smiles_strings:
    smiles_text = smiles_strings.read()
    drugs = re.findall('(ZINC\d+)', smiles_text)
for drug_num in range(1, num_drugs + 1):
    squeue = subprocess.Popen(squeue_command, shell=True, stdout=subprocess.PIPE)
    out, err = squeue.communicate()
    out = out.splitlines()
    count = sum (1 for line in out)
    while count > 500:
        time.sleep(5)
        squeue = subprocess.Popen(squeue_command, shell=True, stdout=subprocess.PIPE)
        out, err = squeue.communicate()
        out = out.splitlines()
        count = sum (1 for line in out)
    compound = drugs[drug_num - 1]
    compound_dir = os.path.join(cwd, compound)
    if not os.path.exists(compound_dir):
        os.makedirs(compound_dir)
    submit_docking_command = (
            'sbatch --output=/dev/null --export=o=' +  compound + ',x=' + str(drug_num) + ',r=' + receptor + ',l=' ' + docking_helper
            )
    docking_job = subprocess.Popen(submit_docking_command, shell=True, cwd=compound_dir, stdout=subprocess.PIPE)


