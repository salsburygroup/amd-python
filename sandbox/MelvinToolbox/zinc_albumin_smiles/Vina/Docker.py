import subprocess
import os
import re

class Docker:
    def __init__(self, ligand, receptor, center_x, center_y, center_z, size_x, size_y, size_z, out_folder):
        self.ligand = ligand
        self.receptor = receptor
        self.center_x = center_x
        self.center_y = center_y
        self.center_z = center_z
        self.size_x = size_x
        self.size_y = size_y
        self.size_z = size_z
        self.out_folder = out_folder
        self.job = []

    def dock(self):
        raise NotImplementedError


class Shell(Docker):
    def dock(self):
        shell_helper = os.path.join(dir, 'src', 'shell_helper.sh')
        submit_command = (
            [shell_helper, self.receptor, self.ligand, self.center_x, self.center_y, self.center_z,
            self.size_x, self.size_y, self.size_z]
        )
        self.job = subprocess.Popen(submit_command, shell=True, cwd=self.out_folder, stdout=subprocess.PIPE)


class Torque(Docker):
    def dock(self):
        raise NotImplementedError


class Slurm(Docker):
    def dock(self):
        dir = os.path.dirname(__file__)
        slurm_helper = os.path.join(dir, 'src', 'slurm_helper.slurm')
        docking_variables = ('l=' + self.ligand
                             + ',r=' + self.receptor
                             + ',cx=' + self.center_x
                             + ',cy=' + self.center_y
                             + ',cz=' + self.center_z
                             + ',sx=' + self.size_x
                             + ',sy=' + self.size_y
                             + ',sz=' + self.size_z
                             )
        submit_command = (
            'sbatch --export=' + docking_variables + ' ' + slurm_helper
        )
        self.job = subprocess.Popen(submit_command, shell=True, cwd=self.out_folder, stdout=subprocess.PIPE)
        out, _ = self.job.communicate()
        jobid = out.split(' ')[3]
        return jobid

class EnergyOnly(Shell):
    def get_dG(self):
        self.dock()
        self.job.wait()
        for file in os.listdir(self.out_folder):
            if re.match(".*_out.pdbqt", file):
                model_details = os.path.join(self.out_folder, file)
        with open(model_details, 'rb') as f2:
            for line in f2:  # Loop over lines in each file
                # creates a list of single tuple [(model,engergy0]
                model = re.findall(r'\s+([1-9])\s+(-[0-9]+\.[0-9]+)',line)
        return model

