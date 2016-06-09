import subprocess


class Prepper():
    def __init__(self, in_file, out_file):
        self.in_file = in_file
        self.out_file = out_file

    def prep(self):
        raise NotImplementedError


class Ligand(Prepper):
    def prep(self):
        ligand_prep_command = (
            '/home/luy/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l '
            + self.in_file + ' -o' + self.out_file
        )
        subprocess.call(ligand_prep_command, shell=True)


class Receptor(Prepper):
    def prep(self):
        raise NotImplementedError
