import subprocess

class converter():
    def __init__(self, in_file, out_file):
        self.in_file = in_file
        self.out_file = out_file

    def convert(self):
        raise NotImplementedError

class smi_to_pdb(converter):
    def __init__(self, in_file, out_file, smi_num=1):
        self.smi_num = smi_num
        super().__init__(in_file, out_file)

    def convert(self):
        babel_command = (
            '/home/salsbufr/babel/bin/babel -i smi {0} -O {1}_{2}.pdb -f {2} -l {2} --gend3d'.format(
                self.in_file, self.out_file, self.smi_num
            )
        )
        subprocess.call(babel_command, shell=True)