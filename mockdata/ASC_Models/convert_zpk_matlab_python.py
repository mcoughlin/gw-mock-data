'''
 Code to convert zpk in .mat to .npz format readable in python

python convert_zpk_matlab_python.py -m Nikhil_Data/SUS_model_PUM_Yaw_to_TST_Yaw.mat -p Nikhil_Data/SUS_model_PUM_Yaw_to_TST_Yaw.npz
'''

import numpy as np
import argparse
from os import system
from scipy.io import loadmat


current_dir = system('pwd')



class helpfulParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = helpfulParser()

parser.add_argument('-m', '--mat', default='ZPK.mat', type=str, nargs='?',
                    help='Matlab zpk filename in .mat format. ZPK.mat ')

parser.add_argument('-p', '--npz', default='ZPK.npz', type=str, nargs='?',
                    help='python  zpk filename in .npz format. Defaults to ZPK.npz')

parser.add_argument('-f', '--filepath', default=current_dir, type=str, nargs='?',
                    help='filepath to read/write data. Defaults to current directory ')

# Get parameters into global namespace
args   = parser.parse_args()

filepath = args.filepath
matlab_filename = args.mat
python_filename = args.npz



zpk = loadmat(matlab_filename)
z = zpk['z'].reshape(len(zpk['z']),)
p = zpk['p'].reshape(len(zpk['p']),)
k = zpk['k'].reshape(len(zpk['k']),)

np.savez(python_filename,z=z,p=p,k=k)

print('{} file has been converted to {}'.format(matlab_filename, python_filename))
