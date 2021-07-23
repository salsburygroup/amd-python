import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot RMSF', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-b',
                    action='store',
                    dest='bins',
                    help='num of bins',
                    type=int,
                    required=True)

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default=os.getcwd().split('/')[-1])

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix for text and png',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Load SASA
a=np.loadtxt(os.path.join(UserInput.out_dir, 'unbinding/Total_sasa.dat'))
b=np.loadtxt(os.path.join(UserInput.out_dir, 'double/Total_sasa.dat'))
c=np.loadtxt(os.path.join(UserInput.out_dir, 'outer/Total_sasa.dat'))
d=np.loadtxt(os.path.join(UserInput.out_dir, 'inner/Total_sasa.dat'))

count1,bin1=np.histogram(a,bins=UserInput.bins)
count2,bin2=np.histogram(b,bins=UserInput.bins)
count3,bin3=np.histogram(c,bins=UserInput.bins)
count4,bin4=np.histogram(d,bins=UserInput.bins)

# Plot 
fig1=plt.figure(1,figsize=(8,4))
l1=plt.hist(bin1[:-1],bins=bin1,weights=count1/np.sum(count1),histtype='step',color='r',alpha=0.9,label='unbinding')
l2=plt.hist(bin2[:-1],bins=bin2,weights=count2/np.sum(count2),histtype='step',color='g',alpha=0.9,label='double')
l3=plt.hist(bin3[:-1],bins=bin3,weights=count3/np.sum(count3),histtype='step',color='b',alpha=0.9,label='outer')
l4=plt.hist(bin4[:-1],bins=bin4,weights=count4/np.sum(count4),histtype='step',color='y',alpha=0.9,label='inner')
plt.xlabel('Total SASA of interest (nm)^2')
plt.ylabel('Probability')
plt.legend()
plt.title('Distribution of Total SASA of {0}'.format(UserInput.title))
fig1.savefig(os.path.join(UserInput.out_dir, 'hist_sasa.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig1.savefig(os.path.join(UserInput.out_dir, 'hist_sasa.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(os.path.join(UserInput.out_dir, 'hist_sasa.pdf'), bbox_inches='tight')

