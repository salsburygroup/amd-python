import matplotlib.pyplot as plt
import networkx as nx
import mdtraj as md
import numpy as np
import argparse
import os

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='PALD clustering method', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-d',
                    action='store',
                    dest='data',
                    help='Input data',
                    type=str,
                    required=True)

inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-f',
                    action='store',
                    dest='frameFile',
                    help='file with frame of interest',
                    type=str,
                    default='noFile')

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='Closest mean distance between Na$\mathregular{^+}$ and 220s loop')

inputs.add_argument('-tm',
                    action='store',
                    dest='timestep',
                    help='time interval(ps)',
                    type=str,
                    required=True)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    required=True)

inputs.add_argument('-on',
                    action='store',
                    dest='outname',
                    help='Outname for data',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Load the matrix of partitioned local depths
C=np.load(UserInput.data)
N=len(C)

# Compute the particular strong number
ps=0.5*np.sum([C[i,i] for i in range(N)])/N

# Compute the cohesion network Gs
GS=np.zeros((N,N))
for i in range(N-1):
    for j in range(i+1,N):
        GS[j,i]=GS[i,j]=np.min([C[i,j],C[j,i]])

# Plot the cohesion network Gs
fig1=plt.figure(1)
plt.pcolor(GS,cmap="jet")
plt.colorbar()
plt.title('cohesion network of ' + UserInput.title, fontsize=16)
fig1.savefig(os.path.join(out_dir, 'cmap1.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
#fig1.savefig(os.path.join(out_dir, 'cmap1.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
#fig1.savefig(os.path.join(out_dir, 'cmap1.pdf'), bbox_inches='tight')

# Compute and plot the local depths
ldepth=np.sum(C,axis=1)
fig2=plt.figure(2,figsize=(16,4))
plt.plot(ldepth, alpha=0.9, linewidth=0.5)
plt.xlabel('Frame(' + UserInput.timestep + ')', fontsize=14)
plt.ylabel('local depths', fontsize=14)
plt.title('local depths of ' + UserInput.title, fontsize=16)
fig2.savefig(os.path.join(out_dir, 'ldepth.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig2.savefig(os.path.join(out_dir, 'ldepth.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig2.savefig(os.path.join(out_dir, 'ldepth.pdf'), bbox_inches='tight')
        
# Compute the cluster network Gs*
GS2=np.where(GS>=ps,GS,0)

# Plot the cluster network Gs*
fig3=plt.figure(3)
plt.pcolor(GS2,cmap="jet")
plt.colorbar()
plt.title('cluster network of ' + UserInput.title, fontsize=16)
fig3.savefig(os.path.join(out_dir, 'cmap2.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
#fig3.savefig(os.path.join(out_dir, 'cmap2.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
#fig3.savefig(os.path.join(out_dir, 'cmap2.pdf'), bbox_inches='tight')

# Convert the matrix to network
E=nx.from_numpy_matrix(GS2)

# Plot network graph
def plot_graph(G, ax=None):
    D=list(nx.connected_components(G))
    D.sort(key=len,reverse=True)
    fig4=plt.figure(4)
    color=['blue','red','green','yellow','orange','pink','lime','olive','purple','m','c','tab:blue','royalblue','ivory','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black']
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
    
    for i in range(len(D)):
        B=G.subgraph(D[i])
        if len(D[i])>1:
            weights = np.real([*nx.get_edge_attributes(B, 'weight').values()])
            nx.draw(B, pos, node_color=color[i], with_labels=True, edge_color=weights, edge_cmap=plt.cm.get_cmap('binary'), cmap=plt.cm.get_cmap('jet'), width=1, ax=ax, node_size=60, font_size=6)
        else:
            weights = np.real([*nx.get_edge_attributes(B, 'weight').values()])
            nx.draw(B, pos, node_color=color[-1], with_labels=True, edge_color=weights, edge_cmap=plt.cm.get_cmap('binary'), cmap=plt.cm.get_cmap('jet'), width=1, ax=ax, node_size=60, font_size=6)
            
    plt.title('network of ' + UserInput.title, fontsize=24)
    fig4.savefig(os.path.join(out_dir, 'network.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
    fig4.savefig(os.path.join(out_dir, 'network.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
    fig4.savefig(os.path.join(out_dir, 'network.pdf'), bbox_inches='tight')

plot_graph(E)

# Save indexes in each cluster
D=list(nx.connected_components(E))
D.sort(key=len,reverse=True)

if not os.path.exists(os.path.join(out_dir,'idx')):
    os.mkdir(os.path.join(out_dir,'idx'))

for i in range(len(D)):
    np.savetxt(os.path.join(out_dir,'idx','idx'+str(i)+'.txt'), list(D[i]), fmt='%i')

# Find structures in each cluster
outname=UserInput.outname
if not os.path.exists(os.path.join(out_dir, 'visualization', outname, 'cluster'+str(i))):
    os.makedirs(os.path.join(out_dir, 'visualization', outname, 'cluster'+str(i)))
   
#representative_frame = md.load_frame(UserInput.trajectory, UserInput.frameFile, top=UserInput.structure)
#representative_frame.save_pdb(os.path.join(out_dir, 'visualization', outname, 'cluster'+str(i), 'rep.pdb'))
#pca_well_sample = np.column_stack((np.zeros(type1[boolean_in_XY].size),type1[boolean_in_XY]))
#pca_well_sample=pca_well_sample.astype(int)

##temp = md.load(UsernIput.trajectory, top=UserInput.structure)
##frame = np.loadtxt(UserInput.frameFile, delimiter=',')
##traj = temp.slice(frame.astype(int)-1, copy = True)

##traj.save(os.path.join(out_dir, 'visualization', outname, 'cluster'+str(i), 'all.dcd'))

#coor.save_traj(inp, pca_well_sample, os.path.join(out_dir, 'visualization', outname, 'cluster'+str(i), 'all.dcd'))


