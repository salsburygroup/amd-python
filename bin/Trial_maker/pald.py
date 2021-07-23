from scipy.spatial import distance
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import MDAnalysis
import argparse
import os

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='PALD clustering method', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
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

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='protein and resid 167:170')

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

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
if not os.path.exists(UserInput.out_dir):
    os.mkdir(UserInput.out_dir)
    
# Load trajectory
u = MDAnalysis.Universe(UserInput.structure, UserInput.trajectory)
selectedGrp=u.select_atoms(UserInput.sel)

# Compute distance matrix
N=len(u.trajectory)
d_traj=np.zeros((N, 3*len(selectedGrp)))
for i in range(N):
    u.trajectory[i]
    d_traj[i]=selectedGrp.positions.flatten()

d=distance.cdist(d_traj,d_traj,'euclidean')

# Save distance matrix
np.savetxt(os.path.join(UserInput.out_dir,'distance_matrix.txt'),d)

# Compute the matrix of partitioned local depths
A=np.zeros((N,N))
for i in range(N-1):
    for j in range(i+1,N):
        U=[k for k in range(N) if d[k,i]<=d[j,i] or d[k,j]<=d[i,j]]
        l=len(U)
        for p in range(l):
            a=U[p]
            if d[a,i]<d[a,j]:
                A[i,a]=A[i,a]+1/l
            elif d[a,i]==d[a,j]:
                A[i,a]=A[i,a]+0.5/l
                A[j,a]=A[j,a]+0.5/l
            else:
                A[j,a]=A[j,a]+1/l
C=A/(N-1)

# Save the matrix of partitioned local depths
np.savetxt(os.path.join(UserInput.out_dir,'partitioned_local_depths_matrix.txt'),d)

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
fig1.savefig(os.path.join(UserInput.out_dir, 'cmap1.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
#fig1.savefig(os.path.join(UserInput.out_dir, 'cmap1.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
#fig1.savefig(os.path.join(UserInput.out_dir, 'cmap1.pdf'), bbox_inches='tight')

# Compute and plot the local depths
ldepth=np.sum(C,axis=1)
fig2=plt.figure(2,figsize=(16,4))
plt.plot(ldepth, alpha=0.9, linewidth=0.5)
plt.xlabel('Frame(' + UserInput.timestep + ')', fontsize=14)
plt.ylabel('local depths', fontsize=14)
plt.title('local depths of ' + UserInput.title, fontsize=16)
fig2.savefig(os.path.join(UserInput.out_dir, 'ldepth.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig2.savefig(os.path.join(UserInput.out_dir, 'ldepth.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig2.savefig(os.path.join(UserInput.out_dir, 'ldepth.pdf'), bbox_inches='tight')
        
# Compute the cluster network Gs*
GS2=np.where(GS>=ps,GS,0)

# Plot the cluster network Gs*
fig3=plt.figure(3)
plt.pcolor(GS2,cmap="jet")
plt.colorbar()
plt.title('cluster network of ' + UserInput.title, fontsize=16)
fig3.savefig(os.path.join(UserInput.out_dir, 'cmap2.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
#fig3.savefig(os.path.join(UserInput.out_dir, 'cmap2.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
#fig3.savefig(os.path.join(UserInput.out_dir, 'cmap2.pdf'), bbox_inches='tight')

# Convert the matrix to network
E=nx.from_numpy_matrix(GS2)

# Plot network graph
def plot_graph(G, ax=None):
    D=list(nx.connected_components(G))
    fig4=plt.figure(4)
    color=['blue','red','green','yellow','orange','pink','lime','olive','purple','m','c','tab:blue','royalblue','ivory','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black']
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
    
    for i in range(len(D)):
        B=G.subgraph(D[i])
        weights = np.real([*nx.get_edge_attributes(B, 'weight').values()])
        nx.draw(B, pos, node_color=color[i], with_labels=True, edge_color=weights, edge_cmap=plt.cm.get_cmap('binary'), cmap=plt.cm.get_cmap('jet'), width=1, ax=ax, node_size=60, font_size=6)
        
    plt.title('network of ' + UserInput.title, fontsize=24)
    fig4.savefig(os.path.join(UserInput.out_dir, 'network.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
    fig4.savefig(os.path.join(UserInput.out_dir, 'network.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
    fig4.savefig(os.path.join(UserInput.out_dir, 'network.pdf'), bbox_inches='tight')

plot_graph(E)
