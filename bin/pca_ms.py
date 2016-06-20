import numpy as np
import mdtraj as md
import sklearn.cluster
import sklearn.decomposition
import argparse
import matplotlib.pyplot as plt
from itertools import cycle
import sys

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Run and score hdbscan clustering', add_help=False) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='name CA')
inputs.add_argument('-min', action='store', dest='min', help='minimum cluster membership',type=int,default=10)
inputs.add_argument('-o', action='store', dest='out_name',help='Output file',type=str,required=True)

# Parse into useful form
UserInput=parser.parse_args()

topology = UserInput.structure
trajectory = UserInput.trajectory
t = md.load(trajectory,top=topology)
sel = t.topology.select(UserInput.sel)
t = t.atom_slice(sel)

pca=sklearn.decomposition.PCA(n_components=10)
projection = pca.fit_transform(t.xyz.reshape(t.n_frames, t.n_atoms * 3))

for i in range(0,4):
    for j in range(i+1,5):
        temp = np.column_stack((projection[:,i], projection[:,j]))

        #bins = int(np.round(np.log2(t.xyz.shape[0])) + 1)
        q75,q25 = np.percentile(temp,[75, 25])
        irq = q75-q25
        bin_width = 2*irq/((t.xyz.shape[0])**(1./3))
        bins=np.ceil((np.max(temp[:,0])-np.min(temp[:,0]))/bin_width)
        #Biny=np.ceil((max(temp[:,1])-min(temp[:,1]))/bin_width)
        #bins=Binx*Biny
        print(bins)
        H, xedges, yedges = np.histogram2d(temp[:,0], temp[:,1], bins=bins)
        H[H==0] = np.nan
        E = -0.6 * np.log(H)
        Em = np.ma.masked_where(np.isnan(E), E)
        
        plt.pcolormesh(xedges, yedges, Em.T)
        plt.savefig(UserInput.out_name + '/bfree_energy{0}-{1}.png'.format(i,j))
        plt.close()
    
        clusterer = sklearn.cluster.MeanShift(n_jobs=-1,cluster_all=True)
        labels = clusterer.fit_predict(temp)
        cluster_centers = clusterer.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters_ = max(labels_unique)+1
        print(labels_unique)
        #print("number of estimated clusters : %d" % n_clusters_)
        
        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
        for k, col in zip(range(n_clusters_), colors):
            my_members = labels == k
            cluster_center = cluster_centers[k]
            plt.plot(projection[my_members, i], projection[my_members, j], col + '.')
            plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                     markeredgecolor='k', markersize=14)
            
        progress = "\r Motif calculation on Frame " + str(j) + " of " + str(j)  #status  
        sys.stdout.write(progress)
        sys.stdout.flush() #report status to terminal output
        plt.title('Estimated number of clusters: %d' % n_clusters_)
        plt.savefig(UserInput.out_name + '/bMeanShift{0}-{1}.png'.format(i,j))
        plt.close()
    
