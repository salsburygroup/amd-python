from clustering import Clustering
import numpy as np
import MDAnalysis as md
import sklearn.metrics
import matplotlib.pyplot as plt
import seaborn
import random
import argparse

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Calculate RGYR time series', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-top', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-traj', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection',type=str,default='all')
inputs.add_argument('-o', action='store', dest='out_name',help='Output file',type=str,required=True)

#Parse into useful form
UserInput=parser.parse_args()

# For calculating gap index, following https://datasciencelab.wordpress.com/2013/12/27/finding-the-k-in-k-means-clustering/

topology = UserInput.structure
trajectory = UserInput.trajectory
t = md.Universe(topology,trajectory)
sel = t.select_atoms(UserInput.sel)

# Helper functions for gap index. 
# Gap index compares to a random distribution, so we need some functions for generating and clustering such a distribution.
def Wk(mu, clusters):
    K = len(mu)
    return sum([np.linalg.norm(mu[i]-c)**2/(2*len(c)) \
               for i in range(K) for c in clusters[i]])

def bounding_box(X):
    xmin, xmax = min(X,key=lambda a:a[0])[0], max(X,key=lambda a:a[0])[0]
    ymin, ymax = min(X,key=lambda a:a[1])[1], max(X,key=lambda a:a[1])[1]
    return (xmin,xmax), (ymin,ymax)

def cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters
 
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu
 
def has_converged(mu, oldmu):
    return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))

def find_centers(X, K):
    # Initialize to K random centers
    oldmu = random.sample(X, K)
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
    return(mu, clusters)


# Format trajectory 
temp = t.trajectory.timeseries(sel,format=u'fac')
frames = len(t.trajectory)
atoms = len(sel.atoms)
data = temp.reshape((frames,atoms*3))

np.seterr(all='raise')
cl = Clustering()
#
#data = np.genfromtxt("iris.csv", delimiter=',')
#
# standardize data
data = cl.my_math.standardize(data)

# Trying to find the optimal p
(xmin,xmax), (ymin,ymax) = bounding_box(data)
p_to_try = np.arange(1.1,5.1,0.1)
silhouette_scores = np.zeros(p_to_try.size)
gap_scores = np.zeros(p_to_try.size)
gap_stde = np.zeros(p_to_try.size) 
for q in range(0, p_to_try.size):
    print(q)
    [u, centroids, weights, ite, dist_tmp] = cl.imwk_means(data, p_to_try[q])
    # Calculate silhouette score
    silhouette_scores[q] = sklearn.metrics.silhouette_score(data,u)
    # Calculate gap score
    k = max(u)
    mu, clusters = find_centers(data,k)
    Wks = np.log(Wk(mu, clusters))
    # Create B reference datasets
    B = 10 
    BWkbs = np.zeros(B)
    for i in range(B):
        Xb = []
        for n in range(len(data)):
            Xb.append([random.uniform(xmin,xmax),
                          random.uniform(ymin,ymax)])
        Xb = np.array(Xb)
        mu, clusters = find_centers(Xb,k)
        BWkbs[i] = np.log(Wk(mu, clusters))
    Wkbs = sum(BWkbs)/B
    sk = np.sqrt(sum((BWkbs-Wkbs)**2)/B)
    sk = sk*np.sqrt(1+1/B)
    gap_scores[q] = Wkbs
    gap_stde[q] = sk



np.savetxt(UserInput.out_name + '/silhouette_scores.ext',silhouette_scores)
np.savetxt(UserInput.out_name + '/ptrials_gap.txt',gap_scores)
np.savetxt(UserInput.out_name + '/ptrials_gap_stde.txt',gap_stde)

plt.figure()
plt.plot(p_to_try,silhouette_scores)
plt.xlabel('Minkowski Weight')
plt.ylabel('Silhouette Score')
plt.title('Finding Minkowski Weight')
plt.savefig(UserInput.out_name + '/silhouette.png')
plt.figure()
plt.errorbar(p_to_try,gap_scores, yerr=gap_stde)
plt.xlabel('Minkowski Weight')
plt.ylabel('Gap Score')
plt.title('Finding Minkowski Weight')
plt.savefig(UserInput.out_name + '/gap.png')

