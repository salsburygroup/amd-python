import hdbscan
import Unsupervised.clustering
import numpy
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import silhouette_score
from sklearn import mixture
import Optimizer


class Clusterer:
    def __init__(self, trajectory_2d):
        self.trajectory_2d = trajectory_2d

    def fit(self):
        raise NotImplementedError()


class HDBSCAN(Clusterer):
    def __init__(self, trajectory_2d, minimum_membership):
        self.minimum_membership = minimum_membership
        super().__init__(trajectory_2d)  # inherits constructor

    def fit(self):
        framework = hdbscan.HDBSCAN(min_cluster_size=self.minimum_membership)
        labels = framework.fit_predict(self.trajectory_2d)
        return labels


class KMeans(Clusterer):
    def fit(self):
        cl = Unsupervised.clustering.Clustering()
        _, max_clusters, _, _, _, = cl.ik_means(data=self.trajectory_2d)
        k_to_try = range(2, max_clusters + 2)  # Have to go +2 for slope
        scores = numpy.zeros(max_clusters)
        for k in k_to_try:
            clusterer = MiniBatchKMeans(n_clusters=k, n_init=5)
            labels = clusterer.fit_predict(self.trajectory_2d)
            scores[k - 2] = silhouette_score(self.trajectory_2d, labels)
        optimizer = Optimizer.Slope(scores, k_to_try)
        num_clusters = optimizer.minimize()
        clusterer = MiniBatchKMeans(n_clusters=num_clusters)
        labels = clusterer.fit_predict(self.trajectory_2d)
        return labels


class GMM(Clusterer):
    def __init__(self, trajectory_2d, max_clusters):
        self.max_clusters = max_clusters
        super().__init__(trajectory_2d)

    def fit(self):
        aic = numpy.zeros(self.max_clusters - 1)
        k_to_try = range(2, self.max_clusters + 1)
        for k in k_to_try:
            clusterer = mixture.GMM(
                n_components=k,
                covariance_type='tied',
                n_init=5)
            cluster = clusterer.fit(self.trajectory_2d)
            labels = clusterer.predict(self.trajectory_2d)
            aic[k - 2] = cluster.aic(self.trajectory_2d)
        return labels
