import hdbscan
import Unsupervised.clustering
import numpy
from sklearn.cluster import MiniBatchKMeans
import Scorer
from sklearn import mixture
import Optimizer
import copy


class Clusterer:
    def __init__(self, trajectory_2d):
        self.trajectory_2d = trajectory_2d
        self.labels = []

    def fit(self):
        raise NotImplementedError()


class HDBSCAN(Clusterer):
    def __init__(self, trajectory_2d, minimum_membership=10):
        self.minimum_membership = minimum_membership
        super().__init__(trajectory_2d)  # inherits constructor

    def fit(self):
        framework = hdbscan.HDBSCAN(min_cluster_size=self.minimum_membership)
        self.labels = framework.fit_predict(self.trajectory_2d)
        return self.labels


class KMeans(Clusterer):
    def __init__(self, trajectory_2d):
        self.trajectory_2d = trajectory_2d
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        cl = Unsupervised.clustering.Clustering()
        labels, _, _, _, _, = cl.ik_means(data=self.trajectory_2d)
        max_clusters = max(labels) + 1
        k_to_try = range(2, max_clusters + 2)  # Have to go +2 for slope
        scores = numpy.zeros(max_clusters)
        for k in k_to_try:
            clusterer = MiniBatchKMeans(
                n_clusters=k, n_init=100)
            labels = clusterer.fit_predict(self.trajectory_2d)
            scores[k - 2] = Scorer.Silhouette(data=self.trajectory_2d, labels=labels).evaluate()
        optimizer = Optimizer.Slope(scores, k_to_try)
        num_clusters = optimizer.minimize()
        clusterer = MiniBatchKMeans(n_clusters=num_clusters)
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.cluster_centers_
        return self.labels, self.centers


class GMM(Clusterer):
    def __init__(self, trajectory_2d, max_clusters=100):
        self.max_clusters = max_clusters
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        aic = numpy.zeros(self.max_clusters - 1)
        k_to_try = range(2, self.max_clusters + 1)
        for k in k_to_try:
            clusterer = mixture.GMM(
                n_components=k,
                covariance_type='tied',
                n_init=5
            )
            cluster = clusterer.fit(self.trajectory_2d)
            _ = clusterer.predict(self.trajectory_2d)
            aic[k - 2] = cluster.aic(self.trajectory_2d)
        optimizer = Optimizer.Optimizer(aic, k_to_try)
        num_clusters = optimizer.minimize()
        clusterer = mixture.GMM(
            n_components=num_clusters,
            covariance_type='tied'
        )
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.means_
        return self.labels, self.centers


class IMWKRescaled(Clusterer):
    def __init__(self, trajectory_2d, minkowski_weight=2):
        self.minkowski_weight = minkowski_weight
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        cl = Unsupervised.clustering.Clustering()
        [labels, _, _, _, _] = cl.ik_means(self.trajectory_2d, p=self.minkowski_weight)
        max_clusters = max(labels) + 1
        silhouette_averages = numpy.zeros(max_clusters - 1)
        k_to_try = numpy.arange(2, max_clusters + 1)
        for k in k_to_try:
            cl = Unsupervised.clustering.Clustering()
            data = copy.copy(self.trajectory_2d)
            [labels, centroids, weights, _, _] = cl.imwk_means(data, p=self.minkowski_weight, k=k)
        # Rescale the data
            for k1 in numpy.arange(0, max(labels)+1):
                data[labels == k1] = numpy.multiply(
                    self.trajectory_2d[labels == k1], numpy.tile(weights[k1], (numpy.sum(labels == k1), 1))
                )
                centroids[k1] = numpy.multiply(centroids[k1], weights[k1])
        # Apply Euclidean KMeans
            kmeans_clusterer = MiniBatchKMeans(n_clusters=k, n_init=5)
            kmeans_clusters = kmeans_clusterer.fit(data)
            labels = kmeans_clusters.labels_
            silhouette_averages[k - 2] = Scorer.Silhouette(labels=labels, data=data).evaluate()

        optimal_k = Optimizer.Optimizer(scores=silhouette_averages, parameter_list=k_to_try).maximize()
        # Do optimal clustering
        data = copy.copy(self.trajectory_2d)
        [labels, centroids, weights, _, _] = cl.imwk_means(data, p=self.minkowski_weight, k=optimal_k)
        # Rescale the data
        for k1 in numpy.arange(0, max(labels)+1):
            data[labels == k1] = numpy.multiply(data[labels == k1],
                                                numpy.tile(weights[k1], (numpy.sum(labels == k1), 1)))
            centroids[k1] = numpy.multiply(centroids[k1], weights[k1])

        # Apply Euclidean KMeans
        kmeans_clusterer = MiniBatchKMeans(n_clusters=optimal_k, n_init=5)
        kmeans_clusters = kmeans_clusterer.fit(data)
        self.labels = kmeans_clusters.labels_
        self.centers = kmeans_clusterer.cluster_centers_
        return self.labels, self.centers


class VBGMM(Clusterer):
    def __init__(self, trajectory_2d, max_clusters=100):
        self.max_clusters = max_clusters
        self.centers = []
        super().__init__(trajectory_2d)

    def fit(self):
        aic = numpy.zeros(self.max_clusters - 1)
        k_to_try = range(2, self.max_clusters + 1)
        for k in k_to_try:
            clusterer = mixture.VBGMM(
                n_components=k,
                covariance_type='tied',
                )
            cluster = clusterer.fit(self.trajectory_2d)
            _ = clusterer.predict(self.trajectory_2d)
            aic[k - 2] = cluster.aic(self.trajectory_2d)
        optimizer = Optimizer.Optimizer(aic, k_to_try)
        num_clusters = optimizer.minimize()
        clusterer = mixture.VBGMM(
            n_components=num_clusters,
            covariance_type='tied'
        )
        self.labels = clusterer.fit_predict(self.trajectory_2d)
        self.centers = clusterer.means_
        return self.labels, self.centers
