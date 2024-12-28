"""
Mapping Protein-Protein Interaction Networks
Robin Anwyl
December 2024
Description: Maps protein-protein interaction networks from the STRING
database.
"""

import requests
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
pd.set_option('display.max_columns', None)

class PPINetworkFile:
    """
    A class that represents a protein-protein interaction network from the
    STRING database.
    """
    def __init__(self, filename, k_):
        """
        Initializes class attributes.
        :param str filename: STRING .tsv file (with reciprocal edges) for
            PPI network
        :param int k_: number of clusters
        """
        self.graph = nx.Graph()
        self.file = filename
        self.k = k_
        self.cluster_labels = None

    def construct_graph(self):
        # Read in file
        data = pd.read_csv(self.file, sep='\t')
        data.columns = ["node1", "node2", "node1_string_id", "node2_string_id",
                        "neighborhood_on_chromosome", "gene_fusion",
                        "phylogenetic_cooccurrence", "homology",
                        "coexpression",
                        "experimentally_determined_interaction",
                        "database_annotated", "automated_textmining",
                        "combined_score"]
        # Add edges between proteins with edge weights = scores
        for index, row in data.iterrows():
            protein1 = row["node1"]
            protein2 = row["node2"]
            score = row["combined_score"]
            self.graph.add_edge(protein1, protein2, weight=score)

    def spectral_clustering(self):
        # Compute matrices (numpy arrays)
        # Adjacency matrix
        adjacency_m = nx.to_numpy_array(self.graph)
        # Degree matrix
        degrees = np.sum(adjacency_m, axis=1)
        degree_m = np.diag(degrees)
        # Laplacian matrix
        laplacian_m = degree_m - adjacency_m

        # Eigen-analysis on Laplacian
        eigenvalues, eigenvectors = np.linalg.eigh(laplacian_m)
        # Construct matrix formed by first k eigenvectors
        #   (corresponding to the k smallest eigenvalues of laplacian_m)
        eig_m = eigenvectors[:, :self.k]
        # Perform k-means clustering on eigenvectors
        cluster_labels = KMeans(n_clusters=self.k).fit_predict(eig_m)

        # Visualize
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(self.graph)
        nx.draw(self.graph, pos, with_labels=True, node_color=cluster_labels,
                cmap=plt.cm.get_cmap('Pastel1'), node_size=500)
        plt.title("Protein-Protein Interaction Network")
        plt.show()

    def map_network(self):
        self.construct_graph()
        self.spectral_clustering()


class PPINetworkURL:
    """
    A class that represents a protein-protein interaction network from the
    STRING database.
    """
    def __init__(self, networkId, k_):
        self.graph = nx.Graph()
        url = "https://string-db.org/api/tsv/network?networkId=" + networkId
        file = requests.get(url).text
        lines = file.split("\n")
        data = [l.split("\t") for l in lines]
        self.df = pd.DataFrame(data[1:-1], columns=data[0])
        self.k = k_
        self.cluster_labels = None

    def construct_graph(self):
        # Add edges between proteins with edge weights = scores
        for index, row in self.df.iterrows():
            protein1 = row["preferredName_A"]
            protein2 = row["preferredName_B"]
            score = row["score"]
            self.graph.add_edge(protein1, protein2, weight=score)

    def spectral_clustering(self, title, cmap_name):
        # Compute matrices (numpy arrays)
        # Adjacency matrix
        adjacency_m = nx.to_numpy_array(self.graph)
        # Degree matrix
        degrees = np.sum(adjacency_m, axis=1)
        degree_m = np.diag(degrees)
        # Laplacian matrix
        laplacian_m = degree_m - adjacency_m

        # Eigen-analysis on Laplacian
        eigenvalues, eigenvectors = np.linalg.eigh(laplacian_m)
        # Construct matrix formed by first k eigenvectors
        #   (corresponding to the k smallest eigenvalues of laplacian_m)
        eig_m = eigenvectors[:, :self.k]
        # Perform k-means clustering on eigenvectors
        cluster_labels = KMeans(n_clusters=self.k).fit_predict(eig_m)

        # Visualize
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(self.graph)
        nx.draw(self.graph, pos, with_labels=True, node_color=cluster_labels,
                cmap=plt.get_cmap(cmap_name), node_size=500)
        plt.title(title)
        plt.show()

    def map_network(self, title, cmap_name):
        self.construct_graph()
        self.spectral_clustering(title, cmap_name)


#G1 = PPINetworkFile("FAA4_S_cerevisiae.tsv", 3)
#G1.map_network()

#G2 = PPINetworkURL("bzZuIbI0NG4C", 3)
#G2.map_network("FAA4", "Pastel1")