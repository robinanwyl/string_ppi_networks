"""
old code to construct graph from locally saved .tsv file
"""

import pandas as pd
import networkx as nx
pd.set_option('display.max_columns', None)


class PPINetworkLocalFile:
    """
    A class that represents a protein-protein interaction network from the
    STRING database. Input: tab-separated values (.tsv) file for the PPI
    network, with reciprocal edges, downloaded from the STRING database and
    saved in the same directory as this file.
    """
    def __init__(self, filename):
        """
        Initializes class attributes.
        :param str filename: STRING .tsv file for PPI network, with reciprocal
            edges
        """
        self.graph = nx.Graph()
        self.file = filename

    def construct_graph(self):
        """
        Reads in a STRING .tsv file for a protein-protein interaction network
        and constructs a representative NetworkX graph.
        """
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


def main():
    faa4 = PPINetworkLocalFile("FAA4_S_cerevisiae.tsv")