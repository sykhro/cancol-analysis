#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import networkx as nx
from dataclasses import dataclass, field

@dataclass
class Pathway:
    name: str
    graph: nx.Graph
    measures: dict = field(default_factory=dict)

    def calculate_measure(self, function, with_complexes=False):
        """Returns a series with weights for a specific function.
            Indexed by biomarker (*not* ID)"""
        # Cache results
        if function in self.measures and with_complexes in self.measures[function]:
            return self.measures[function][with_complexes]

        gene_names = nx.get_node_attributes(self.graph, "label")
        weights = function(self.graph)
        # Index by biomarker
        weights = {gene_names[k]: weights[k] for k, _ in weights.items()}

        # Separate paths as optimization
        self.measures[function] = {}
        if not with_complexes:
            self.measures[function][with_complexes] = pd.Series(weights)
            return self.measures[function][with_complexes]
        else:
            gene_weights = nx.get_node_attributes(self.graph, "famcomw")
            self.measures[function][with_complexes] = pd.Series(weights).mul(
                pd.Series(
                    {gene_names[k]: gene_weights[k] for k, v in gene_weights.items()}
                )
            )
            return self.measures[function][with_complexes]

    def get_genes(self):
        return nx.get_node_attributes(self.graph, "label").values()


def make_pathway_from_thres(threshold, coexpression):
    coexpression = coexpression.applymap(lambda x: 1 if abs(x) > threshold else 0)
    np.fill_diagonal(coexpression.values,0)
    graph = nx.Graph()
    for gene in coexpression:
        for other in coexpression.columns:
            if coexpression[gene][other] == 1:
                if gene not in graph:
                    graph.add_node(gene, label=gene)
                if other not in graph:
                    graph.add_node(other, label=other)

                graph.add_edge(gene, other)
    
    return Pathway("GPL570-{}".format(threshold), graph)
