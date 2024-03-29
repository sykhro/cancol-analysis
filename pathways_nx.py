import logging as log
from collections import namedtuple
from dataclasses import dataclass, field

import networkx as nx
import pandas as pd

AliasItem = namedtuple("AliasItem", ["parent", "genes"])
GeneElem = namedtuple("GeneElem", ["name", "id", "parent"])


@dataclass
class Pathway:
    name: str
    graph: nx.DiGraph
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


def __make_node_aliases(data: list[str]):
    """Alias a genes ID to their families
    in order to build edges between them"""
    famcom = {}
    elems = [tokens for tokens in data if tokens[2] in ["FAMILY", "COMPLEX"]]
    # Add all (gene) containers first
    for tokens in elems:
        famcom[tokens[1]] = AliasItem(tokens[3], [])

    log.debug(famcom)
    elems = [tokens for tokens in data if tokens[2] == "GENE"]
    for tokens in elems:
        # Add gene to its parent
        famcom[tokens[3]].genes.append(GeneElem(tokens[0], tokens[1], tokens[3]))

    return famcom


def pathway_to_nx(path: str) -> Pathway:
    """Parses a pathwaymapper file, returning a NetworkX graph"""
    g = nx.DiGraph()

    with open(path) as pwfile:
        # First line in file is pathway name
        title = pwfile.readline().strip()
        # Second line is empty
        pwfile.readline()
        # Third line is pathway description, ignore for now
        pwfile.readline()
        # Fourth line is empty, fifth line is genes header
        pwfile.readline()
        pwfile.readline()

        # Now parse the genes
        cur_str = pwfile.readline()
        stash = []
        while cur_str.strip() != "":
            toks = cur_str.split("\t")

            # We are interested in the first 4 columns:
            # NAME, ID, TYPE, PARENT_ID

            # In the first pass, add top-level nodes by id;
            # store their name for convenience
            if toks[3] == "-1" and toks[2] == "GENE":
                g.add_node(toks[1], label=toks[0], famcomw=1)
                log.debug(f"Node added: {toks[0]}, {toks[1]}")
            else:
                stash.append(toks)
                log.debug(f"Process later: {toks[0]}, {toks[1]}")

            cur_str = pwfile.readline()

        aliases = __make_node_aliases(stash)
        # Add all genes to the top-level graph
        for item in aliases:
            to_add = [g for g in aliases[item].genes if not g.id in aliases]
            for gene in to_add:
                parent = aliases[gene.parent]
                famcomsize = len(parent.genes)
                # Take into account the size of all the ancestor containers
                # e.g.: if gene A is in a complex of size 2 inside of a complex of size 3,
                # it should have a final weight of 1/6
                while parent.parent != "-1":
                    parent = aliases[parent.parent]
                    famcomsize *= len(parent.genes)
                g.add_node(
                    gene.id,
                    label=gene.name,
                    famcomw=(1 if famcomsize == 0 else 1 / famcomsize),
                )
                log.debug(f"Node added: {gene.name}, {gene.id}, 1/{famcomsize}")

        # Hit an empty line, edge definitions should follow after the csv header
        pwfile.readline()
        cur_str = pwfile.readline()
        while cur_str.strip() != "":
            toks = cur_str.split("\t")

            source_nodes = []
            if toks[1] in aliases:
                deepen = [toks[1]]
                while deepen:
                    gid = deepen.pop(0)
                    if gid in aliases:
                        deepen.extend([g.id for g in aliases[gid].genes])
                    else:
                        source_nodes.append(gid)
            else:
                source_nodes = [toks[1]]

            target_nodes = []
            if toks[2] in aliases:
                deepen = [toks[2]]
                while deepen:
                    gid = deepen.pop(0)
                    if gid in aliases:
                        deepen.extend([g.id for g in aliases[gid].genes])
                    else:
                        target_nodes.append(gid)
            else:
                target_nodes = [toks[2]]

            for source in source_nodes:
                if source not in g:
                    continue

                for target in target_nodes:
                    # Don't link to processes for now
                    if target not in g:
                        continue
                    g.add_edge(source, target, label=toks[3])

            log.debug(f"Edge added: {source_nodes}-{toks[3]}->{target_nodes}")

            cur_str = pwfile.readline()

    return Pathway(title, g)
