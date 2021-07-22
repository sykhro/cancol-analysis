from enum import Enum

import networkx as nx


# Nodes in a pathway can be genes or clusters of genes
# Genes inside a cluster are CHILDGENEs so that they can be treated
# independently from their parent nodes
class PathwayElement(Enum):
    GENE = 1
    COMPLEX = 2
    FAMILY = 3
    PROCESS = 4
    CHILDGENE = 5


# Transforms a pathwaymapper file into a graph
def pathway_to_graph(path: str):
    g = nx.Graph()

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
        genes = []
        complexes = {}
        while cur_str.strip() != "":
            toks = cur_str.split("\t")

            # An edge can connect to more than one vertex because of gene complexes or families,
            # this is modelled with an intermediate node
            ptype = (
                PathwayElement.CHILDGENE if toks[3] != "-1" else PathwayElement[toks[2]]
            )
            g.add_node(toks[1], label=toks[0], ptype=ptype)
            if ptype == PathwayElement.CHILDGENE:
                print(toks[0] + " " + toks[1] + " -> " + toks[3])
                g.add_edge(toks[1], toks[3])
            cur_str = pwfile.readline()

        # Hit an empty line, edge definitions should follow after the csv header
        pwfile.readline()
        cur_str = pwfile.readline()
        while cur_str.strip() != "":
            toks = cur_str.split("\t")
            if toks[1] in g.nodes() and toks[2] in g.nodes():
                print(toks[1] + " -> " + toks[2])
                g.add_edge(toks[1], toks[2], eid=toks[0])

            cur_str = pwfile.readline()

    return (title, g)


# Counts genes, complexes and families
def total_gene_groups(pathway: nx.Graph) -> int:
    gattrs = nx.get_node_attributes(pathway, "ptype")
    return sum(
        value in [PathwayElement.GENE, PathwayElement.FAMILY, PathwayElement.COMPLEX]
        for value in gattrs.values()
    )


# Counts all the genes present in the pathway
def total_genes(pathway: nx.Graph) -> int:
    gattrs = nx.get_node_attributes(pathway, "ptype")
    return sum(
        value in [PathwayElement.GENE, PathwayElement.CHILDGENE]
        for value in gattrs.values()
    )


def get_gene_names(pathway: nx.Graph):
    return nx.get_node_attributes(pathway, "label").values()


res = pathway_to_graph("./pathways/HIPPO.txt")
output = nx.nx_agraph.to_agraph(res[1])
output.layout("dot")
output.draw("test.png")

print("Total genes: " + str(total_gene_groups(res[1])))
