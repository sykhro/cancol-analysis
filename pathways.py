from collections import namedtuple
from enum import Enum
import os
import logging

Node = namedtuple('Node', ['name', 'parent', 'children', 'type'])

class PType(Enum):
    GENE = 0
    FAMILY = 1
    COMPLEX = 2
    PROCESS = 3

def __construct_node(toks):
    parent = None if toks[3] == "-1" else toks[3]
    return Node(toks[0], parent, {}, PType[toks[2]])

def __attempt_add_child(g, parent, nid, node):
    if parent in g:
        g[parent][0].children[nid] = [node, []]
        logging.debug(f"Added {node} to {g[parent]}")
        return True

    for subg in g:
        if __attempt_add_child(g[subg][0].children, parent, nid, node):
            return True

    return False

def __process_stash(g, stash):
    for toks in stash[:]:
        logging.debug(f"{toks[0]} :looking for parent in stash {toks[3]}")
        if toks[3] in g:
            g[toks[3]][0].children[toks[1]] = [__construct_node(toks), []]
            stash.remove(toks)
        else:
            if __attempt_add_child(g, toks[3], toks[1], __construct_node(toks)):
                stash.remove(toks)

def __find_node(g, nid):
    if not g:
        return {}

    logging.debug(f"Looking for node {nid}")
    logging.debug(f"Graph is {g}")
    if nid in g:
        logging.debug("Element located")
        return g[nid]

    logging.debug(f"Choices: {g.keys()}")
    for node in g:
        logging.debug(f"Now looking at {node}")
        res = __find_node(g[node][0].children, nid)
        if(res):
            return res

    return {}

# Transforms a pathwaymapper file into a graph
def parse_pathway(path: str):
    g = {}

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

            # In the first pass, add top-level nodes
            if toks[3] == "-1":
                g[toks[1]] = [__construct_node(toks), []]
            else:
                stash.append(toks)

            cur_str = pwfile.readline()

        while len(stash) > 0:
            __process_stash(g, stash)

        # Hit an empty line, edge definitions should follow after the csv header
        pwfile.readline()
        cur_str = pwfile.readline()
        while cur_str.strip() != "":
            toks = cur_str.split("\t")

            __find_node(g, toks[1])[1].append(toks[2])
            logging.debug(f"Edge added: {toks[1]}  ->{toks[2]}")

            cur_str = pwfile.readline()

    return (title, g)

def grouped_genes_size(pathway):
    return sum(node[0].type != PType.PROCESS for node in pathway.values())

def get_genes(pathway):
    genes = set()
    for node in pathway.values():
        if node[0].type == PType.GENE:
            genes.add(node[0].name)
        genes.update(get_genes(node[0].children))

    return genes

logging.basicConfig(level=logging.DEBUG, filename='pathway_parser.log')

pathways = []
for pw in os.listdir('./pathways'):
    pathway = parse_pathway('./pathways/' + pw)
    pathways.append(pathway)
    print(pathway[0], "has toplevel nodes:", len(pathway[1]),
          "grouped ", grouped_genes_size(pathway[1]))
    print("And the following genes", get_genes(pathway[1]))
    print()
