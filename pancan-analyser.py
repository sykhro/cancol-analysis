#!/usr/bin/python3

import analysis_nx as anx
import pathways_nx as pnx
import networkx as nx
import logging as log
import pandas as pd
import os

PATHWAYS_DIRECTORY = "./pathways/"
log.basicConfig(level=log.INFO)

pathways = [
    pnx.pathway_to_nx(PATHWAYS_DIRECTORY + pw) for pw in os.listdir(PATHWAYS_DIRECTORY)
]
pathways.sort(key=lambda pw: pw.name)
log.info(f"Loaded {len(pathways)} pathways from {PATHWAYS_DIRECTORY}")

patients_log = pd.read_csv("TRIBE2_db.csv")
log.info(f"Loaded {len(patients_log)} patients")
mutations_data = pd.read_csv("TRIBE2_seq_res.csv")
log.info("Loaded gene mutations")
log.debug(f"Got {len(mutations_data)}")

arm0_df_indeg = anx.process_patients_with_f(
    patients_log[patients_log["arm"] == 0]["PatientFirstName"],
    nx.in_degree_centrality,
    pathways,
    mutations_data,
)
print(arm0_df_indeg.describe())

log.info("Goodbye")
