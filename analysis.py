#!/usr/bin/env python

import os

import numpy as np
import pandas
import seaborn as sns

from pathways import *


def calculate_patient_mutations(pid, seq_data, pathways):
    patient_data = seq_data[
        (seq_data["PatientFirstName"] == pid)
        & (seq_data["Technology"] == "NGS Q3")
        # We only care about variants and pathogenic mutations
        & (seq_data["TestResult"].isin(["variantdetected", "Mutated, Pathogenic"]))
    ]

    patient_data = patient_data[["Biomarker", "NGS_PercentMutated"]]

    # Realistically, should never happen
    if patient_data.empty:
        return {}

    results = {}
    for pw in pathways:
        pw_genes = get_genes(pw[1])
        pathway_mutations = patient_data[patient_data["Biomarker"].isin(pw_genes)]
        if pathway_mutations.empty:
            results[pw[0]] = np.float64(0.0)
            continue

        perc_mutation = pathway_mutations.groupby("Biomarker").max()[
            "NGS_PercentMutated"
        ].sum() / grouped_genes_size(pw[1])
        results[pw[0]] = perc_mutation

    return results


def process_patients(patients, mutations_data, pathways):
    results = {}
    for patient in patients:
        results[patient] = calculate_patient_mutations(
            patient, mutations_data, pathways
        )

    return (
        pandas.DataFrame.from_dict(results, orient="index")
        .rename_axis("PatientFirstName")
        .reset_index()
    )


if __name__ == "__main__":
    pathways = []
    for pw in os.listdir("./pathways"):
        pathway = parse_pathway("./pathways/" + pw)
        pathways.append(pathway)
    pathways.sort()

    mutations_data = pandas.read_csv("TRIBE2_seq_res.csv")
    patients_log = pandas.read_csv("TRIBE2_db.csv")
    columns = ["PatientFirstName"] + [pw[0] for pw in pathways]

    print(f"Loaded {len(patients_log)} patients and {len(mutations_data)} mutations")
    
    out = pandas.ExcelWriter("TRIBE2_avgs.xlsx", engine="openpyxl")

    final = process_patients(patients_log[patients_log["arm"] == 0]["PatientFirstName"], mutations_data, pathways)
    print("ARM0 processed")
    if not os.path.exists("arm0"):
        os.mkdir("arm0")
    for pw in pathways:
        name = pw[0]
        sns_plot = sns.displot(final[name])
        sns_plot.savefig(f"arm0/{name}.png")
    arm0d = final.describe()
    arm0d.to_excel(out, "Summary (arm 0)")
    print(arm0d.to_latex())
    df = final.join(patients_log.set_index("PatientFirstName"), on="PatientFirstName")
    df.to_excel(out, "Average mutations (arm 0)", index=False)

    final = process_patients(patients_log[patients_log["arm"] == 1]["PatientFirstName"], mutations_data, pathways)
    print("ARM1 processed")
    if not os.path.exists("arm1"):
        os.mkdir("arm1")
    for pw in pathways:
        name = pw[0]
        sns_plot = sns.displot(final[name])
        sns_plot.savefig(f"arm1/{name}.png")
    final.describe().to_excel(out, "Summary (arm 1)")
    print(final.describe().to_latex())
    df = final.join(patients_log.set_index("PatientFirstName"), on="PatientFirstName")
    df.to_excel(out, "Average mutations (arm 1)", index=False)

    print("Saving results to TRIBE2_avgs.xlsx")
    out.save()
