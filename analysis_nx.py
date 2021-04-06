import networkx as nx
import numpy as np
import pandas
import logging as log


def calculate_patient_mutations_with_f(pid, seq_data, pathways, f, factor_famcom=False):
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
        pathway_mutations = patient_data[patient_data["Biomarker"].isin(pw.get_genes())]
        if pathway_mutations.empty:
            results[pw.name] = np.float64(0.0)
            continue

        weights = pw.calculate_measure(f, factor_famcom)
        patient_mutations = pathway_mutations.groupby("Biomarker").max()[
            "NGS_PercentMutated"
        ]
        total_weights = weights.sum()
        if total_weights != 0:
            perc_mutation = (
                weights.mul(patient_mutations, fill_value=np.float64(0.0)).sum()
                / total_weights
            )
        else:
            perc_mutation = 0
        results[pw.name] = perc_mutation

    return results


def process_patients_with_f(patients, f, pathways, mutations_data, complexes=False):
    results = {}
    for patient in patients:
        results[patient] = calculate_patient_mutations_with_f(
            patient, mutations_data, pathways, f, complexes
        )

    return (
        pandas.DataFrame.from_dict(results, orient="index")
        .rename_axis("PatientFirstName")
        .reset_index()
    )
