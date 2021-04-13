import networkx as nx
import numpy as np
import pandas
import logging as log
import analysis as lan
from collections import namedtuple

PathwayConfig = namedtuple("PathwayConfig", ["measure", "hierarchy"])


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


def calculate_patient_mutations_with_config(
    pid, seq_data, pathways, legacy_pathways, config
):
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
    legacy_to_compute = []
    for pw in pathways:
        pathway_mutations = patient_data[patient_data["Biomarker"].isin(pw.get_genes())]
        if pathway_mutations.empty:
            results[pw.name] = np.float64(0.0)
            continue

        if config[pw.name].measure == "baseline":
            legacy_to_compute.append(pw.name)
            continue

        weights = pw.calculate_measure(
            config[pw.name].measure, config[pw.name].hierarchy
        )
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

    legacy_results = lan.calculate_patient_mutations(
        pid, seq_data, [pw for pw in legacy_pathways if pw[0] in legacy_to_compute]
    )
    results = {**results, **legacy_results}
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


def process_patients_with_config(
    patients, pathways, legacy_pathways, mutations_data, config
):
    results = {}
    for patient in patients:
        results[patient] = calculate_patient_mutations_with_config(
            patient, mutations_data, pathways, legacy_pathways, config
        )

    return (
        pandas.DataFrame.from_dict(results, orient="index")
        .rename_axis("PatientFirstName")
        .reset_index()
    )
