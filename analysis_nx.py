import networkx as nx
import numpy as np
import pandas

def calculate_patient_mutations_with_f(pid, seq_data, pathways, f, factor_famcom=False):
    patient_data = seq_data[
        (seq_data['PatientFirstName'] == pid)
        & (seq_data['Technology'] == 'NGS Q3')
        # We only care about variants and pathogenic mutations
        & (seq_data['TestResult'].isin(['variantdetected', 'Mutated, Pathogenic']))
    ]

    # Realistically, should never happen
    if(patient_data.empty):
        return {}

    results = {}
    for pw in pathways:
        pw_names = nx.get_node_attributes(pw[1], 'label')
        pw_famcom = nx.get_node_attributes(pw[1], 'famcomw')
        pathway_mutations = patient_data[patient_data['Biomarker'].isin(pw_names.values())]
        if(pathway_mutations.empty):
           results[pw[0]] = np.float64(0.0)
           continue
        
        weights = f(pw[1])
        # We have them stored by ID, need them by biomarker name
        weights = {pw_names[k]: weights[k] * (pw_famcom[k] if factor_famcom else 1) for k, v in pw_names.items()}
        weights = pandas.Series(weights)
        patient_mutations = pathway_mutations.groupby(['Biomarker']).max()['NGS_PercentMutated']
        perc_mutation = weights.mul(patient_mutations, fill_value=np.float64(0.0)).sum() / weights.sum()
        results[pw[0]] = perc_mutation

    return results

def process_patients_with_f(patients, f, pathways, mutations_data, complexes=False):
    results = {}
    for patient in patients:
        results[patient] = calculate_patient_mutations_with_f(patient, mutations_data, pathways, f, complexes)

    return pandas.DataFrame.from_dict(results, orient='index').rename_axis('PatientFirstName').reset_index()
