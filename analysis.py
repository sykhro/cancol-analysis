from pathways import *
import pandas
import numpy as np
logging.basicConfig(level=logging.DEBUG, filename='pathway_parser.log')

pathways = []
for pw in os.listdir('./pathways'):
    pathway = parse_pathway('./pathways/' + pw)
    pathways.append(pathway)

patients_log = pandas.read_excel('TRIBE2_db.xlsx')
patients_log = patients_log[patients_log['arm'] == 0]
print('Patients list ready')

mutations_data = pandas.read_excel('TRIBE2_seq_res.xlsx')
print('Sequencing results ready')

patients = patients_log['PatientFirstName']

results = {}
columns = ['PatientFirstName'] + [pw[0] for pw in pathways]

print(f'About to process {len(patients)} patients, this may take a while')
for patient in patients:
    results[patient] = {}
    patient_data = mutations_data[
        (mutations_data['PatientFirstName'] == patient)
        & (mutations_data['Technology'] == 'NGS Q3')
        # We only care about variants and pathogenic mutations
        & (mutations_data['TestResult'].isin(['variantdetected', 'Mutated, Pathogenic']))
    ]

    for pw in pathways:
        pw_genes = get_genes(pw[1])
        pathway_mutations = patient_data[patient_data['Biomarker'].isin(pw_genes)]
        if(pathway_mutations.empty):
           results[patient][pw[0]] = np.float64(0.0)
           continue

        perc_mutation = pathway_mutations.groupby(['Biomarker']).max()['NGS_PercentMutated'].sum() / len(pw_genes)
        results[patient][pw[0]] = perc_mutation

final = pandas.DataFrame.from_dict(results, orient='index')
out = pandas.ExcelWriter('TRIBE2_avgs.xlsx', engine='openpyxl')
final.describe().to_excel(out, 'Summary')
final.to_excel(out, 'Average mutations')
out.save()

