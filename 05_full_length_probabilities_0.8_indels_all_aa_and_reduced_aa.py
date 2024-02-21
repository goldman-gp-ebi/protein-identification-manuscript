import os
import shutil
from itertools import product
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("paper", font_scale=2.0)
sns.set_style('whitegrid')

from functions import *

from pandarallel import pandarallel

pandarallel.initialize(progress_bar=False)

with SequenceFile("data/uniprot-9606.fasta", digital=True, alphabet=alphabet) as seq_file:
    sequences = list(seq_file)
    

uniprot = fasta_reader('data/uniprot-9606.fasta')



aa_knowns = ['LSEAG', 'CKY', 'WMCHY', 'ACDEFGHIKLMNPQRSTVWY']

prob_range = [0.8,]


all_indels = [i for i in product([0, 10, 20, 30, 40, 50, 60], [0, 10, 20, 30, 40, 50, 60])]


repeats = 1


random = 0

for a, b in enumerate(all_indels):
    ins_rate = b[0]
    del_rate = b[1]

    for aa_known in aa_knowns:
        for rep in range(repeats):

            random = 12345
            for prob in prob_range:
#                 fname = f'{rep}_full_length_{aa_known}_max_prob_{prob}_ins_{ins_rate}_del_{del_rate}_1_ins.pkl.gz'
#                 if not os.path.isfile('results/scan_results/'+fname):
                temp_df = uniprot.copy()
                temp_df['temp'] = temp_df.values.tolist() 
                temp_df['posteriors'] =  temp_df.temp.parallel_apply(\
                                     lambda x: gen_reads_indels(seq=x[1], aa_known=aa_known, max_prob=prob, \
                                                               ins_rate=ins_rate, del_rate=del_rate, random=random, ))


                temp_df['temp1'] = temp_df[['Accession', 'posteriors']].values.tolist() 
                temp_df['hmm'] = temp_df.temp1.parallel_apply(\
                         lambda x: hmm_build(x[1], x[0], f'{rep}_full_length_{x[0]}_{aa_known}_max_prob_{prob}_ins_{ins_rate}_del_{del_rate}'))


                hmms = np.array_split(temp_df.hmm, 10)
                hmms = [i.tolist() for i in hmms]

                results = []
                for i, v in enumerate(hmms):
                    results.append(score(v, sequences, background))
                #             print(f'\t\t\t\tdone: {i}', end='\r')
                fname = f'{rep}_full_length_{aa_known}_max_prob_{prob}_ins_{ins_rate}_del_{del_rate}_1_ins.pkl.gz'

                res_df = pd.concat(results)
                res_df.to_pickle('results/scan_results/' + fname)
                del temp_df, res_df, hmms
                shutil.rmtree('temp/', ignore_errors=True)
                try:
                    os.makedirs('temp/')
                except FileExistsError:
                    pass