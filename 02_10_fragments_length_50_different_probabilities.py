import shutil
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

aa_known = 'ACDEFGHIKLMNPQRSTVWY'
frag_length = 50
num_frag = 10

prob_range = [0.9, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05]


def fragment(seq, frag_len, random):
    if len(seq) <= frag_len:
        return seq
    rnd = np.random.default_rng(random)
    idx = rnd.integers(0, len(seq) - frag_len + 1) 
    return seq[idx:(idx+frag_len)]


# generate 10 random fragments for each seq
for frag in range(num_frag):
    for prob in prob_range:
        temp_df = uniprot.copy()
        temp_df['Sequence'] = uniprot['Sequence'].apply(lambda x: \
                                fragment(seq=x, frag_len=frag_length, random=frag))

        temp_df['temp'] = temp_df.values.tolist() 
        temp_df['posteriors'] =  temp_df.temp.parallel_apply(\
                             lambda x: generate_reads(seq=x[1], aa_known=aa_known, max_prob=prob))


        temp_df['temp1'] = temp_df[['Accession', 'posteriors']].values.tolist() 
        temp_df['hmm'] = temp_df.temp1.parallel_apply(\
                 lambda x: hmm_build(x[1], x[0], f'{frag}_frag_length_{frag_length}_{x[0]}_{aa_known}_max_prob_{prob}'))


        hmms = np.array_split(temp_df.hmm, 10)
        hmms = [i.tolist() for i in hmms]

        results = []
        for i, v in enumerate(hmms):
            results.append(score(v, sequences, background))
        #             print(f'\t\t\t\tdone: {i}', end='\r')
        fname = f'{frag}_frag_length_{frag_length}_{aa_known}_max_prob_{prob}_del.pkl.gz'

        res_df = pd.concat(results)
        res_df.to_pickle('results/scan_results/' + fname)
        del temp_df, res_df, hmms
        shutil.rmtree('temp/', ignore_errors=True)
        try:
            os.makedirs('temp/')
        except FileExistsError:
            pass

#         ! rsync -a --delete empty/ temp/
