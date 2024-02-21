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

from functions1 import *

from pandarallel import pandarallel

pandarallel.initialize(progress_bar=False)

with SequenceFile("data/uniprot-9606.fasta", digital=True, alphabet=alphabet) as seq_file:
    sequences = list(seq_file)
    

uniprot = fasta_reader('data/uniprot-9606.fasta')


frag_length = 100
num_frag = 1
aa_knowns = ['LSEAG', 'CKY', 'WMCHY', 'ACDEFGHIKLMNPQRSTVWY']

prob_range = [0.8,]

all_indels = [i for i in product([0, 10, 20, 30, 40, 50, 60], [0, 10, 20, 30, 40, 50, 60])]

def fragment(seq, frag_len, random):
    if len(seq) <= frag_len:
        return seq
    rnd = np.random.default_rng(random)
    idx = rnd.integers(0, len(seq) - frag_len + 1) 
    return seq[idx:(idx+frag_len)]


for a, b in enumerate(all_indels):
    ins_rate = b[0]
    del_rate = b[1]
    for frag in range(num_frag):
        for aa_known in aa_knowns:
            random = ins_rate + del_rate + len(aa_knowns) + frag
            for prob in prob_range:
                fname = f'{frag}_frag_length_100_{aa_known}_max_prob_{prob}_ins_{ins_rate}_del_{del_rate}_1_ins.pkl.gz'
                if not os.path.isfile('results/scan_results/'+fname):
                    temp_df = uniprot.copy()
                    temp_df['Sequence'] = uniprot['Sequence'].apply(lambda x: \
                                            fragment(seq=x, frag_len=frag_length, random=random))
                    temp_df['temp'] = temp_df.values.tolist() 
                    temp_df['posteriors'] =  temp_df.temp.parallel_apply(\
                                         lambda x: gen_reads_indels(seq=x[1], aa_known=aa_known, max_prob=prob, \
                                                                   ins_rate=ins_rate, del_rate=del_rate, random=random, ))


                    temp_df['temp1'] = temp_df[['Accession', 'posteriors']].values.tolist() 
                    temp_df['hmm'] = temp_df.temp1.parallel_apply(\
                             lambda x: hmm_build(x[1], x[0], f'{frag}_frag_length_100_{x[0]}_{aa_known}_max_prob_{prob}_ins_{ins_rate}_del_{del_rate}'))


                    hmms = np.array_split(temp_df.hmm, 10)
                    hmms = [i.tolist() for i in hmms]

                    results = []
                    for i, v in enumerate(hmms):
                        results.append(score(v, sequences, background))
                    #             print(f'\t\t\t\tdone: {i}', end='\r')
                    fname = f'{frag}_frag_length_100_{aa_known}_max_prob_{prob}_ins_{ins_rate}_del_{del_rate}_1_ins.pkl.gz'

                    res_df = pd.concat(results)
                    res_df.to_pickle('results/scan_results/' + fname)
                    del temp_df, res_df, hmms
                    shutil.rmtree('temp1/', ignore_errors=True)
                    try:
                        os.makedirs('temp1/')
                    except FileExistsError:
                        pass