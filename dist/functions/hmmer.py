import os
import fileinput
from subprocess import run, PIPE
import numpy as np
import pandas as pd
import pyhmmer
from pyhmmer.plan7 import HMM, Background, Pipeline, HMMFile
from pyhmmer.easel import SequenceFile


aa = 'ACDEFGHIKLMNPQRSTVWY'


def fasta_reader(file):
    '''Converts .fasta to a pandas dataframe with accession as index
    and sequence in a column 'sequence'
    '''
    fasta_df = pd.read_csv(file, sep='>', lineterminator='>', header=None)
    fasta_df[['Accession', 'Sequence']] = fasta_df[0].str.split('\n', 1, \
                                        expand=True)
    fasta_df['Accession'] = fasta_df['Accession'].str.split('\s').apply(lambda x: x[0])
    fasta_df['Sequence'] = fasta_df['Sequence'].replace('\n', '', regex=True).\
                            astype(str).str.upper().str.replace('U', 'C')
    fasta_df.drop(0, axis=1, inplace=True)
    fasta_df = fasta_df[(fasta_df.Sequence != '') & (fasta_df.Sequence != 'NONE') & \
                        fasta_df['Sequence'].str.isalpha()]
    fasta_df = fasta_df[(~fasta_df.Sequence.str.contains('X')) & (~fasta_df.Sequence.str.contains('Z'))].copy()
    final_df = fasta_df.dropna()

    return final_df


def generate_reads(seq, aa_known, max_prob):
    '''
    Return uniform prob if AA where AA is unknown
    Else return max_prob for AA known, divide rem 
    prob equally    
    '''
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    probs = []
    remaining_prob = np.round((1-max_prob)/19, 2)
    uniform_dist = [1/20]*20
    for i, v in enumerate(seq):
        p = []
        if v not in aa_known:
            p = uniform_dist
        else:
            p = [remaining_prob]*20
            p[aa.index(v)] = max_prob
        probs.append(p)
    return probs


cwd = os.getcwd()

dirs_to_make = ['temp', 'data', 'results', ]


tmp, data, results = [os.path.join(cwd, i) for i in dirs_to_make]


for d in [tmp, data, results,]:
    try:
        os.makedirs(d)
    except FileExistsError:
        pass

alphabet = pyhmmer.easel.Alphabet.amino()
background = Background(alphabet)
alphabet = pyhmmer.easel.Alphabet.amino()

with SequenceFile("data/uniprot-9606.fasta", digital=True, alphabet=alphabet) as seq_file:
    sequences = list(seq_file)

def write_hmm(readings, accession, outfname, transition_prob=None, insert_emissions=None,
             alphabet=pyhmmer.easel.Alphabet.amino()):

    M = len(readings)

#     Since the first row corresponds to the entry probabilities, 
#     the emissions are unused. By convention, it should still 
#     contain valid probabilities, so it will always be set as 
#     follow with 1 probability for the first symbol, and 0 for 
#     the rest.
    
    readings =  np.concatenate([[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0]], 
                         readings],)
    
    hmm = HMM(M, alphabet)
    
    # transition probabilites
    # m->m     m->i     m->d     i->m     i->i     d->m     d->d
    # Start transition probability (Always a match)
    start_tp = [[1, 0.  , 0.  , 0.54, 0.46, 1., 0 ]]
    if transition_prob == None:
        tp = [[0.8, 0.1, 0.1, 0.54, 0.46, 0.77, 0.23]] * M
    else:
        tp = transition_prob

    all_tp = start_tp + tp

    
    for idx, val in enumerate(all_tp):
        for idx1, trans_p in enumerate(val):
            hmm.transition_probabilities[idx, idx1] = trans_p
    
    # Insert emissions
    if insert_emissions == None:
        ins_em =  [background.residue_frequencies]*(M+1)
       
    else:
        ins_em = insert_emissions
    

    for idx, val in enumerate(ins_em):
        for idx1, ins_p in enumerate(val):
            hmm.insert_emissions[idx, idx1] = ins_p


    # Match probabilites
    match_prob = readings
    for idx, val in enumerate(match_prob):
        for idx1, match_p in enumerate(val):
            hmm.match_emissions[idx, idx1] = match_p

    hmm.set_composition()
    hmm.name = accession.encode()
    output_fname = 'temp/{}.hmm'.format(outfname)
    with open(output_fname, 'wb') as output_file:
        hmm.write(output_file)
    return output_fname


def calibrate_hmm(hmmfile):
    with fileinput.FileInput(hmmfile, inplace=True, backup=None) as file:
        for line in file:
            if 'CONS  no' in line:
                print(line.replace('CONS  no', 'CONS  yes'), end='')
            else:
                print(line, end='')

    vit_args = ['hmmsim', '--vit', '--fast', '--seed', '12345', hmmfile]
    vit_res = run(vit_args, stdout=PIPE, stderr=PIPE,)

    msv_args = ['hmmsim', '--msv', '--fast', '--seed', '12345', hmmfile]
    msv_res = run(msv_args, stdout=PIPE, stderr=PIPE,)

    fwd_args = ['hmmsim', '--fwd', '--fast', '--seed', '12345', hmmfile]
    fwd_res = run(fwd_args, stdout=PIPE, stderr=PIPE,)

    calibration_names = ['STATS LOCAL MSV', 
                         'STATS LOCAL VITERBI',
                         'STATS LOCAL FORWARD'
                        ]

    calibration_strings = []


    for i, v in enumerate([msv_res, vit_res, fwd_res]):
        res = [i for i in v.stdout.decode('ascii').split('\n#')[0].split(' ') if i != '']
        calibration_strings.append('{} {}  {}'.format(calibration_names[i], 
                                                  res[-3], res[-2]
                                                  ))

    cal = '\n'.join(calibration_strings)
    with fileinput.FileInput(hmmfile, inplace=True, backup=None) as file:
        for line in file:
            if 'MAP   no' in line:
                print(line.replace('MAP   no\n', f"MAP   no\n{cal}\n"), end='')
            else:
                print(line, end='')


def hmm_build(readings, accession, outfname, transition_prob=None, insert_emissions=None,
             alphabet=pyhmmer.easel.Alphabet.amino()):
    output = write_hmm(readings, accession, outfname, transition_prob=transition_prob,
              insert_emissions=insert_emissions, alphabet=pyhmmer.easel.Alphabet.amino())
    
    calibrate_hmm(output)
    
    with HMMFile(output) as hmm_file:
        hmm = next(hmm_file)
    return hmm


def score(hmms, sequences: "easel.SequenceFile", background=background):
#     pipeline = Pipeline(alphabet, background=background)
#     hits = pipeline.search_hmm(query=hmm, sequences=sequences)
    done_seqs = 0
    def progress(hmm, total):
        nonlocal done_seqs
        done_seqs+= 1
        print(f'{round(done_seqs * 100/total)} %', end='\r')
    
    if type(hmms) is not list:
        hmms = [hmms]

    temp_hits = pyhmmer.hmmer.hmmsearch(hmms, sequences, callback=progress,\
                                       bias_filter=False)
    all_hits = [next(temp_hits) for _ in range(len(hmms))]

    temp_list = []
    for h in all_hits:
        temp_list.extend([(h.query_name.decode(), i.name.decode(), i.evalue, i.score) for i in h])

    results = pd.DataFrame(temp_list, columns=['Query', 'Hit', 'E-value', 'Score'])
    return results

def gen_stats(hmms, result_df):
    tt = result_df.groupby(['Query', 'Hit'], \
                  sort=False,)[['E-value', 'Score', ]].max()

    tt = tt.loc[tt.groupby(['Query'], sort=False)['E-value'].idxmin()].reset_index()
    total_queries = len(hmms)
    total_hits = tt.shape[0]
    total_identified = tt[tt.Query == tt.Accession].shape[0]
    
    return {'Total queries': total_queries, 'Hits': total_hits, 'Identified': total_identified}
