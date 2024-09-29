import os
import sys
import argparse
from functions.hmmer import *

'''
Bikash Kumar Bhandari, Nick Goldman, 
A generalized protein identification method for novel and diverse sequencing 
technologies, NAR Genomics and Bioinformatics, Volume 6, Issue 3, September 2024,
lqae126, https://doi.org/10.1093/nargab/lqae126
'''

def check_arg(args=None):
    '''arguments.
    '''
    parser = argparse.ArgumentParser(prog='Program',
                     description='Protein identification using novel sequencing devices',
                     epilog='bkb3')
    parser.add_argument('-v', '--version',
                    action='version',
                    version='%(prog)s ' + '1',
                    help="Show program's version number and exit.")
    parser.add_argument('-r', '--readings',
                    type=str,
                    help='Input decoded readings file',
                    required=True)

    results = parser.parse_args(args)
    return results.readings

def main():
    '''Main func
    '''
    file_name = os.path.basename(r)
    readings = pd.read_csv(r).values

    hmm = hmm_build(readings, file_name, file_name)

    results = score(hmm, sequences, background)
    results.to_csv(f'results/{file_name}_results.csv', \
                                            index=None)
    print(results)

if __name__ == '__main__':
    r = check_arg(sys.argv[1:])
    main()

