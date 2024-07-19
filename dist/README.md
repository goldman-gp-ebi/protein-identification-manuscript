This directory contains the scripts suitable for protein sequencing experiments which will identify proteins from a database  `data/uniprot-9606.fasta`. Given a decoded readout from a sequencing device (eg. `sample-reads-P30419.csv`), `main.py` will construct a HMM and scan the database to identify the protein corresponding to the readouts. Further instructions on setting the environment and running the program are given below. 

## Usage
#### Environment
python 

```
Python 3.9.7 (default, Sep 16 2021, 13:09:58) 
[GCC 7.5.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.

```

python dependencies

```
pandas 1.3.4
numpy 1.22.3
pyhmmer 0.6.3
```
[HMMER](http://hmmer.org)

```
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# Copyright (C) 2020 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: hmmsearch [options] <hmmfile> <seqdb>

Basic options:
  -h : show brief help on version and usage 
```

#### Input files 
The input files are the decoded readouts from the sequencing device i.e.. they are the posterior probabilites of the signal after decoding. Therfore, this will be of size read length x 20. Few sample reads are provided as csv. For eg, this reading has length 434, therefore is a 434 x 20 matrix. 

|     |    0 |    1 |    2 |    3 |    4 |    5 |    6 |    7 |    8 |    9 |   10 |   11 |   12 |   13 |   14 |   15 |   16 |   17 |   18 |   19 |
|----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
|  0  | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
|  1  | 0.80 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
|  2  | 0.01 | 0.01 | 0.01 | 0.80 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
|  3  | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.80 | 0.01 | 0.01 | 0.01 | 0.01 |
|  4  | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |
| ... | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  | ...  |
| 431 | 0.80 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| 432 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.80 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| 433 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 | 0.05 |

#### Usage
```python3 main.py -r <your_readings_file>```

This will print and save the results of the found proteins inside the `results` directory as a csv file. For eg: 
```
user@host:~$ python3 main.py -r sample-reads-P30419.csv
                     Query     Hit       E-value       Score
0  sample-reads-P30419.csv  P30419  1.638253e-48  166.687073
1  sample-reads-P30419.csv  O60551  2.344075e-19   70.342026
2  sample-reads-P30419.csv  Q9UF47  3.997793e+00    6.790119

```
The sequence with highest score is the likely protein for the readouts. 

#### Extension to other databases
The current database is ~20K human proteins as a fasta file inside `data` directory. If you prefer a different database, please edit the database location inside `functions/hmmer.py` file to the location of your database. 

If your readouts contains low amount of errors, you can also adjust the transition probabilites inside the `functions/hmmer.py` file.

### Funding
[EU Horizon 2020 grant agreement no. 964363](https://cordis.europa.eu/project/id/964363)

#### Citation
If you use this tool, please cite: 

- 

