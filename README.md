This repository contains all the Jupyter notebooks and scripts to reproduce the results of the paper [Paper name](https://doi.org/). 

If you wish to use our method in your protein identification experiments, the [dist](dist/) directory contains a cleaned up version of the necessary files, a program implementation of our method, sample data and instructions to get you started.  

## Environment
python 

```
Python 3.9.7 (default, Sep 16 2021, 13:09:58) 
[GCC 7.5.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.

```

python dependencies

```
pandas 1.3.4
seaborn 0.12.2
matplotlib 3.4.3
numpy 1.22.3
pyhmmer 0.6.3
pandarallel 1.6.1
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

## Running the jupyter notebooks (*.ipynb)
The code was written in Python v3.9.7.
The notebooks depend upon the data generated from other notebooks and scripts for eg. to generate figures. Hence, they are ordered using numeric prefix in their order of execution. 

1. 00_database_statistics.ipynb
2. Please run other ```.py``` files now. They share the same ```temp``` directory so please run one file at a time to avoid conflicts. The scripts will generate results to be used by the following notebooks. 
3. 01-data-analysis.ipynb
4. 02_plots.ipynb
5. 03_combined_result_from_10_fragments.ipynb

It is recommended to run these files in a __HPC environment__ with sufficient access to disk space, memory (200 - 300 GiB) and cores (~50). While the protein identification for a single sequence is fast, many of the scripts will attempt identification of each sequences in the database (N=20,181) for different combinations of parameters. Thus, some of the resulting files will be quite big and the process will take a long time. The scripts will also create several directories for temp files. There will be many temp files in those directores, but are cleared once the execution completes. This step will also take some time.

### Funding
[EU Horizon 2020 grant agreement no. 964363](https://cordis.europa.eu/project/id/964363)


### Citation
 - 
