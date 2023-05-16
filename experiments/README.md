# Manuscript Experiments

This sub-repository contains all scripts to used to generate the results 
presented in [Koehler et al.](). Please read the instructions, especially
on how to [download the required data](#data).

## General Notes

To run the experiments you need to install the `pymantra` package.
All packages required to reproduce the experiments are also dependencies of
`pymantra`.
Nevertheless, all packages including the versions used to produce the results
shown in the paper are show in pipenv.txt (coming from pip freeze).

For generating the metabolic networks, we recommend using the
[mantraAPI](https://github.com/lipitum/pymantraAPI) locally (see instructions
on how to start the API as a docker application)

Tow major differences to the package are that the experiments require
**Python >= 3.9** and **matplotlib >= 3.5**.

Python files in this folder follow a scheme:

* generate_*.py generates `mantra`-compatible metabolic networks
* pre_process_*.py handles all pre-processing (normalization,
  transformation etc.)
* run_*.py runs the actual analyses incl. enrichment and produces the
  figures shown in the results section of the paper

All other python files contain utility functions.

## Reproducing Paper Results

The parameters used in the paper are stored in `run_experiments.sh`.
While this is a bash file, the commands are generally the same on Windows as
well.

We recommend running them in a clean environment using Python 3.10 
(originally run with 3.10.4). There might be some errors with numpy >= 1.23.0.

The data required to produce the results in the paper is not contained in the
sub-folders in this repository. Please refer to the next section for 
instructions on how to [download the data](#data).

## Data

Both datasets used in the experiments are publicly available.

The data for the pure metabolome analysis is coming from 
[Xiao et al.](https://www.nature.com/articles/s41422-022-00614-0) 
(Supplementary Table). The downloaded should have be named "FUSCCTNBC.xlsx"
and be placed in the "Metabolome" folder.

The data for the metabolome-microbiome analysis is coming from 
[Franzosa et al.](https://www.nature.com/articles/s41564-018-0306-4) 
(Supplementary Datasets 1, 2 and 4).
The files should be named "41564_2018_306_MOESM3_ESM.xlsx", 
"41564_2018_306_MOESM4_ESM.xlsx" and "41564_2018_306_MOESM6_ESM.xlsx" and be
placed in the "MetabolomeMicrobiome" folder.
