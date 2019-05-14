# abpred

## Introduction

This repository contains all of the code required to run the [protein-sol abpred](https://protein-sol.manchester.ac.uk) machine learning models trained on the [Jain 2017 dataset](https://doi.org/10.1073/pnas.1616408114) of 137 clinical stage therapeutic antibodies sequences and performance on 12 biophysical characterisation platforms.

The paper describing this work has been submitted for consideration and the preprint is available [here](https://www.biorxiv.org/content/10.1101/625830v1).

This code can run predictions on multiple proteins.

## Installation

## Docker
There is a docker image preconfigured with all necessary software and dependencies to run this software [here](https://cloud.docker.com/repository/docker/maxhebditch/abpred).

## Conda virtual enviroment
For this project, we have supplied a list of yaml configuration file (`abpred_env.yml`) that contains a list of all packages available for installation from [conda](https://conda.io/en/latest/).
If you already have `conda` installed you can run the `install.sh` script that sets up the conda environment and also installs the required [`bestNormalize`](https://cran.r-project.org/web/packages/bestNormalize/vignettes/bestNormalize.html) package for the mathematical transformations.

## Manual
The software is designed to work with the following languages and versions but may work with other versions.

`bash`, `perl 5`, `R 3.5.1`, `python3.6`

For these languages we also require the following `R` packages

`randomforest`, `caret`, `glmnet`, `kernlab`, `bestNormalize`

and `pandas` and its dependencies for python.

# Running code
The `run_preds.sh` shell script is a simple pipeline for running the required code. It requires two arguments to run, the fasta file and an optional identification name. For example

`run_preds.sh test.fasta job1`

# Outputs
In your current directory, the script will return a file for each of the biophysical characterisation platforms as well as a file containing sequence features for each of the input Fv.
