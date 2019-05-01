#!/bin/bash

curDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

conda env create --name abpred --file abpred_env.yml
conda activate $curDir/abpred
echo "install.packages(\"bestNormalize\", repos=\"https://www.stats.bris.ac.uk/R/\")" | ./env/abpred/bin/R --no-save
