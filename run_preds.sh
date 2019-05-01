#!/bin/bash

curDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

$curDir/env/abpred/bin/python $curDir/bin/calc-comp.py --FASTA $1 --calcID $2
$curDir/env/abpred/bin/Rscript $curDir/bin/scale-and-pred.r $2 $curDir
rm run.log
