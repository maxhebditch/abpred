#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
FASTA_in=$1
TimeStamp=$2

cp $FASTA_in reformat.in
perl $DIR/fasta_seq_reformat_export.pl > run_log

cp $FASTA_in $FASTA_in"_ORIGINAL"
mv reformat_out $FASTA_in 

cp $FASTA_in "composition.in";

perl $DIR/seq_compositions_perc_pipeline_export.pl >> run_log

mv composition_all.out seq_composition.txt

rm reformat.in composition.in blah.txt $FASTA_in"_ORIGINAL"

mv run_log run.log

mv seq_composition.txt $TimeStamp"_seq_composition.txt"
