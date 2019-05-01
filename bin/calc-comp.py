#!/usr/bin/env python
# coding: utf-8

import subprocess
import csv
import os
import pandas as pd
from optparse import OptionParser

usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("--FASTA", type="string", dest="fasta")
parser.add_option("--calcID", type="string", dest="calcID")
(options, args) = parser.parse_args()

py_path = os.path.dirname(os.path.realpath(__file__))

fasta_file = options.fasta
calcID = options.calcID

header = ["antibody","scaled_sol","K - R","D - E","Number of amino acids","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","KpR","DpE","PmN","PpN","F+W+Y","pI","kyte doolittle","absolute charge","Fold propensity","disorder","entropy","beta proprensity"]

run_comp   = py_path+"/calc-comp/run_model.sh"

comp_file  = calcID+"_seq_composition.txt" 
outp_file  = calcID+"_table_composition.csv" 

subprocess.call([run_comp,fasta_file,calcID])

comp_raw = []
names = []

for line in open(comp_file,"r"):
    if line.startswith("WHOLE-SEQ,>"):
        clean = line.replace(" ","")
        clean = clean.replace("\n","")
        comp_line = (clean.split(","))
        names.append(comp_line[1])
        comp_vals = [float(i) for i in comp_line[2:]]
        del comp_vals[3]
        comp_raw.append(comp_vals)
        
with open(outp_file,"w") as outp_table_file:
    writer = csv.writer(outp_table_file)
    writer.writerow(header)
    
    i = 0
    for protein_vals in comp_raw:
        csv_vals = [names[i],"Na"] + protein_vals
        writer.writerow(csv_vals)
        i=i+1
