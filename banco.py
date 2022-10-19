#!/usr/bin/python
# vim: set fileencoding= utf8 :

import pandas as pd
from os import listdir
from os.path import isfile, join

amp_files = [f for f in listdir("/home/idengene/CNV/GenRef/banco") if isfile(join("/home/idengene/CNV/GenRef/banco", f))]

with open('all.amp.bed', 'w') as outfile:
    for fname in amp_files:
      with open(fname) as infile:  
        outfile.write(infile.read())

amp = pd.read_csv('all.amp.bed', sep ='\t')
amp_series = amp.groupby(['#Chrom', 'Start', 'End', 'SV_type']).size()
amp_df = amp_series.to_frame(name = 'N_Samples').reset_index()
amp_df = amp_df.drop(0)
amp_df = amp_df[["#Chrom", "Start", "End", "N_Samples", "SV_type"]]
amp_df.to_csv('/home/idengene/CNV/GenRef/banco.bed', sep='\t', index=False)

