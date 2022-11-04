#!/usr/bin/env python3

import sys
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

for file_tab in input.all_out :
    if os.stat(file_tab).st_size != 0 :
        with open(file_tab, "rt") as r_file :
            header = r_file.readline().rstrip().replace("#", "")
            break

with open(str(output), "wt") as w_file :
    w_file.write(header+"\n")

    for file_out in snakemake.input.all_out : 
        pd.read_table(
            file_out, 
            names=header.split(), 
            comment="#", 
            skiprows=1).to_csv(w_file, index=False, header=False, sep="\t")