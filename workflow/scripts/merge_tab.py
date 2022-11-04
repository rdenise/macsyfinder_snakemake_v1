#!/usr/bin/env python3

import sys, os
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

list_tab = []

for file_out in snakemake.input.all_out : 
    if os.stat(file_out).st_size != 0 :
        with open(file_out, "rt") as r_file :
            columns = r_file.readline().rstrip().replace("#", "").split()

        dtypes = {columns[0]:"string"}
        dtypes.update({column:int for column in columns[1:]})

        list_tab.append(
            pd.read_table(
                file_out, 
                names=columns, 
                comment="#",
                dtype=dtypes,
                skiprows=1
            )
        )

pd.concat(list_tab).fillna(0).to_csv(snakemake.output.out, index=False, sep="\t")