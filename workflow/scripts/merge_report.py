#!/usr/bin/env python3

import sys
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

columns = [
    "Hit_Id",
    "Replicon_name",
    "Position",
    "Sequence_length",
    "Gene",
    "Reference_system",
    "Predicted_system",
    "System_Id",
    "System_status",
    "Gene_status",
    "i-evalue",
    "Score",
    "Profile_coverage",
    "Sequence_coverage",
    "Begin_match",
    "End_match"
    ]

header = "\t".join(columns)

with open(snakemake.output.out, "wt") as w_file :
    w_file.write(header+"\n")

    for file_out in input.all_out : 
        pd.read_table(
            file_out, 
            names=columns, 
            comment="#", 
            skiprows=1).to_csv(w_file, index=False, header=False, sep="\t")