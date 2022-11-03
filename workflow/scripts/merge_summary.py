#!/usr/bin/env python3

import sys
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

columns = [
    "Replicon_name",
    "System_Id",
    "Reference_system",
    "System_status",
    "Nb_loci",
    "Nb_Ref_mandatory",
    "Nb_Ref_accessory",
    "Nb_Ref_Genes_detected_NR",
    "Nb_Genes_with_match",
    "System_length",
    "Nb_Mandatory_NR",
    "Nb_Accessory_NR",
    "Nb_missing_mandatory",
    "Nb_missing_accessory",
    "List_missing_mandatory",
    "List_missing_accessory",
    "Loci_positions",
    "Occur_Mandatory",
    "Occur_Accessory",
    "Occur_Forbidden"
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