#!/usr/bin/env python3

import sys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

with open(snakemake.output.out, "wt") as w_file :
    with open(snakemake.input.all_out[0], "rt") as r_file :
        for line in r_file :
            if line.startswith("*") :
                break
            else :
                w_file.write(line)
    
    Replicon_nohit = "--- Replicons with no hits: ---\n\n"

    for file_out in snakemake.input.all_out :
        nohit_begin = True
        hit_begin = False

        with open(file_out, "rt") as r_file :
            for line in r_file :
                if line.startswith("*") and hit_begin == False :
                    hit_begin = True
                    w_file.write(line)
                elif "Replicons" in line :
                    hit_begin = False
                    nohit_begin = False
                elif hit_begin :
                    w_file.write(line)

            if nohit_begin :
                Replicon_nohit += file_out.split("/")[-2]+"\n"

    w_file.write(Replicon_nohit)