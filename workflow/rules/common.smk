##########################################################################

import os, sys
from snakemake.utils import validate
import pandas as pd

##########################################################################


def get_final_output():
    """
    Generate final output name
    """
    final_output = [
        os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.report"
        ),
        os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.out"
        ),
        os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.summary"
        ),
        os.path.join(
            OUTPUT_FOLDER, 
            "macsyfinder.tab"
        ),
    ]

    return final_output


##########################################################################

# GLOBAL VARIABLES

OUTPUT_FOLDER = config["output_folder"]

##########################################################################

FASTA_FOLDER = config["fasta_folder"]

if not config["extension"].startswith("."):
    EXT_FILE = config["extension"]
else:
    EXT_FILE = config["extension"][1:]

# (REPLICON_NAME,) = glob_wildcards(os.path.join(FASTA_FOLDER, f"{{replicon_name}}.{EXT_FILE}"))

##########################################################################

if config["models_tsv"]:
    dict_models = pd.read_table(config["models_tsv"], index_col=0)
    dict_models = dict_models.iloc[:,0].to_dict()
else: 
    dict_models = {}

REPLICON_NAME = [i for i in dict_models.keys() if os.path.isfile(os.path.join(FASTA_FOLDER, f"{i}.{EXT_FILE}"))] 

