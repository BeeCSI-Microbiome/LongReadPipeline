######## Example run #################################

# snakemake --core 2 --use-conda --rerun-incomplete --conda-frontend conda --use-singularity
# this is used for if you don't have conda.

#######################################################

# Import Libraries - check which are necessary and whgich are not. 
import os
import sys 
import csv
import pandas as pd
import numpy as np

from snakemake.io import load_configfile
from snakemake.utils import update_config as snakemake_update_config

configfile: "./config.yaml"

sampleTable = pd.read_csv("./samples.csv", header = 0)
SAMPLES = list(sampleTable['name'].to_numpy())

##########################################################################
# LongRead - Start
##########################################################################

# Snakefile includes
include: "rules/longReadAssembly.smk"
include: "rules/longReadBinning.smk"

localrules: # may be better to do this based on the various sequencing technologies. 
    all,
    binning,


rule binning:
    input:
        "finished_binning",
    output:
        touch("binning_complete"),


rule all: 
    input:
        "finished_pipeline",
    output:
        touch("pipeline_complete")

