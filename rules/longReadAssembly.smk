import pandas as pd
import numpy as np
from snakemake.io import load_configfile
from snakemake.utils import update_config as snakemake_update_config

sample = SAMPLES
fastqfolder = ""
configfile:"./config.yaml"


######## Improvements ################
# 1. make a results directory in config and greatly truncate the scripts. 
# 2. Make a sample directory and "simple" sample loading system

#######################################



rule initializeLongRead:
    input: # Fix this input so that it works for a specific sample file. 
        #fastq = ""
        "samples/{sample}.fastq"
        #unpack(sampleTable['name' == {sample}, 'location']),
    output:
        "results/{sample}/trimmed_fasta/{sample}_trimmed.fastq.gz"
    priority: 80
    params:
        # extras = config.get("porechop_extras")
    container:
        "./rules/envs/porechop:0.2.4--py39hc16433a_3"  
    log:
        "results/logs/Assembly/{sample}_initialization.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell: 
        "porechop -i samples/{wildcards.sample}.fastq -o results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz -t 20 2> {log}"  


rule longReadAssembly:
    input:
        "results/{sample}/trimmed_fasta/{sample}_trimmed.fastq.gz"
    output:
        "results/{sample}/flye_assembly/assembly.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/flye:2.9--py39h6935b12_1"
    log:
        "results/logs/Assembly/{sample}_assembly.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "flye --nano-raw results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz --meta --genome-size 4m --out-dir results/{wildcards.sample}/flye_assembly -i 0 -t 20 2> {log} "

 

rule longReadPolishing1_minimap:
    input:
       "results/{sample}/flye_assembly/assembly.fasta",
        # fastq = "samples/{sample}.fastq"
    output:
        "results/{sample}/flye_assembly/{sample}.sam"
    priority: 80
    params:
    container:
        "./rules/envs/minimap2:2.9--1"
    log:
        "results/logs/Assembly/{sample}_Polishing1_minimap.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "minimap2 -ax map-ont results/{wildcards.sample}/flye_assembly/assembly.fasta results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz > results/{wildcards.sample}/flye_assembly/{wildcards.sample}.sam 2> {log}"



rule longReadPolishing1_racon:
    input:
        "results/{sample}/flye_assembly/{sample}.sam",
    output:
        "results/{sample}/flye_assembly/{sample}_racon1.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/racon:1.5.0--h7ff8a90_0"
    log:
        "results/logs/Assembly/{sample}_Polishing1_racon.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"], 
    shell: 
        "racon -t 20 results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz results/{wildcards.sample}/flye_assembly/{wildcards.sample}.sam results/{wildcards.sample}/flye_assembly/assembly.fasta > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon1.fasta 2>> {log}"


rule longReadPolishing1_sed:
    input:
       "results/{sample}/flye_assembly/{sample}_racon1.fasta",
    output:
        "results/{sample}/flye_assembly/{sample}_racon1.mod.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/samtools:1.9--h91753b0_8"
    log:
        "results/logs/Assembly/{sample}_Polishing1_sed.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],
    shell:
        "sed -n '/^>/,$p' results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon1.fasta | sed 's/\s.*$//g' > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon1.mod.fasta 2>> {log}"

rule longReadPolishing2_minimap:
    input:
       "results/{sample}/flye_assembly/{sample}_racon1.mod.fasta",
    output:
        "results/{sample}/flye_assembly/{sample}_racon1.map.sam"
    priority: 80
    params:
    container:
        "./rules/envs/minimap2:2.9--1"
    log:
        "results/logs/Assembly/{sample}_Polishing2_minimap.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],
    shell: 
        "minimap2 -ax map-ont results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon1.mod.fasta results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon1.map.sam 2>> {log}"

rule longReadPolishing2_racon:
    input:
       "results/{sample}/flye_assembly/{sample}_racon1.map.sam",
    output:
        "results/{sample}/flye_assembly/{sample}_racon2.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/racon:1.5.0--h7ff8a90_0"
    log:
        "results/logs/Assembly/{sample}_Polishing2_racon.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"], 
    shell:
        "racon -t 20 results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon1.map.sam results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon1.mod.fasta > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon2.fasta 2>> {log}"


rule longReadPolishing2_sed:
    input:
       "results/{sample}/flye_assembly/{sample}_racon2.fasta",
    output:
        "results/{sample}/flye_assembly/{sample}_racon2.mod.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/samtools:1.9--h91753b0_8"
    log:
        "results/logs/Assembly/{sample}_Polishing2_sed.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"], 
    shell:
        "sed -n '/^>/,$p' results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon2.fasta | sed 's/\s.*$//g' > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon2.mod.fasta 2>> {log}"

rule longReadPolishing3_minimap:
    input:
        "results/{sample}/flye_assembly/{sample}_racon2.mod.fasta"
    output:
        "results/{sample}/flye_assembly/{sample}_racon2.map.sam"
    priority: 80
    params:
    container:
        "./rules/envs/minimap2:2.9--1"
    log:
        "results/logs/Assembly/{sample}_Polishing3_minimap.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "minimap2 -ax map-ont results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon2.mod.fasta results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon2.map.sam 2> {log} "

rule longReadPolishing3_racon:
    input:
        "results/{sample}/flye_assembly/{sample}_racon2.map.sam"
    output:
        "results/{sample}/flye_assembly/{sample}_racon3.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/racon:1.5.0--h7ff8a90_0"
    log:
        "results/logs/Assembly/{sample}_Polishing3_racon.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "racon -t 20 results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon2.map.sam results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon2.mod.fasta > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon3.fasta 2>> {log} "


rule longReadPolishing3_sed:
    input:
        "results/{sample}/flye_assembly/{sample}_racon3.fasta"
    output:
        "results/{sample}/flye_assembly/{sample}_racon3.mod.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/samtools:1.9--h91753b0_8"
    log:
        "results/logs/Assembly/{sample}_Polishing3_sed.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "sed -n '/^>/,$p' results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon3.fasta | sed 's/\s.*$//g' > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon3.mod.fasta 2>> {log} "


rule longReadPolishing4_minimap:
    input:
        "results/{sample}/flye_assembly/{sample}_racon3.mod.fasta"
    output:
        "results/{sample}/flye_assembly/{sample}_racon3.map.sam"
    priority: 80
    params:
    container:
        "./rules/envs/minimap2:2.9--1"
    log:
        "results/logs/Assembly/{sample}_Polishing4_minimap.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "minimap2 -ax map-ont results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon3.mod.fasta results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon3.map.sam 2>> {log} "


rule longReadPolishing4_racon:
    input:
        "results/{sample}/flye_assembly/{sample}_racon3.map.sam"
    output:
        "results/{sample}/flye_assembly/{sample}_racon4.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/racon:1.5.0--h7ff8a90_0"
    log:
        "results/logs/Assembly/{sample}_Polishing4_racon.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "racon -t 20 results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon3.map.sam results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon3.mod.fasta > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon4.fasta 2>> {log} "



rule longReadPolishing4_sed:
    input:
        "results/{sample}/flye_assembly/{sample}_racon4.fasta"
    output:
        "results/{sample}/flye_assembly/{sample}_racon4.mod.fasta"
    priority: 80
    params:
    container:
        "./rules/envs/samtools:1.9--h91753b0_8"
    log:
        "results/logs/Assembly/{sample}_Polishing4_sed.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "sed -n '/^>/,$p' results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon4.fasta | sed 's/\s.*$//g' > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon4.mod.fasta 2>> {log} "

