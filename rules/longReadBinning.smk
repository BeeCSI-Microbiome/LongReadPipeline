import pandas as pd
import numpy as np


##################  ASSEMBLY    ###################
rule longReadMedakaPolishing:
    input:
        "results/{sample}/flye_assembly/{sample}_racon4.mod.fasta"
    output:
        "results/{sample}/flye_assembly/{sample}_Medaka_polish4/consensus.fasta",
    priority: 80
    params:
    container:
        "./rules/envs/medaka:1.6.0--py38h84d2cc8_0"
    log:
        "results/{sample}/logs/longRead/MedakaPolishing.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "medaka_consensus -i results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz -d results/{wildcards.sample}/flye_assembly/{wildcards.sample}_racon4.mod.fasta -t 25 -m r941_min_high_g303 -o results/{wildcards.sample}/flye_assembly/{wildcards.sample}_Medaka_polish4 2> {log} "
 
 ########################################

rule longReadGenerateDepthFile_minimap:
    input:
        "results/{sample}/flye_assembly/{sample}_Medaka_polish4/consensus.fasta"
    output:
        "results/{sample}/flye_assembly/{sample}_polished_seqs.sam",
    priority: 80
    params:
    container:
        "./rules/envs/minimap2:2.9--1"
    log:
        "results/{sample}/logs/longRead/DepthFile_minimap.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "minimap2 -ax map-ont results/{wildcards.sample}/flye_assembly/{wildcards.sample}_Medaka_polish4/consensus.fasta results/{wildcards.sample}/trimmed_fasta/{wildcards.sample}_trimmed.fastq.gz > results/{wildcards.sample}/flye_assembly/{wildcards.sample}_polished_seqs.sam 2> {log} "



rule longReadGenerateDepthFile_samtools:
    input:
        "results/{sample}/flye_assembly/{sample}_polished_seqs.sam",
    output:
        "results/{sample}/flye_assembly/{sample}_polished_seqs_sort.bam",
    priority: 80
    params:
    container:
        "./rules/envs/samtools:1.9--h91753b0_8"  
    log:
        "results/{sample}/logs/longRead/DepthFile_samtools.log"
    threads: 8
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "samtools view -b results/{wildcards.sample}/flye_assembly/{wildcards.sample}_polished_seqs.sam -o results/{wildcards.sample}/flye_assembly/{wildcards.sample}_polished_seqs.bam 2>> {log} "
        ";"
        "samtools sort -o results/{wildcards.sample}/flye_assembly/{wildcards.sample}_polished_seqs_sort.bam results/{wildcards.sample}/flye_assembly/{wildcards.sample}_polished_seqs.bam 2>> {log} "
        ";"
        "samtools index results/{wildcards.sample}/flye_assembly/{wildcards.sample}_polished_seqs_sort.bam 2>> {log} "


rule longReadMetabatDepthFile:
    input:
        "results/{sample}/flye_assembly/{sample}_polished_seqs_sort.bam",
    output:
        "results/{sample}/{sample}_depth.txt"
    priority: 80
    params:
    container:
        "./rules/envs/metabat2:2.15--h986a166_1"  
    log:
        "results/{sample}/logs/longRead/{sample}_MetabatDepthFile.log"
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth results/{wildcards.sample}/{wildcards.sample}_depth.txt results/{wildcards.sample}/flye_assembly/*sort.bam 2>> {log} "  

rule metabatBinning:
    input:
        "results/{sample}/{sample}_depth.txt",
    output:
        "results/{sample}/{sample}_binning_complete"
    priority: 80
    params:
    container:
        "./rules/envs/metabat2:2.15--h986a166_1"  
    log:
        "results/{sample}/logs/longRead/{sample}_MetabatBinning.log"
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "metabat -i results/{wildcards.sample}/flye_assembly/{wildcards.sample}_Medaka_polish4/consensus.fasta -a results/{wildcards.sample}/{wildcards.sample}_depth.txt -o results/{wildcards.sample}/Binning/{wildcards.sample} -v"
        " ; "
        "metabat -i results/{wildcards.sample}/flye_assembly/{wildcards.sample}_Medaka_polish4/consensus.fasta -o results/{wildcards.sample}/Binning/{wildcards.sample} -v"
        " ; "
        "touch results/{wildcards.sample}/{wildcards.sample}_binning_complete"

rule binning_complete: 
    input:
        expand("results/{sample}/{sample}_binning_complete", sample = SAMPLES),
    output:
        "finished_binning",
        "finished_pipeline",
    priority: 80
    params:
    container:
        "./rules/envs/metabat2:2.15--h986a166_1"  
    log:
        "results/logs/longRead/binning_complete.log"
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["simplejob"],    
    shell:
        "touch finished_binning finished_pipeline "
