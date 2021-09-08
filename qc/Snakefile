# ********************************
# * Snakefile for metqc pipeline *
# ********************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Define logic ****
#qc = config["qc_only"]

#def all_input_reads(qc):
#    if config["qc_only"]:
#        return expand(config["path"]+"{sample}"+config["for"]+".fastq", sample=SAMPLES)
#    else:
#        if config["run_bbmap"]:
#            return expand("output/bbmap/{sample}_bbmapped_1.fastq", sample=SAMPLES)
#        else:
#            if config["run_bmtagger"]:
#                return expand("output/bmtagger/{sample}_bmtagged_1.fastq", sample=SAMPLES)
#            else:
#                return expand("output/prinseq/{sample}_filtered_1.fastq", sample=SAMPLES)

# **** Rules ****

rule all:
    input:
        "results/multiqc_report.html"
#        "results/multiqc_report.html" if config["qc_only"] else "results/multiqc_report_filtered.html",
#        all_input_reads

rule fastqc:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        r1 = "output/qc/{sample}"+"_R1_001_fastqc.html",
        r2 = "output/qc/{sample}"+"_R2_001_fastqc.html"
    conda: "utils/envs/fastqc_env.yaml"
    shell: "fastqc -o output/qc {input.r1} {input.r2}"

rule cutadapt:
    input:
        r1 = config["path"]+"{sample}"+config["for"],
        r2 = config["path"]+"{sample}"+config["rev"]
    output:
        r1 = "output/cutadapt/{sample}_r1_cutadapt.fastq.gz",
        r2 = "output/cutadapt/{sample}_r2_cutadapt.fastq.gz"
    params: pre = "{sample}"
    conda: "utils/envs/cutadapt_env.yaml"
    shell:
            "cutadapt -q 15,10 -m {config[minlength]} --max-n {config[maxn]} -a {config[fwd_adapter]} "
            "-A {config[rev_adapter]} -o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2} > output/cutadapt/logs/{params.pre}.cutadapt.log"

rule index:
    input:
        r1 = expand("output/cutadapt/{sample}_r1_cutadapt.fastq.gz", sample=SAMPLES),
        r2 = expand("output/cutadapt/{sample}_r2_cutadapt.fastq.gz", sample=SAMPLES)
    output: touch("index.done")
    conda: "utils/envs/bbmap_env.yaml"
    shell: "bbmap.sh ref={config[bbmap_ref]} -Xmx24g"

rule bbmap:
    input: "index.done"
    output:
        r1 = "output/final/{sample}_microbe_1.fastq.gz",
        r2 = "output/final/{sample}_microbe_2.fastq.gz"
    params:
        i = "output/cutadapt/{sample}_r#_cutadapt.fastq.gz",
        u = "output/final/{sample}_microbe_#.fastq.gz",
        m = "output/bbmap/{sample}_host_#.fastq.gz",
        pre = "{sample}"
    conda: "utils/envs/bbmap_env.yaml"
    shell:
        "bbmap.sh in={params.i} outu={params.u} outm={params.m} scafstats=output/bbmap/stats/{params.pre}_scafstats.txt ihist=output/bbmap/stats/{params.pre}_ihist.txt statsfile=output/bbmap/stats/{params.pre}_statsfile.txt -Xmx24g covstats=output/bbmap/stats/{params.pre}_constats.txt"

rule multiqc:
    input:
        r1 = expand("output/final/{sample}_microbe_1.fastq.gz", sample=SAMPLES),
        r2 = expand("output/final/{sample}_microbe_2.fastq.gz", sample=SAMPLES),
        r3 = expand("output/qc/{sample}"+"_R1_001_fastqc.html", sample=SAMPLES)
    output: "results/multiqc_report.html"
    conda: "utils/envs/multiqc_env.yaml"
    shell: "multiqc -f output -o results -n multiqc_report.html"
