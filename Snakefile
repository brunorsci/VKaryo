# THE SNAKE IS GOING TO SMOKE #
# Config file
configfile: "config.yaml"
# To get the samples names form sample dir
(SAMPLES,) = glob_wildcards(config["input"]["fastq"] + "{sample}.fastq")


# To import the rules for each tool
#include: "rules/fastqc.snakefile"
#include: "rules/fastp.snakefile"
include: "rules/genefuse.snakefile"



# rule all for all output files
rule all:
    input:
        #expand("results/QC/{sample}_fastqc.html", sample=SAMPLES),     # output from fastQC
        #expand("results/QC/{sample}_fastp.html", sample=SAMPLES),      # output from fastp
        #expand("results/QC/{sample}_fastp.json", sample=SAMPLES),      # output from fastp
        #expand("results/trimmed/{sample}_fastp.fastq", sample=SAMPLES),# output from fastp
        expand("output/genefuse/{sample}_genefuse_report.html", sample=SAMPLES),
        expand("output/genefuse/{sample}_fusions.txt", sample=SAMPLES)



