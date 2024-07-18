configfile: "../config.yaml"

(SAMPLES,) = glob_wildcards(config["input"]["samples"] + "{sample}.fastq")

rule all:
    input:
        "output/genefuse/{sample}_genefuse_report.html",
        "output/genefuse/{sample}_fusions.txt",
        # genefuse
        #config["genefuse"]["dir"],

rul/ genefuse:
    input:
        fastq1="output/genefuse/{sample}.fastq"
        #fastq1="output/genefuse/{sample}_R1.fastq",
        #fastq2="output/genefuse/{sample}_R2.fastq",
        #config=config["genefuse"]["dir"] + "genes.csv",
        reference="/media/prospecmol/disk4/VKaryo/dev/ref/GRCh38_genome.fa",
    output:
        html="output/genefuse/{sample}_genefuse_report.html",
        fusions="output/genefuse/{sample}_fusions.txt",
    log:
        "logs/{sample}_genefuse.log",
    params:
        # optional parameters
        extra="",
    shell:
        """
        
        """