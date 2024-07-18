rule fastqc_before_single:
    input:
        "samples/{sample}.fastq"
    output:
        "results/QC/{sample}_fastqc.html"
    shell:
        "fastqc -o results/QC {input}"