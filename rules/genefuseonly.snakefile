configfile: "config.yaml"

SAMPLES, = glob_wildcards(config["input"]["fastq"] + "{sample}.fastq")

rule all:
    input:
        expand("output/genefuse/{sample}_genefuse_report.html", sample=SAMPLES),
        expand("output/genefuse/{sample}_fusions.json", sample=SAMPLES)

rule genefuse:
    input:
        fastq1="samples/R1.fastq",
        reference="reference_genome/hg38.fa"
    output:
        html="output/genefuse/{sample}_genefuse_report.html",
        fusions="output/genefuse/{sample}_fusions.json"
    log:
        "output/genefuse/{sample}_genefuse.log"
    params:
        config="reference_genome/cancer.hg38.csv",
        extra=""  # Parâmetros extras opcionais
    shell:
        """
        genefuse --read1 {input.fastq1} --fusion {params.config} --ref  {input.reference}  --html {output.html}  --json {output.fusions}
        """
