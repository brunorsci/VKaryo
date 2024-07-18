configfile: "../config.yaml"

# Obter os nomes das amostras do diretório de amostras
SAMPLES, = glob_wildcards(config["input"]["samples"] + "{sample}.fastq")

# Função para determinar se os dados são paired-end
def determine_input_type(sample):
    fastq_files = glob(f"output/genefuse/{sample}*.fastq")
    if len(fastq_files) == 1:
        return f"-1 {fastq_files[0]}"  # Single-end
    elif len(fastq_files) == 2:
        return f"-1 {fastq_files[0]} -2 {fastq_files[1]}"  # Paired-end
    else:
        raise ValueError(f"Número inválido de arquivos FASTQ para {sample}")

# Regra principal que define os arquivos de saída finais
rule all:
    input:
        expand("output/genefuse/{sample}_genefuse_report.html", sample=SAMPLES),
        expand("output/genefuse/{sample}_fusions.txt", sample=SAMPLES)

# Regra para rodar o Genefuse
rule genefuse:
    input:
        # Determina se é paired-end ou single-end
        input_files=lambda wildcards: determine_input_type(wildcards.sample),
        reference="/media/prospecmol/disk4/VKaryo/dev/ref/GRCh38_genome.fa"
    output:
        html="output/genefuse/{sample}_genefuse_report.html",
        fusions="output/genefuse/{sample}_fusions.txt"
    log:
        "logs/{sample}_genefuse.log"
    params:
        config=config["genefuse"]["dir"] + "genes.csv",
        extra=""  # Parâmetros extras opcionais
    shell:
        """
        genefuse -c {params.config} -r {input.reference} {input.input_files} -o {output.html} > {output.fusions} 2> {log}
        """
