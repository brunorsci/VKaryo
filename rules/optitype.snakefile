rule optitype_extract_hla_reads:
    """
    Extrai reads que mapeiam para região HLA
    """
    input:
        bam = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        bai = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bai"
    output:
        fq1 = config["output"] + "output_vkaryo/optitype/{sample}_hla_1.fastq",
        fq2 = config["output"] + "output_vkaryo/optitype/{sample}_hla_2.fastq"
    params:
        # Região HLA chr6:28,510,120-33,480,577 (hg38)
        region = "chr6:28510120-33480577",
        outdir = lambda wildcards, output: os.path.dirname(output.fq1)
    log:
        config["output"] + "logs/optitype_extract/{sample}.log"
    threads: 4
    conda:
        "../envs/optitype.yaml"
    shell:
        """
        # Ensure output directory exists
        mkdir -p {params.outdir}

        # Extract HLA reads and convert to FASTQ
        samtools view -b {input.bam} {params.region} | \
        samtools sort -n - | \
        samtools fastq -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -n \
        2>&1 | tee {log}

        # Verify output files were created
        if [ ! -s {output.fq1} ] || [ ! -s {output.fq2} ]; then
            echo "WARNING: No HLA reads extracted for sample {wildcards.sample}" | tee -a {log}
            # Create empty placeholder files to avoid downstream failure
            touch {output.fq1} {output.fq2}
        fi
        """


rule optitype_run:
    """
    HLA typing com OptiType
    """
    input:
        fq1 = config["output"] + "output_vkaryo/optitype/{sample}_hla_1.fastq",
        fq2 = config["output"] + "output_vkaryo/optitype/{sample}_hla_2.fastq"
    output:
        result = config["output"] + "output_vkaryo/optitype/{sample}/{sample}_result.tsv",
        coverage = config["output"] + "output_vkaryo/optitype/{sample}/{sample}_coverage_plot.pdf"
    params:
        outdir = config["output"] + "output_vkaryo/optitype/{sample}",
        prefix = "{sample}",
        sequencing = "dna"  # dna ou rna
    log:
        config["output"] + "logs/optitype/{sample}.log"
    threads: 4
    conda:
        "../envs/optitype.yaml"
    shell:
        """
        OptiTypePipeline.py \
            -i {input.fq1} {input.fq2} \
            --{params.sequencing} \
            -v \
            -o {params.outdir} \
            -p {params.prefix} \
            2>&1 | tee {log}
        """


rule optitype_format_result:
    """
    Formata resultado do OptiType para relatório
    """
    input:
        result = config["output"] + "output_vkaryo/optitype/{sample}/{sample}_result.tsv"
    output:
        formatted = config["output"] + "output_vkaryo/optitype/{sample}_hla_typing.txt"
    log:
        config["output"] + "logs/optitype_format/{sample}.log"
    run:
        import pandas as pd

        # Lê resultado OptiType
        df = pd.read_csv(input.result, sep='\t')

        with open(output.formatted, 'w') as f:
            f.write("="*60 + "\n")
            f.write("HLA TYPING - OptiType\n")
            f.write("="*60 + "\n\n")
            f.write(f"Sample: {wildcards.sample}\n\n")

            # Pega melhor solução
            best = df.iloc[0]

            f.write("HLA Alleles (4-digit resolution):\n")
            f.write(f"  HLA-A: {best['A1']}, {best['A2']}\n")
            f.write(f"  HLA-B: {best['B1']}, {best['B2']}\n")
            f.write(f"  HLA-C: {best['C1']}, {best['C2']}\n")
            f.write(f"\nReads: {int(best['Reads'])}\n")
            f.write(f"Objective: {best['Objective']:.2f}\n\n")

            f.write("="*60 + "\n")
            f.write("RELEVÂNCIA CLÍNICA\n")
            f.write("="*60 + "\n")
            f.write("Uso para:\n")
            f.write("- Matching de doador para transplante alogênico\n")
            f.write("- Risco de GvHD (Graft-versus-Host Disease)\n")
            f.write("- Farmacogenética (HLA-B*15:02, etc.)\n")
