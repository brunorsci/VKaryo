rule multiqc:
    """
    MultiQC do QC de FASTQ (exoma) -- consolida os relatórios fastp de todas
    as amostras num único report.html.

    O output declarado é o arquivo final (multiqc_report.html), não o
    diretório -- declarar só directory(...) faz o Snakemake considerar a
    regra "pronta" assim que a pasta existe, mesmo vazia, deixando passar
    despercebida uma execução que falhou silenciosamente.
    --no-version-check evita uma checagem de rede no import do MultiQC que,
    neste servidor, adiciona ~100s de latência sem nenhum benefício aqui.
    --no-ai desliga o "AI toolbox" embutido do MultiQC -- mantém o relatório
    determinístico por padrão; qualquer camada de IA no VKaryo é uma feature
    separada e explícita, não algo embutido silenciosamente numa ferramenta
    de terceiros.
    """
    input:
        expand(config["output"] + "output_vkaryo/fastp/exoma/{sample}_fastp.html", sample=SAMPLES),
        expand(config["output"] + "output_vkaryo/fastp/exoma/{sample}_fastp.json", sample=SAMPLES),
    output:
        report=config["output"] + "output_vkaryo/multiqc/exoma/multiqc_report.html",
    params:
        outdir=config["output"] + "output_vkaryo/multiqc/exoma/",
        indir=config["output"] + "output_vkaryo/fastp/exoma/",
    log:
        config["output"] + "output_vkaryo/multiqc/exoma/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    threads: 5
    shell:
        """
        # $CONDA_PREFIX/bin/multiqc explícito, não "multiqc" cru: uma
        # instalação pip do sistema (ex. ~/.local/bin/multiqc) pode ter
        # precedência no PATH mesmo dentro do ambiente conda ativado --
        # o caminho explícito garante que roda o multiqc do ambiente conda.
        $CONDA_PREFIX/bin/multiqc --force --no-version-check --no-ai \
            --outdir {params.outdir} {params.indir} \
            > {log} 2>&1
        """


rule multiqc_rnaseq:
    """MultiQC do QC de FASTQ (RNA-seq) -- mesma lógica da regra acima,
    cobrindo fastp/rnaseq/ (não existia antes)."""
    input:
        expand(config["output"] + "output_vkaryo/fastp/rnaseq/{sample}_fastp.html", sample=SAMPLES_RNA),
        expand(config["output"] + "output_vkaryo/fastp/rnaseq/{sample}_fastp.json", sample=SAMPLES_RNA),
    output:
        report=config["output"] + "output_vkaryo/multiqc/rnaseq/multiqc_report.html",
    params:
        outdir=config["output"] + "output_vkaryo/multiqc/rnaseq/",
        indir=config["output"] + "output_vkaryo/fastp/rnaseq/",
    log:
        config["output"] + "output_vkaryo/multiqc/rnaseq/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    threads: 5
    shell:
        """
        # $CONDA_PREFIX/bin/multiqc explícito, não "multiqc" cru: uma
        # instalação pip do sistema (ex. ~/.local/bin/multiqc) pode ter
        # precedência no PATH mesmo dentro do ambiente conda ativado --
        # o caminho explícito garante que roda o multiqc do ambiente conda.
        $CONDA_PREFIX/bin/multiqc --force --no-version-check --no-ai \
            --outdir {params.outdir} {params.indir} \
            > {log} 2>&1
        """
