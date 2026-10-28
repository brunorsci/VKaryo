# Expressão de CNVs -- integra CNVkit (.cns) com expressão RNA-seq (STAR
# ReadsPerGene.out.tab) de TODAS as amostras DNA+RNA pareadas via o pacote R
# CNVExpression (github.com/brunorsci/CNVExpression_v1, v0.99.0, já
# instalado e testado neste servidor -- usa a versão já instalada, não
# reinstala de lugar nenhum). Opcional: só roda se
# tools.cnv_expression=True no config.yaml (padrão False).
#
# Por que precisa de várias amostras concatenadas, não uma por vez: o
# pacote calcula log fold-change de expressão por gene ENTRE grupos de
# amostras com estados de número de cópias diferentes (CNVRanger::cnvEQTL,
# min_samples=3 por padrão) -- com 1 amostra só não há desvio/média pra
# calcular. Extensão futura (não implementada agora, por decisão do
# usuário): análise por amostra usando expressão média de controles/banco
# de referência em vez de comparação entre amostras do próprio lote.

rule stage_cnv_expression_input:
    """Copia .cns (CNVkit) + ReadsPerGene.out.tab (STAR) de cada amostra
    DNA+RNA pareada pra uma pasta de input dedicada -- snapshot auditável
    do que entrou na análise, separado das pastas gerais do pipeline."""
    input:
        cns=expand(config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns", sample=SAMPLES_BOTH),
        expr=expand(config["output"] + "output_vkaryo/star_aligned/{rna}_ReadsPerGene.out.tab",
                    rna=[DNA_TO_RNA[d] for d in SAMPLES_BOTH]),
    output:
        manifest=config["output"] + "output_vkaryo/cnv_expression/input/.manifest.txt",
    params:
        output_dir=config["output"] + "output_vkaryo",
        sample_mapping="sample_mapping.tsv",
        staging_dir=config["output"] + "output_vkaryo/cnv_expression/input",
    log:
        config["output"] + "logs/cnv_expression/stage_input.log",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/stage_cnv_expression_input.py \
            --output-dir {params.output_dir} \
            --sample-mapping {params.sample_mapping} \
            --staging-dir {params.staging_dir} \
            2>&1 | tee {log}
        ls {params.staging_dir} > {output.manifest}
        """


rule prepare_cnv_expression_tables:
    """Concatena os arquivos copiados no formato esperado pelo pacote
    (sample/Chrom/Start/End/value pra CNV; genes x amostras pra expressão)."""
    input:
        manifest=config["output"] + "output_vkaryo/cnv_expression/input/.manifest.txt",
    output:
        cnv=config["output"] + "output_vkaryo/cnv_expression/cnv_concatenated.tsv",
        expression=config["output"] + "output_vkaryo/cnv_expression/expression_concatenated.tsv",
    params:
        staging_dir=config["output"] + "output_vkaryo/cnv_expression/input",
    log:
        config["output"] + "logs/cnv_expression/prepare_tables.log",
    conda:
        "../envs/python_analysis.yaml"
    shell:
        """
        python scripts/prepare_cnv_expression_input.py \
            --staging-dir {params.staging_dir} \
            --out-cnv {output.cnv} \
            --out-expression {output.expression} \
            2>&1 | tee {log}
        """


rule run_cnv_expression_analysis:
    """Roda cnv_expression_pipeline() do pacote CNVExpression (versão
    instalada no R deste servidor -- NÃO usa conda, o R do sistema já tem
    todas as dependências Bioconductor testadas). Pode demorar alguns
    minutos: carregar o pacote sozinho leva ~40s neste servidor (disco
    lento, dependências Bioconductor profundas), biomaRt pode levar mais
    alguns segundos pra resolver coordenadas de genes via rede."""
    input:
        cnv=config["output"] + "output_vkaryo/cnv_expression/cnv_concatenated.tsv",
        expression=config["output"] + "output_vkaryo/cnv_expression/expression_concatenated.tsv",
    output:
        filtered=config["output"] + "output_vkaryo/cnv_expression/cnv_expression_results_filtered.tsv",
        merged=config["output"] + "output_vkaryo/cnv_expression/cnv_expression_results_merged.tsv",
        raw=config["output"] + "output_vkaryo/cnv_expression/cnv_expression_results_raw.tsv",
    params:
        logfc_cutoff=config.get("cnv_expression", {}).get("logfc_cutoff", 1),
        pvalue_cutoff=config.get("cnv_expression", {}).get("pvalue_cutoff", 0.05),
    log:
        config["output"] + "logs/cnv_expression/run_analysis.log",
    threads: 4
    shell:
        """
        Rscript scripts/run_cnv_expression.R \
            {input.cnv} {input.expression} \
            {output.filtered} {output.merged} {output.raw} \
            {params.logfc_cutoff} {params.pvalue_cutoff} \
            2>&1 | tee {log}
        """
