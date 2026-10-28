rule cnv_scatter_plot:
    """
    Scatter de CNV (log2 por bin) do cromossomo inteiro onde está o gene de
    interesse daquela amostra, no estilo cnvkit.py scatter, mas com segmento
    colorido por estado (azul=ganho, vermelho=deleção) e rótulo do gene
    ancorado no segmento que de fato o sobrepõe -- ver
    scripts/plot_cnv_scatter_colored.py para o porquê dessas duas escolhas.
    Genes de interesse por amostra vêm de reference/genes_interesse_por_amostra.tsv
    (uma linha por par amostra-gene; uma amostra pode ter mais de um gene).
    """
    input:
        cnr = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr",
        cns = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns",
        gene_coords = config["cnv"]["gene_coords"]
    output:
        png = config["output"] + "output_vkaryo/cnv_scatter/{sample}_{gene}_scatter.png"
    params:
        chrom = lambda wc: GENE_TO_CHROM[(wc.sample, wc.gene)],
        title = lambda wc: f"{GENE_TO_CELLLINE.get((wc.sample, wc.gene), wc.sample)} ({wc.sample}) -- {GENE_TO_CHROM[(wc.sample, wc.gene)]}",
        # mesmo corte ja usado em todo o resto do pipeline (config.yaml, tambem
        # usado por scripts/plot_aneuploidy_ideogram.py) -- nao hardcoded aqui
        gain_threshold = config["cnv"]["min_log2_ratio_amp"],
        loss_threshold = config["cnv"]["max_log2_ratio_del"]
    log:
        config["output"] + "logs/cnv_scatter/{sample}_{gene}.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        # NOTA: o site-packages do matplotlib deste ambiente conda foi compilado
        # contra uma libstdc++ mais nova (CXXABI_1.3.15) do que a do sistema
        # (so tem ate 1.3.13) -- sem isso, "import matplotlib" falha com
        # "version CXXABI_1.3.15 not found", mesmo com o ambiente certo ativado,
        # porque o linker dinamico acha a libstdc++ do sistema primeiro. Aponta
        # explicitamente pra lib do proprio ambiente ja ativado (via sys.prefix).
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        python3 scripts/plot_cnv_scatter_colored.py \
            --cnr {input.cnr} \
            --cns {input.cns} \
            --chrom {params.chrom} \
            --gene {wildcards.gene} \
            --gene-coords {input.gene_coords} \
            --gain-threshold {params.gain_threshold} \
            --loss-threshold {params.loss_threshold} \
            --output {output.png} \
            --title "{params.title}" \
            2>&1 | tee {log}
        """


def _genome_scatter_gene_args(wc):
    genes = SAMPLE_TO_GENES.get(wc.sample, [])
    return " ".join(f"--gene {g}" for g in genes)


rule cnv_genome_scatter_plot:
    """
    Cariótipo virtual -- todos os 24 cromossomos de uma amostra num grid só,
    mesmo estilo de ponto+segmento colorido da cnv_scatter_plot, mas pensado
    pra achado GRANDE (braço/cromossomo inteiro: deleção de braço,
    monossomia, trissomia), não pico focal -- ver
    scripts/plot_cnv_genome_scatter.py. Roda pra todas as amostras
    DNA+RNA do painel, não só as 70 com alteração-chave documentada; quando a
    amostra tem gene(s) de interesse (reference/genes_interesse_por_amostra.tsv),
    eles saem destacados também.

    Corte de cor usa config.yaml:cnv.karyotype_gain_threshold/
    karyotype_loss_threshold -- PARÂMETRO SEPARADO do min_log2_ratio_amp/
    max_log2_ratio_del usado no achado focal por gene (cnv_scatter_plot
    acima) e no ideograma: aqui o objetivo é alteração estrutural/
    cromossômica, não CNV pontual (já existe tabela de CNV pra isso), então
    não faz sentido herdar o mesmo corte. Valor pré-definido ainda pendente
    de validação contra amostras com cariótipo anotado -- por isso é
    parâmetro fácil de ajustar no config, não hardcoded.
    """
    input:
        cnr = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr",
        cns = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns",
        gene_coords = config["cnv"]["gene_coords"]
    output:
        png = config["output"] + "output_vkaryo/cnv_scatter/{sample}_genoma_completo.png"
    params:
        title = lambda wc: f"{SAMPLE_TO_CELLLINE.get(wc.sample, wc.sample)} ({wc.sample}) -- cariótipo virtual",
        # corte PRÓPRIO (não o min_log2_ratio_amp/max_log2_ratio_del do achado
        # focal por gene) -- ver comentário em config.yaml:cnv.karyotype_*:
        # regime "germline samples" do próprio CNVkit, ainda pendente de
        # validação contra amostras com cariótipo anotado
        gain_threshold = config["cnv"]["karyotype_gain_threshold"],
        loss_threshold = config["cnv"]["karyotype_loss_threshold"],
        gene_args = _genome_scatter_gene_args
    log:
        config["output"] + "logs/cnv_scatter/{sample}_genoma_completo.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        export LD_LIBRARY_PATH="$(python3 -c 'import sys; print(sys.prefix)')/lib:${{LD_LIBRARY_PATH:-}}"
        python3 scripts/plot_cnv_genome_scatter.py \
            --cnr {input.cnr} \
            --cns {input.cns} \
            --gene-coords {input.gene_coords} \
            {params.gene_args} \
            --gain-threshold {params.gain_threshold} \
            --loss-threshold {params.loss_threshold} \
            --output {output.png} \
            --title "{params.title}" \
            2>&1 | tee {log}
        """
