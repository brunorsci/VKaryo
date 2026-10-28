"""
VKaryo - Virtual Karyotyping Workflow for Leukemia
====================================================

Comprehensive genomic analysis pipeline for leukemia samples combining:
- RNA-seq: Fusion detection (Arriba, GeneFuse) + Gene expression
- WES: SNV/Indel calling (Mutect2 + VarDict), CNV detection, SV calling

WORKFLOW STRUCTURE:
1. Quality Control & Alignment
   - FASTP (QC + trimming)
   - BWA-MEM (DNA) / STAR (RNA)
   - MarkDuplicates + BQSR

2. Variant Detection
   - SNVs/Indels: Mutect2 (primary) + VarDict (MRD)
   - Fusions: Arriba + GeneFuse
   - SVs: MANTA + DELLY → SURVIVOR consensus
   - CNVs: CNVkit + PureCN

3. Annotation & Filtering
   - Germline filtering (gnomAD)
   - COSMIC annotation
   - Driver gene prioritization

4. Reporting
   - Per-sample HTML reports
   - Cohort dashboard
   - JSON exports

Version: 2.0
"""

print("##############################################")
print("-----------    RUNNING VKARYO    ------------" )
print("------- The Snake Is Going To Smoke! --------" )
print("##############################################")

configfile: "config.yaml"

import pandas as pd

#(SAMPLES,) = glob_wildcards("samples")
(SAMPLES,) = glob_wildcards(config["input"]["exome"] + "{sample}_1.fastq.gz")
(SAMPLES_RNA,) = glob_wildcards(config["input"]["rnaseq"] + "{sample}_1.fastq.gz")
#SAMPLES_RNA=['SAMPLE_RNA_1', 'SAMPLE_RNA_2']

print(SAMPLES_RNA)

# ============================================================
# Sample pairing: DNA and RNA samples have different IDs
# Read mapping from sample_mapping.tsv
# ============================================================
sample_mapping_file = "sample_mapping.tsv"
try:
    sample_map_df = pd.read_csv(sample_mapping_file, sep="\t")
    # Create dictionaries for bidirectional mapping
    DNA_TO_RNA = dict(zip(sample_map_df['DNA_sample'], sample_map_df['RNA_sample']))
    RNA_TO_DNA = dict(zip(sample_map_df['RNA_sample'], sample_map_df['DNA_sample']))
    SAMPLE_TO_CELLLINE = dict(zip(sample_map_df['DNA_sample'], sample_map_df['Cell_line']))

    # SAMPLES_BOTH: DNA samples that have corresponding RNA data available
    SAMPLES_BOTH = [dna for dna in SAMPLES if dna in DNA_TO_RNA and DNA_TO_RNA[dna] in SAMPLES_RNA]

    print(f"INFO: Loaded sample mapping from {sample_mapping_file}")
    print(f"INFO: {len(DNA_TO_RNA)} DNA-RNA pairs defined")
    print(f"INFO: {len(SAMPLES_BOTH)} samples with both DNA and RNA data available")
except FileNotFoundError:
    print(f"WARNING: Sample mapping file {sample_mapping_file} not found. Falling back to name-based matching.")
    DNA_TO_RNA = {}
    RNA_TO_DNA = {}
    SAMPLE_TO_CELLLINE = {}
    SAMPLES_BOTH = list(set(SAMPLES) & set(SAMPLES_RNA))
except Exception as e:
    print(f"WARNING: Error loading sample mapping: {e}. Falling back to name-based matching.")
    DNA_TO_RNA = {}
    RNA_TO_DNA = {}
    SAMPLE_TO_CELLLINE = {}
    SAMPLES_BOTH = list(set(SAMPLES) & set(SAMPLES_RNA))

if config["RNA"]:
    include: "rules/star_aligment.snakefile"
    include: "rules/index_BAM.snakefile"
    include: "rules/arriba.snakefile"
    # Anotação de fusões
    include: "rules/annotate_arriba.snakefile"

if config["DNA"]:
    include: "rules/fastp.snakefile"
    include: "rules/multiqc.snakefile"
    include: "rules/bwa_mem.snakefile"
    include: "rules/indexBAM.snakefile"
    include: "rules/markduplicates.snakefile"
    include: "rules/recal_table.snakefile"
    include: "rules/recal_apply.snakefile"
    include: "rules/manta.snakefile"
    include: "rules/fix_mantaVCF.snakefile"
    include: "rules/annotSV.snakefile"
    include: "rules/delly_sv.snakefile"
    include: "rules/delly_filter.snakefile"
    ## SURVIVOR merge (Manta + Delly)
    include: "rules/survivor.snakefile"
    include: "rules/filtering_vcf.snakefile"
    ##CNVs-aneuploidias
    # NOTA: incluído antes de cnvkitRef, que consome seus outputs
    # (target.bed/antitarget.bed).
    include: "rules/cnvkitTargets.snakefile"
    include: "rules/cnvkitRef.snakefile"
    include: "rules/cnvkitBatch"
    include: "rules/cnvFocalEvents.snakefile"
    include: "rules/cnvArmLevel.snakefile"
    # A figura das duas etapas acima. Fora da rule all de proposito --
    # ver o cabecalho do arquivo.
    include: "rules/cariotipoVirtual.snakefile"
    include: "rules/cnv_scatter.snakefile"
    ## Mutect2 (GATK padrão para cancer) - Primary SNV caller
    include: "rules/mutect2.snakefile"
    ## Filtragem de variantes germinais
    include: "rules/filter_germline.snakefile"
    ## Sample and cohort reports
    include: "rules/sample_reports.snakefile"

    ## ============================================================
    ## Ferramentas opcionais (controladas por config.yaml -> tools)
    ## ============================================================
    # PureCN - Estimativa de pureza tumoral e ploidia
    if config.get("tools", {}).get("purecn", False):
        include: "rules/purecn.snakefile"

    # OptiType - Tipagem HLA
    if config.get("tools", {}).get("optitype", False):
        include: "rules/optitype.snakefile"

    # VarDict - Detecção de MRD (usa BED de genes drivers)
    if config.get("tools", {}).get("vardict", False):
        include: "rules/vardict.snakefile"
        include: "rules/merge_variant_callers.snakefile"







    include: "rules/annotSVDUO.snakefile"

# Novas análises integradas (DNA + RNA)
if config["RNA"] and config["DNA"]:
    include: "rules/validate_cnv_expression.snakefile"
    include: "rules/plot_aneuploidy.snakefile"
    include: "rules/classify_leukemia.snakefile"
    include: "rules/generate_pdf_report.snakefile"
    # Expressão de CNVs (pacote CNVExpression) -- opcional, só roda se
    # tools.cnv_expression=True no config.yaml. Precisa de DNA+RNA pareados
    # (usa CNVkit + STAR de todas as amostras concatenadas).
    if config.get("tools", {}).get("cnv_expression", False):
        include: "rules/cnv_expression.snakefile"

# Database e relatórios consolidados
if config["DNA"]:
    include: "rules/database.snakefile"
    include: "rules/consolidated_reports.snakefile"
    # "Parte 1": concatena os resultados de cada programa de TODAS as
    # amostras num TSV só -- regra de verdade (não script solto), para
    # funcionar rodando o workflow sozinho em outro projeto. Ver
    # rules/consolidate_results.snakefile.
    include: "rules/consolidate_results.snakefile"

# Define amostras com AMBOS DNA e RNA (usado em múltiplos lugares)
# NOTE: SAMPLES_BOTH is now defined earlier in the file using the sample_mapping.tsv file
# SAMPLES_BOTH = list(set(SAMPLES) & set(SAMPLES_RNA)) if config["RNA"] and config["DNA"] else []

# ============================================================
# Outputs opcionais (controlados por config.yaml -> tools)
# ============================================================
PURECN_OUTPUTS = []
if config.get("tools", {}).get("purecn", False) and config.get("DNA", False):
    PURECN_OUTPUTS = [
        expand(config["output"] + "output_vkaryo/purecn/{sample}.rds", sample=SAMPLES),
        expand(config["output"] + "output_vkaryo/purecn/{sample}.pdf", sample=SAMPLES),
        expand(config["output"] + "output_vkaryo/purecn/{sample}_purity_ploidy.txt", sample=SAMPLES),
    ]

OPTITYPE_OUTPUTS = []
if config.get("tools", {}).get("optitype", False) and config.get("DNA", False):
    OPTITYPE_OUTPUTS = [
        expand(config["output"] + "output_vkaryo/optitype/{sample}_hla_typing.txt", sample=SAMPLES),
    ]

VARDICT_OUTPUTS = []
if config.get("tools", {}).get("vardict", False) and config.get("DNA", False):
    VARDICT_OUTPUTS = [
        expand(config["output"] + "output_vkaryo/vardict/{sample}_mrd_filtered.vcf.gz", sample=SAMPLES),
        expand(config["output"] + "output_vkaryo/vardict/{sample}_mrd_summary.txt", sample=SAMPLES),
        expand(config["output"] + "output_vkaryo/merged_variants/{sample}_merged_snvs.vcf.gz", sample=SAMPLES),
        expand(config["output"] + "output_vkaryo/merged_variants/{sample}_high_confidence.vcf.gz", sample=SAMPLES),
    ]

CNV_EXPRESSION_OUTPUTS = []
if config.get("tools", {}).get("cnv_expression", False) and config.get("RNA", False) and config.get("DNA", False):
    CNV_EXPRESSION_OUTPUTS = [
        config["output"] + "output_vkaryo/cnv_expression/cnv_expression_results_filtered.tsv",
        config["output"] + "output_vkaryo/cnv_expression/cnv_expression_results_merged.tsv",
    ]

RNA_OUTPUTS = [
        ### FASTP ###

        ### GENEFUSE ###

        ### STAR ALIGN ###
        expand(config["output"] + "output_vkaryo/star_aligned/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES_RNA),
        expand(config["output"] + "output_vkaryo/star_aligned/{sample}_Aligned.sortedByCoord.out.bai"  , sample=SAMPLES_RNA),

        ### ARRIBA ###
        #expand(config["output"] + "output_vkaryo/arriba/{sample}_fusions.tsv", sample=SAMPLES_RNA),
        #expand(config["output"] + "output_vkaryo/arriba/{sample}_fusions_discarded.tsv", sample=SAMPLES_RNA),

        ### MANTA ###

        ### ARRIBA ANOTADO ###
        expand(config["output"] + "output_vkaryo/arriba_annotated/{sample}_fusions_annotated.tsv", sample=SAMPLES_RNA),
]

DNA_OUTPUTS = [
    # ============================================================
    # TESTE INCREMENTAL: Ativando etapas uma por uma
    # Descomente cada bloco conforme a etapa anterior completar 100/100
    # ============================================================

    ## ETAPA 1: Mutect2 filter - COMPLETA (100/100) ✓
    expand(config["output"] + "output_vkaryo/mutect2/{sample}_filtered.vcf.gz", sample=SAMPLES),

    ## ETAPA 2: Mutect2 pass_only - COMPLETA (100/100) ✓
    expand(config["output"] + "output_vkaryo/mutect2/{sample}_somatic_pass.vcf.gz", sample=SAMPLES),

    ## ETAPA 3: Mutect2 ANNOVAR annotation - COMPLETA (100/100) ✓
    expand(config["output"] + "output_vkaryo/mutect2/annovar/{sample}_mutect2.hg38_multianno.txt", sample=SAMPLES),

    ## ETAPA 4 (ATUAL): Mutect2 summary
    expand(config["output"] + "output_vkaryo/mutect2/{sample}_somatic_summary.txt", sample=SAMPLES),

    ## ETAPA 5: REMOVIDA - classificação somática/germinativa incorporada na ETAPA 4

    ## --- DESATIVADO TEMPORARIAMENTE (reativar após Mutect2 completo) ---

    ## ETAPA 5: SURVIVOR (Merge Manta + Delly SVs) - COMPLETA (100/100) ✓
    expand(config["output"] + "output_vkaryo/survivor/{sample}_merged.vcf.gz", sample=SAMPLES),
    expand(config["output"] + "output_vkaryo/survivor/{sample}_high_confidence.vcf.gz", sample=SAMPLES),
    expand(config["output"] + "output_vkaryo/survivor/{sample}_sv_summary.txt", sample=SAMPLES),
    expand(config["output"] + "output_vkaryo/survivor/{sample}_comparison_matrix.txt", sample=SAMPLES),

    ## ETAPA 6: CNVKit
    expand(config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr", sample=SAMPLES),
    expand(config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns", sample=SAMPLES),
    # Eventos focais de CNV apoiados em evidencia de cobertura
    # (rules/cnvFocalEvents.snakefile). Necessarios porque a razao log2 do
    # CNVkit com referencia plana nao e quantitativamente valida no regime de
    # perdas profundas -- ver o cabecalho daquele arquivo.
    expand(config["output"] + "output_vkaryo/cnv_focal/{sample}_focal_events.tsv", sample=SAMPLES),
    config["output"] + "output_vkaryo/cnv_focal/all_samples_focal_events.tsv",
    # Desequilibrio em nivel de braco (rules/cnvArmLevel.snakefile) -- cobre o
    # regime hemizigotico, que a etapa focal nao alcanca.
    expand(config["output"] + "output_vkaryo/cnv_arm/{sample}_arm_level.tsv", sample=SAMPLES),
    config["output"] + "output_vkaryo/cnv_arm/all_samples_arm_level.tsv",

    ## Sample reports (HTML + JSON)
    #expand(config["output"] + "output_vkaryo/reports/{sample}_vkaryo_report.html", sample=SAMPLES),
    #expand(config["output"] + "output_vkaryo/reports/{sample}_summary.json", sample=SAMPLES),

    ## Cohort dashboard
    #config["output"] + "output_vkaryo/reports/cohort_dashboard.html",

    ## Database centralizado
    #config["output"] + "output_vkaryo/database/vkaryo_results.db",

    ## Relatórios consolidados
    #expand(config["output"] + "output_vkaryo/reports_consolidated/{sample}_vkaryo_consolidated.html",
    #       sample=SAMPLES_BOTH if (config["RNA"] and config["DNA"] and len(SAMPLES_BOTH) > 0) else SAMPLES),
    #config["output"] + "output_vkaryo/reports_consolidated/cohort_summary.html",

    ## Dashboard ready flag
    #config["output"] + "output_vkaryo/dashboard/.dashboard_ready",
    ]

# Outputs integrados DNA+RNA (quando ambos estão ativados)
DNA_RNA_INTEGRATED = []
if config["RNA"] and config["DNA"] and len(SAMPLES_BOTH) > 0:
    # SAMPLES_BOTH já definido no início do arquivo
    # Só gera outputs integrados se houver amostras com ambos DNA e RNA
    DNA_RNA_INTEGRATED = [
        ### CNVs validados com expressão ###
        expand(config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_validated.tsv", sample=SAMPLES_BOTH),
        expand(config["output"] + "output_vkaryo/cnv_validated/{sample}_cnv_high_confidence.tsv", sample=SAMPLES_BOTH),

        ### Plots de aneuploidias ###
        expand(config["output"] + "output_vkaryo/plots/{sample}_aneuploidy_ideogram.pdf", sample=SAMPLES_BOTH),
        expand(config["output"] + "output_vkaryo/plots/{sample}_cnv_circular.pdf", sample=SAMPLES_BOTH),
        expand(config["output"] + "output_vkaryo/plots/{sample}_aneuploidy_summary.txt", sample=SAMPLES_BOTH),

        ### Classificação ###
        expand(config["output"] + "output_vkaryo/classification/{sample}_classification.json", sample=SAMPLES_BOTH),

        ### Relatório PDF final ###
        expand(config["output"] + "output_vkaryo/reports/{sample}_vkaryo_report.pdf", sample=SAMPLES_BOTH),
    ]
else:
    if config["RNA"] and config["DNA"]:
        print(f"WARNING: No samples with both DNA and RNA data (SAMPLES_BOTH is empty). Skipping integrated outputs.")

# Scatter de CNV colorido (azul=ganho/vermelho=delecao) por amostra x gene de
# interesse -- genes vem de reference/genes_interesse_por_amostra.tsv (uma
# linha por par amostra-gene; amostras sem alteracao-chave documentada nao
# entram, nao ha gene de interesse pra destacar nelas).
import os as _os
GENE_TO_CHROM = {}
GENE_TO_CELLLINE = {}
SAMPLE_TO_GENES = {}       # dna_sample -> [gene, ...] (pode ser vazio)
SAMPLE_TO_CELLLINE = {}    # dna_sample -> cell_line (só as 70 com alteração-chave documentada)
CNV_SCATTER_OUTPUTS = []
_genes_interesse_path = config.get("cnv", {}).get("genes_interesse")
if _genes_interesse_path and _os.path.exists(_genes_interesse_path):
    import csv as _csv
    with open(_genes_interesse_path) as _f:
        for _row in _csv.DictReader(_f, delimiter="\t"):
            _key = (_row["dna_sample"], _row["gene"])
            GENE_TO_CHROM[_key] = _row["chrom"]
            GENE_TO_CELLLINE[_key] = _row["cell_line"]
            SAMPLE_TO_GENES.setdefault(_row["dna_sample"], []).append(_row["gene"])
            SAMPLE_TO_CELLLINE[_row["dna_sample"]] = _row["cell_line"]
            CNV_SCATTER_OUTPUTS.append(
                config["output"] + f"output_vkaryo/cnv_scatter/{_row['dna_sample']}_{_row['gene']}_scatter.png"
            )
else:
    print(f"WARNING: {_genes_interesse_path} não encontrado -- pulando outputs de cnv_scatter_plot.")

# Cariótipo virtual (todos os 24 cromossomos, uma amostra por figura) -- para
# as 100 amostras DNA+RNA do painel, não só as 70 com alteração-chave
# documentada (as 70 ganham destaque do(s) gene(s) de interesse; as demais
# saem sem destaque, só o panorama).
CNV_GENOME_SCATTER_OUTPUTS = [
    config["output"] + f"output_vkaryo/cnv_scatter/{_s}_genoma_completo.png" for _s in SAMPLES_BOTH
] if config["RNA"] and config["DNA"] and len(SAMPLES_BOTH) > 0 else []

rule all:
    resources:
        mem_mb=54000,
        cores=30
    input:
        (RNA_OUTPUTS if config["RNA"] == True else []) +
        (DNA_OUTPUTS if config["DNA"] == True else []) +
        DNA_RNA_INTEGRATED +
        PURECN_OUTPUTS +
        OPTITYPE_OUTPUTS +
        VARDICT_OUTPUTS +
        CNV_EXPRESSION_OUTPUTS +
        CNV_SCATTER_OUTPUTS +
        CNV_GENOME_SCATTER_OUTPUTS
