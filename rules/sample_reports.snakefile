# Define samples with RNA data for conditional inputs
# Note: DNA and RNA samples have different IDs, so we need to use the DNA_TO_RNA mapping
_samples_with_rna = set(SAMPLES_BOTH) if config.get("RNA") and config.get("DNA") else set()


def get_fusion_input(wildcards):
    """Return fusion file path if sample has RNA data, else empty list

    Uses DNA_TO_RNA mapping to find the correct RNA sample ID
    """
    dna_sample = wildcards.sample

    # Check if this DNA sample has corresponding RNA data
    if dna_sample in SAMPLES_BOTH:
        # Get the RNA sample ID from the mapping
        rna_sample = DNA_TO_RNA.get(dna_sample, dna_sample)
        return config["output"] + f"output_vkaryo/arriba_annotated/{rna_sample}_fusions_annotated.tsv"
    return []


def get_purecn_input(wildcards):
    """Return PureCN purity/ploidy file if PureCN is enabled, else empty list"""
    if config.get("tools", {}).get("purecn", False):
        return config["output"] + f"output_vkaryo/purecn/{wildcards.sample}_purity_ploidy.txt"
    return []


rule generate_sample_report:
    """
    Generate comprehensive HTML report for individual sample
    Includes SNVs, fusions (if RNA available), CNVs, QC metrics
    """
    input:
        # SNVs (Mutect2 only - VarDict disabled)
        driver_vcf = config["output"] + "output_vkaryo/mutect2/{sample}_driver_mutations.vcf.gz",
        somatic_vcf = config["output"] + "output_vkaryo/mutect2/{sample}_somatic_filtered.vcf.gz",

        # Fusions (conditional - only if sample has RNA data)
        fusions_tsv = get_fusion_input,

        # CNVs
        cnr_file = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cnr",
        cns_file = config["output"] + "output_vkaryo/cnvkit/{sample}_markduplicates_recal.cns",

        # QC
        contamination_file = config["output"] + "output_vkaryo/mutect2/{sample}_contamination.table",
        purity_file = get_purecn_input,  # Conditional - only if PureCN enabled

        # BAM (for coverage stats)
        bam_file = config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam"
    output:
        html = config["output"] + "output_vkaryo/reports/{sample}_vkaryo_report.html",
        json = config["output"] + "output_vkaryo/reports/{sample}_summary.json"
    params:
        template = "templates/sample_report_template.html",
        leukemia_type = lambda wildcards: "Leukemia",
        has_rna = lambda wildcards: wildcards.sample in _samples_with_rna
    log:
        config["output"] + "logs/reports/{sample}_report.log"
    conda:
        "../envs/python_analysis.yaml"
    script:
        "../scripts/generate_sample_report.py"


rule generate_cohort_dashboard:
    """
    Generate interactive cohort dashboard with all samples
    Aggregates statistics, creates visualizations, and provides filtering
    """
    input:
        summaries = lambda wildcards: expand(
            config["output"] + "output_vkaryo/reports/{sample}_summary.json",
            sample=SAMPLES
        )
    output:
        dashboard = config["output"] + "output_vkaryo/reports/cohort_dashboard.html"
    log:
        config["output"] + "logs/reports/cohort_dashboard.log"
    conda:
        "../envs/python_analysis.yaml"
    script:
        "../scripts/generate_cohort_dashboard.py"
