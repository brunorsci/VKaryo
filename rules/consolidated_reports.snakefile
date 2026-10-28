rule generate_consolidated_report:
    """
    Gera relatório consolidado HTML/PDF por amostra
    Integra todos os resultados em um único documento
    """
    input:
        # Use integrated DB only if this specific sample has both DNA and RNA data
        db = lambda wildcards: (
            config["output"] + f"output_vkaryo/database/{wildcards.sample}_results_integrated.db"
            if config["RNA"] and config["DNA"] and wildcards.sample in SAMPLES_BOTH
            else config["output"] + f"output_vkaryo/database/{wildcards.sample}_results.db"
        ),
        # Plots opcionais (conditional on tool availability)
        purecn_plot = lambda wildcards: (
            config["output"] + f"output_vkaryo/purecn/{wildcards.sample}.pdf"
            if config.get("tools", {}).get("purecn", False)
            else []
        ),
        # mrd_plot disabled - vardict_track_mrd rule requires ANNOVAR
        # mrd_plot = lambda wildcards: [],
        # Aneuploidy plot only for samples with both DNA and RNA
        aneuploidy_plot = lambda wildcards: (
            config["output"] + f"output_vkaryo/plots/{wildcards.sample}_aneuploidy_ideogram.pdf"
            if config["RNA"] and config["DNA"] and wildcards.sample in SAMPLES_BOTH
            else []
        )
    output:
        html = config["output"] + "output_vkaryo/reports_consolidated/{sample}_vkaryo_consolidated.html",
        pdf = config["output"] + "output_vkaryo/reports_consolidated/{sample}_vkaryo_consolidated.pdf"
    log:
        config["output"] + "logs/consolidated_reports/{sample}.log"
    conda:
        "../envs/reporting.yaml"
    script:
        "../scripts/generate_consolidated_report.py"


rule generate_cohort_summary:
    """
    Gera relatório sumário de todo o cohort
    Overview de todas as amostras analisadas
    """
    input:
        db = config["output"] + "output_vkaryo/database/vkaryo_results.db",
        # Todos os relatórios individuais
        # Usa apenas amostras com dados completos (DNA+RNA se ambos ativos, senão só DNA)
        individual_reports = lambda wildcards: (
            expand(config["output"] + "output_vkaryo/reports_consolidated/{sample}_vkaryo_consolidated.html",
                   sample=list(set(SAMPLES) & set(SAMPLES_RNA)))
            if config["RNA"] and config["DNA"]
            else expand(config["output"] + "output_vkaryo/reports_consolidated/{sample}_vkaryo_consolidated.html",
                       sample=SAMPLES)
        )
    output:
        cohort_html = config["output"] + "output_vkaryo/reports_consolidated/cohort_summary.html",
        cohort_pdf = config["output"] + "output_vkaryo/reports_consolidated/cohort_summary.pdf"
    log:
        config["output"] + "logs/cohort_summary/cohort.log"
    conda:
        "../envs/reporting.yaml"
    script:
        "../scripts/generate_cohort_summary.py"


rule launch_dashboard:
    """
    Script helper para lançar dashboard Streamlit
    Não faz parte do DAG principal, executar manualmente
    """
    input:
        db = config["output"] + "output_vkaryo/database/vkaryo_results.db"
    output:
        flag = config["output"] + "output_vkaryo/dashboard/.dashboard_ready"
    conda:
        "../envs/reporting.yaml"
    shell:
        """
        echo "Dashboard pronto para ser lançado!" > {output.flag}
        echo ""
        echo "=========================================="
        echo "  Para lançar o dashboard interativo:"
        echo "=========================================="
        echo ""
        echo "  conda activate reporting"
        echo "  streamlit run scripts/dashboard_app.py"
        echo ""
        echo "  Database: {input.db}"
        echo ""
        echo "=========================================="
        """
