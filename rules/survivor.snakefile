rule survivor_prepare_file_list:
    """
    Prepara lista de VCFs para SURVIVOR merge
    Lista contém: Manta candidateSV e Delly SVs
    NOTA: SURVIVOR 1.0.7 tem problemas lendo VCFs comprimidos, então descomprimimos antes
    """
    input:
        manta_sv = config["output"] + "output_vkaryo/manta/{sample}/results/variants/candidateSV.vcf.gz",
        delly_sv = config["output"] + "output_vkaryo/delly/svs/{sample}_sv.vcf.gz"
    output:
        file_list = config["output"] + "output_vkaryo/survivor/{sample}_sv_files.txt",
        manta_vcf = temp(config["output"] + "output_vkaryo/survivor/{sample}_manta.vcf"),
        delly_vcf = temp(config["output"] + "output_vkaryo/survivor/{sample}_delly.vcf")
    shell:
        """
        # Descomprimir VCFs (SURVIVOR 1.0.7 tem bug lendo VCFs comprimidos)
        zcat {input.manta_sv} > {output.manta_vcf}
        zcat {input.delly_sv} > {output.delly_vcf}

        # Cria lista de VCFs descomprimidos para merge
        echo {output.manta_vcf} > {output.file_list}
        echo {output.delly_vcf} >> {output.file_list}
        """


rule survivor_merge:
    """
    SURVIVOR merge de variantes estruturais de múltiplos callers
    Parâmetros:
    - max_distance: 1000bp (distância máxima entre breakpoints)
    - min_callers: 1 (mínimo de callers suportando SV)
    - take_type: 1 (usar tipo de SV de caller com maior suporte)
    - take_strand: 1 (usar strand de caller com maior suporte)
    - estimate_distance: -1 (não estimar distância entre SVs)
    - min_size: 50 (tamanho mínimo de SV)
    """
    input:
        file_list = config["output"] + "output_vkaryo/survivor/{sample}_sv_files.txt",
        manta_vcf = config["output"] + "output_vkaryo/survivor/{sample}_manta.vcf",
        delly_vcf = config["output"] + "output_vkaryo/survivor/{sample}_delly.vcf"
    output:
        merged_vcf = config["output"] + "output_vkaryo/survivor/{sample}_merged.vcf",
        merged_vcf_gz = config["output"] + "output_vkaryo/survivor/{sample}_merged.vcf.gz",
        merged_idx = config["output"] + "output_vkaryo/survivor/{sample}_merged.vcf.gz.tbi"
    params:
        max_distance = 1000,      # 1kb max distance between breakpoints
        min_callers = 1,          # Minimum 1 caller (will filter later)
        take_type = 1,            # Use SV type from caller with most support
        take_strand = 1,          # Use strand from caller with most support
        estimate_distance = -1,   # Don't estimate distance
        min_size = 50             # Minimum SV size 50bp
    log:
        config["output"] + "logs/survivor_merge/{sample}.log"
    conda:
        "../envs/survivor.yaml"
    threads: 4
    shell:
        """
        echo "##############################################"
        echo "-----   SURVIVOR Merge SVs      ------------"
        echo "##############################################"

        # SURVIVOR merge
        SURVIVOR merge {input.file_list} \
            {params.max_distance} \
            {params.min_callers} \
            {params.take_type} \
            {params.take_strand} \
            {params.estimate_distance} \
            {params.min_size} \
            {output.merged_vcf}.unsorted \
            2>&1 | tee {log}

        # Sort VCF (SURVIVOR output pode não estar ordenado)
        bcftools sort -O v -o {output.merged_vcf}.sorted {output.merged_vcf}.unsorted
        rm {output.merged_vcf}.unsorted

        # Fix SUPP field type (SURVIVOR incorretamente define como String ao invés de Integer)
        sed 's/ID=SUPP,Number=1,Type=String/ID=SUPP,Number=1,Type=Integer/' {output.merged_vcf}.sorted > {output.merged_vcf}
        rm {output.merged_vcf}.sorted

        # Compress and index
        bgzip -c {output.merged_vcf} > {output.merged_vcf_gz}
        tabix -f -p vcf {output.merged_vcf_gz}
        """


rule survivor_filter_multi_caller:
    """
    Filtra SVs suportadas por múltiplos callers (high confidence)
    Gera dois VCFs:
    - high_confidence: Suportado por 2+ callers (Manta + Delly)
    - single_caller: Suportado por apenas 1 caller
    """
    input:
        merged_vcf = config["output"] + "output_vkaryo/survivor/{sample}_merged.vcf.gz"
    output:
        high_conf = config["output"] + "output_vkaryo/survivor/{sample}_high_confidence.vcf.gz",
        high_conf_idx = config["output"] + "output_vkaryo/survivor/{sample}_high_confidence.vcf.gz.tbi",
        single_caller = config["output"] + "output_vkaryo/survivor/{sample}_single_caller.vcf.gz",
        single_caller_idx = config["output"] + "output_vkaryo/survivor/{sample}_single_caller.vcf.gz.tbi"
    log:
        config["output"] + "logs/survivor_filter/{sample}.log"
    conda:
        "../envs/survivor.yaml"
    shell:
        """
        echo "##############################################"
        echo "-----   Filtering Multi-Caller SVs  --------"
        echo "##############################################"

        # Check if input VCF is valid
        if ! bcftools view -h {input.merged_vcf} > /dev/null 2>&1; then
            echo "WARNING: Invalid or corrupted VCF file, creating empty outputs" | tee {log}
            # Create empty VCFs with proper headers
            bcftools view -h {input.merged_vcf} 2>/dev/null | bgzip > {output.high_conf} || \
                echo -e "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO" | bgzip > {output.high_conf}
            tabix -p vcf {output.high_conf} || tabix -f -p vcf {output.high_conf}
            cp {output.high_conf} {output.single_caller}
            tabix -p vcf {output.single_caller} || tabix -f -p vcf {output.single_caller}
            exit 0
        fi

        # Extract SVs supported by 2+ callers (high confidence)
        bcftools view -i 'INFO/SUPP>=2' {input.merged_vcf} -Oz -o {output.high_conf} 2>&1 | tee {log} || true
        tabix -f -p vcf {output.high_conf}

        # Extract SVs supported by single caller
        bcftools view -i 'INFO/SUPP==1' {input.merged_vcf} -Oz -o {output.single_caller} 2>&1 | tee -a {log} || true
        tabix -f -p vcf {output.single_caller}

        # Log statistics
        echo "High confidence (2+ callers): $(bcftools view -H {output.high_conf} 2>/dev/null | wc -l)" | tee -a {log}
        echo "Single caller: $(bcftools view -H {output.single_caller} 2>/dev/null | wc -l)" | tee -a {log}
        """


rule survivor_stats:
    """
    Gera estatísticas dos SVs merged
    Usa SURVIVOR stats para análise detalhada
    """
    input:
        high_conf = config["output"] + "output_vkaryo/survivor/{sample}_high_confidence.vcf.gz",
        merged = config["output"] + "output_vkaryo/survivor/{sample}_merged.vcf.gz"
    output:
        stats = config["output"] + "output_vkaryo/survivor/{sample}_sv_stats.txt",
        summary = config["output"] + "output_vkaryo/survivor/{sample}_sv_summary.txt"
    params:
        min_size = 50,
        max_size = 100000000  # 100Mb
    log:
        config["output"] + "logs/survivor_stats/{sample}.log"
    conda:
        "../envs/survivor.yaml"
    shell:
        """
        echo "##############################################"
        echo "-----   SURVIVOR Statistics      -----------"
        echo "##############################################"

        # Basic statistics usando SURVIVOR stats (pode falhar, não é crítico)
        SURVIVOR stats {input.merged} \
            {params.min_size} \
            {params.max_size} \
            {output.stats} \
            2>&1 | tee {log} || true

        # Ensure stats file exists (create empty if SURVIVOR stats failed)
        touch {output.stats}

        # Summary manual
        {{
            echo "========================================"
            echo "SURVIVOR SV MERGE SUMMARY"
            echo "Sample: {wildcards.sample}"
            echo "========================================"
            echo ""

            total_merged=$(bcftools view -H {input.merged} 2>/dev/null | wc -l || echo 0)
            high_conf_count=$(bcftools view -H {input.high_conf} 2>/dev/null | wc -l || echo 0)

            echo "Total merged SVs: $total_merged"
            echo "High confidence (2+ callers): $high_conf_count"
            echo ""

            if [ "$high_conf_count" -gt 0 ]; then
                echo "SV Types in High Confidence:"
                bcftools query -f '%SVTYPE\n' {input.high_conf} 2>/dev/null | sort | uniq -c | awk '{{print "  " $2 ": " $1}}' || echo "  N/A"
                echo ""
                echo "SV Size Distribution (High Confidence):"
                bcftools query -f '%SVLEN\n' {input.high_conf} 2>/dev/null | awk '{{
                    len = ($1 < 0) ? -$1 : $1
                    if (len < 100) small++
                    else if (len < 1000) medium++
                    else if (len < 10000) large++
                    else xlarge++
                }} END {{
                    print "  <100bp: " (small ? small : 0)
                    print "  100bp-1kb: " (medium ? medium : 0)
                    print "  1kb-10kb: " (large ? large : 0)
                    print "  >10kb: " (xlarge ? xlarge : 0)
                }}' || echo "  N/A"
                echo ""
                echo "Chromosome Distribution:"
                bcftools query -f '%CHROM\n' {input.high_conf} 2>/dev/null | sort | uniq -c | awk '{{print "  " $2 ": " $1}}' || echo "  N/A"
            else
                echo "No high confidence SVs found"
            fi
        }} > {output.summary}

        cat {output.summary} | tee -a {log}
        """


rule survivor_comparison_matrix:
    """
    Gera matriz de comparação entre callers
    Mostra overlap e concordância entre Manta e Delly
    """
    input:
        merged = config["output"] + "output_vkaryo/survivor/{sample}_merged.vcf.gz"
    output:
        matrix = config["output"] + "output_vkaryo/survivor/{sample}_comparison_matrix.txt"
    log:
        config["output"] + "logs/survivor_comparison/{sample}.log"
    conda:
        "../envs/survivor.yaml"
    shell:
        """
        echo "##############################################"
        echo "-----   Caller Comparison Matrix   ---------"
        echo "##############################################"

        # Check variant count first
        variant_count=$(zcat {input.merged} 2>/dev/null | grep -v "^#" | wc -l || echo 0)

        {{
            echo "========================================"
            echo "CALLER COMPARISON MATRIX"
            echo "Sample: {wildcards.sample}"
            echo "========================================"
            echo ""
            echo "Support Patterns (SUPP_VEC):"
            echo "  11 = Both Manta and Delly"
            echo "  10 = Manta only"
            echo "  01 = Delly only"
            echo ""

            if [ "$variant_count" -eq 0 ]; then
                echo "No variants found in merged VCF"
                echo ""
                echo "Total SVs: 0"
                echo ""
                echo "Breakdown:"
                echo "  Both callers: 0 (0.0%)"
                echo "  Manta only: 0 (0.0%)"
                echo "  Delly only: 0 (0.0%)"
            else
                # Extract SUPP_VEC values and count them
                echo "Support Distribution:"
                zcat {input.merged} | grep -v "^#" | grep -oP 'SUPP_VEC=\K[^;]+' | sort | uniq -c | while read count vec; do
                    echo "  $vec: $count"
                done

                echo ""

                # Count totals
                total=$(zcat {input.merged} | grep -v "^#" | wc -l)
                both=$(zcat {input.merged} | grep -v "^#" | grep -oP 'SUPP_VEC=\K[^;]+' | grep -c "^11$" || echo 0)
                manta=$(zcat {input.merged} | grep -v "^#" | grep -oP 'SUPP_VEC=\K[^;]+' | grep -c "^10$" || echo 0)
                delly=$(zcat {input.merged} | grep -v "^#" | grep -oP 'SUPP_VEC=\K[^;]+' | grep -c "^01$" || echo 0)

                echo "Total SVs: $total"
                echo ""
                echo "Breakdown:"
                echo "  Both callers: $both"
                echo "  Manta only: $manta"
                echo "  Delly only: $delly"
            fi
        }} > {output.matrix} 2>&1

        cat {output.matrix} | tee {log}
        """
