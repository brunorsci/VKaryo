rule delly_filter:
    """
    Convert Delly SV BCF to VCF and filter for PASS variants
    Note: CNV processing removed - CNVkit handles CNV detection
    """
    input:
        filtered_sv_bcf=config["output"] + "output_vkaryo/delly/svs/{sample}_sv.bcf",
    output:
        filtered_sv_vcf=config["output"] + "output_vkaryo/delly/svs/{sample}_sv.vcf.gz",
    params:
        reference=config["refgenome"]["ref"],
        mappability=config["refgenome"]["map"]
    conda:
        "../envs/bcftools.yaml"
    threads: 5
    shell:
        """
        echo "##############################################"
        echo "--------   Running Filter Delly SV   --------"
        echo "##############################################"

        # Convert SV BCF to VCF and filter for PASS
        bcftools view -i 'FILTER="PASS"' {input.filtered_sv_bcf} -Oz -o {output.filtered_sv_vcf} --threads {threads}
        bcftools index {output.filtered_sv_vcf} --threads {threads}
        """


rule delly_full_vcf:
    """
    Converte o BCF bruto do DELLY (TODOS os eventos, PASS e não-PASS) pra VCF,
    sem filtrar -- alimenta o AnnotSV usado na consolidação (rule annotsvDUO),
    em paralelo ao delly_filter (PASS-only) acima, que SURVIVOR/sansa
    continuam usando sem mudança nenhuma. Mesmo raciocínio do
    mutect2_funcregion_scan/mutect2_exonic_subset_vcf (rules/mutect2.snakefile):
    não descartar não-PASS antes da anotação funcional, porque isso pode
    esconder achado real (caso SOCS1/L-1236 no Mutect2; caso K-562 aqui, onde
    195 eventos PASS do DELLY não apareciam na tabela consolidada só por
    filtro de painel de genes -- ver PLANO_correcao_artigo_e_banco_dados.md).
    """
    input:
        sv_bcf=config["output"] + "output_vkaryo/delly/svs/{sample}_sv.bcf",
    output:
        full_vcf=config["output"] + "output_vkaryo/delly/svs/{sample}_sv_full.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    threads: 5
    shell:
        """
        bcftools view {input.sv_bcf} -Oz -o {output.full_vcf} --threads {threads}
        bcftools index -f {output.full_vcf} --threads {threads}
        """
