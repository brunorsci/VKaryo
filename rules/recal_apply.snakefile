rule recal_apply:
    input:
        bam=config["output"] + "output_vkaryo/marked/{sample}_markduplicates.bam",
        recal_table=config["output"] + "output_vkaryo/recal_table/{sample}_recal.grp",
    output:
        bam=config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        bai=config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bai",
    params: ref=config["refgenome"]["ref"],
    conda:
        "../envs/gatk.yaml"
    threads: 4
    benchmark: config["output"] + "output_vkaryo/benchmarks/recal_apply/{sample}_benchmark.txt",

    shell:
        """
        echo "##############################################" 
        echo "-------    Running GATK AplyBQSR    ---------" 
        echo "##############################################"

        gatk ApplyBQSR \
        -I {input.bam} \
        -R {params.ref} \
        --bqsr-recal-file {input.recal_table} \
        -O {output.bam}
        """