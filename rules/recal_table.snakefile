rule recal_table:
    input:
        bam=config["output"] + "output_vkaryo/marked/{sample}_markduplicates.bam",
        ref=config["refgenome"]["ref"],
        dict=config["refgenome"]["dict"],
        known_var=config["refgenome"]["known_var"],
    output:
        recal_table=config["output"] + "output_vkaryo/recal_table/{sample}_recal.grp",
    log:
        config["output"] + "output_vkaryo/recal_table/log/{sample}.log",
    conda:
        "../envs/gatk.yaml"
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=25000,
    threads: 4
    benchmark: config["output"] + "output_vkaryo/benchmarks/recal_table/recal_table_{sample}_benchmark.txt",

    shell:
        """
        echo "##############################################" 
        echo "-----   Running GATK BaseRecalibrator  -------" 
        echo "##############################################"

        gatk BaseRecalibrator \
        -I {input.bam} \
        -R {input.ref} \
        --known-sites {input.known_var} \
        -O {output.recal_table}
        """

