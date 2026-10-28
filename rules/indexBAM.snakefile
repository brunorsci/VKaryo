rule indexbam:
    input:
        config["output"] + "output_vkaryo/marked/{sample}_markduplicates.bam"
    output:
        config["output"] + "output_vkaryo/marked/{sample}_markduplicates.bai"
    conda:
        "../envs/picard.yaml"
    threads: 5
    shell:
        """
        echo "##############################################" 
        echo "-----    Running PICARD BAM Index    ---------" 
        echo "##############################################"

        picard BuildBamIndex  --INPUT {input} --OUTPUT {output}
        """