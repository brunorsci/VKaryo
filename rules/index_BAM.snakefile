rule indexbam_STAR:
    input:
        bam=config["output"] + "output_vkaryo/star_aligned/{sample}_Aligned.sortedByCoord.out.bam",
        control=config["output"] + "output_vkaryo/star_aligned/{sample}_Chimeric.out.junction",
    output:
        config["output"] + "output_vkaryo/star_aligned/{sample}_Aligned.sortedByCoord.out.bai" 
    threads: 5
    conda:
        "../envs/picard.yaml" 
    shell:
        """
        echo "##############################################" 
        echo "-----Running PICARD BAM Index for STAR--------" 
        echo "##############################################"

        picard BuildBamIndex  --INPUT {input.bam} --OUTPUT {output}
        """