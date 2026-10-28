rule star_alignment:
    input:
        r1=config["input"]["rnaseq"] + "{sample}_1.fastq.gz",
        r2=config["input"]["rnaseq"] + "{sample}_2.fastq.gz",
    output:
        bam=config["output"] + "output_vkaryo/star_aligned/{sample}_Aligned.sortedByCoord.out.bam",
        chimeric=config["output"] + "output_vkaryo/star_aligned/{sample}_Chimeric.out.junction",
        # NOTA: --quantMode GeneCounts produz esse arquivo como efeito
        # colateral -- precisa ser declarado como output para o Snakemake
        # saber que esta regra o produz (senão alvos que dependem dele, como
        # validate_cnv_expression, só funcionam se o arquivo já existir "por
        # acaso" de uma execução anterior).
        reads_per_gene=config["output"] + "output_vkaryo/star_aligned/{sample}_ReadsPerGene.out.tab",
    params:
        sample_prefix=config["output"] + "output_vkaryo/star_aligned/{sample}_",
        index=config["refgenome"]["star_index"],
    benchmark: config["output"] + "output_vkaryo/benchmarks/star_align/star_align_{sample}_benchmark.txt",
    conda:
        "../envs/star_align.yaml" 
    threads: 15
    shell:
        """
        echo "##############################################"
        echo "--------  Running STAR-ALiGNMENT ------------"
        echo "##############################################"  
        STAR --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.sample_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --chimOutType Junctions WithinBAM \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --chimScoreMin 1 \
            --chimScoreDropMax 30 \
            --chimScoreSeparation 5 \
            --chimMultimapNmax 10 \
            --chimNonchimScoreDropMin 10\
            --twopassMode Basic


        """
