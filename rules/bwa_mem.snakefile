rule bwa_mem:
    input:
       r1=config["output"] + "output_vkaryo/fastp/exoma/{sample}_fastp_1.fastq.gz",
       r2=config["output"] + "output_vkaryo/fastp/exoma/{sample}_fastp_2.fastq.gz",
    output:
        bam=config["output"] + "output_vkaryo/aligned/{sample}_aligned.bam",
    conda:
        "../envs/bwa.yaml"
    threads: 13

    params:
        threads=13,
        ref=config["refgenome"]["ref"],
        dir=config["output"] + "output_vkaryo/aligned/",
        temp=config["output"] + "output_vkaryo/aligned/temp/",
        out=config["output"] + "output_vkaryo/aligned/{sample}_aligned.bam",
        sample_name="{sample}",
        read_group="'@RG\\tID:{sample}\\tSM:{sample}\\tLB:SureSelectV6\\tPL:ILLUMINA'"

    benchmark: config["output"] + "output_vkaryo/benchmarks/aligned/bwa_{sample}_benchmark.txt",

    shell:
        """
        echo "##############################################" 
        echo "-----------    Running BWA MEM    ------------" 
        echo "##############################################"

        ## running bwa
        mkdir -p {params.dir}
        mkdir -p {params.temp}

        bwa mem \
        -M \
        -t {params.threads} \
        -v 3 \
        -R {params.read_group} \
        {params.ref} \
        {input.r1} {input.r2} | samtools sort --threads {params.threads} -o {params.out} -T {params.temp}

        """