rule manta:
    input:
        bam=config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        bai=config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bai",
    output:
        f1=config["output"] + "output_vkaryo/manta/{sample}/results/variants/candidateSV.vcf.gz",
    params:
        ref=config["refgenome"]["ref"],
        output_dir = config["output"] + "output_vkaryo/manta/{sample}",
    conda:
        "../envs/manta.yaml"
    threads: 5
    benchmark: config["output"] + "output_vkaryo/benchmarks/manta/manta_{sample}_benchmark.txt",
    shell:
        """

        # TO check if dir alredy exist and exclude it to avoid error with manta
        if [ -d "{params.output_dir}" ]; then
            rm -rf "{params.output_dir}"
        fi

        echo "##############################################" 
        echo "-----------    Running MANTA      ------------" 
        echo "##############################################"
        echo "Verificando arquivos de entrada"
        
        ls -lh {input.bam}
        ls -lh {params.ref}
        # Criação do diretório de execução do Manta
        mkdir -p {params.output_dir}        
        configManta.py --bam {input.bam} --referenceFasta \
        {params.ref} --runDir {params.output_dir} --exome

        # Executa o Manta no diretório gerado
        ( cd {params.output_dir} && ./runWorkflow.py --mode local -j {threads})

        """
