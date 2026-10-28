rule annotsv:
    input:
        vcf=config["output"] + "output_vkaryo/manta/{sample}/{sample}_fixed.vcf"
    output:
        annotated=config["output"] + "output_vkaryo/annotsv/{sample}_annotated.tsv",
    params:
        output_dir=config["output"] + "output_vkaryo/annotsv/",
        db="reference/AnnotSV_DB",  # Diretório contendo bases de dados necessárias
    threads: 60
    conda:
        "../envs/annotsv.yaml"
    shell:
        """
        echo "##############################################" 
        echo "-----------    Running AnnotSV    ------------" 
        echo "##############################################"           


        AnnotSV \
            -SVinputFile {input.vcf} \
            -outputFile {output.annotated} \
            -annotationsDir {params.db}/AnnotSV_annotations  \
            -genomeBuild GRCh38
        """
