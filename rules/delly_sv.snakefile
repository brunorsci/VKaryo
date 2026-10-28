rule delly_sv:
    input:
        bam=config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bam",
        bai=config["output"] + "output_vkaryo/recal_apply/{sample}_markduplicates_recal.bai",
    output:
        sv_bcf=config["output"] + "output_vkaryo/delly/svs/{sample}_sv.bcf",
    params:
        reference=config["refgenome"]["ref"],
    threads: 2
    conda:
        "../envs/delly.yaml"
    shell:
        """   
        echo "##############################################"
        echo "-------    Running   DELLY SV      ----------"
        echo "##############################################"  

        delly call -g {params.reference} -o {output.sv_bcf} {input.bam}
        
        """
