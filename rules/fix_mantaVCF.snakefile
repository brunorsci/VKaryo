rule fix_vcf_format:
    input:
        vcf=config["output"] + "output_vkaryo/manta/{sample}/results/variants/candidateSV.vcf.gz"
    output:
        vcf_fixed=config["output"] + "output_vkaryo/manta/{sample}/{sample}_fixed.vcf",
    threads: 2

    shell:
        """
        echo "##############################################" 
        echo "-------  Fixing Manta VCF FORMAT   -----------" 
        echo "##############################################"

        scripts/fix_vcf_format.sh -i {input.vcf} -o {output.vcf_fixed}
        
        """
