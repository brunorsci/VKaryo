rule arriba_index:
    input:
        # Entradas se houver
    output:
        blacklist="reference/blacklist_hg38_GRCh38_v2.5.0.tsv.gz",
        fusions_ref="reference/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz",
        protein_domains="reference/protein_domains_hg38_GRCh38_v2.5.0.gff3",
    params:
        genoma_version="GRCh38+RefSeq",
        out_dir="reference/"
    conda:
        "../envs/arriba.yaml"
    threads: 10
    shell:
        """
        if [ ! -f {output.blacklist} ] || [ ! -f {output.fusions_ref} ] || [ ! -f {output.protein_domains} ]; then
            echo "##############################################" 
            echo "-------   Running ARRIBA INDEX   ------------" 
            echo "##############################################"
            scripts/download_references.sh {params.genoma_version} {params.out_dir}
        else
            echo "ARRIBA index files already exist. Skipping rule."
        fi
        """
