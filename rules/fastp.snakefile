rule fastp:
    input:
        "samples/{wildcards.sample}.fastq"
    output:
        "results/trimmed/{sample}_fastp.fastq",
        "results/QC/{sample}_fastp.html",
        "results/QC/{sample}_fastp.json"
    shell:
        """
        fastp -i {input} -o {output[0]} --html {output[1]} --json {output[2]}
        """
