### DOWNLOAD OF DATABASE (takes a long long time)
vep_install -a cf -s homo_sapiens -y GRCh38

### RUNNING VEP WITH DELLY RESULTS
vep -i COL3_del.bcf  -o COL3_del_annotated.vcf --format vcf --database
vep -i COL3_dup.bcf  -o COL3_dup_annotated.vcf --format vcf --database
vep -i COL3_inv.bcf  -o COL3_inv_annotated.vcf --format vcf --database
vep -i COL3_tra.bcf  -o COL3_tra_annotated.vcf --format vcf --database

### RUNNING WITH CN>MOPS VCF RESULTS
vep -i resultscnmops_cnvs.vcf  -o resultscnmops_cnvs_annotated.vcf --format vcf --database

## HERE WITHOUT THE  FLAG --DATABASE (TO TEST)
vep -i resultscnmops_cnvs.vcf \
    -o resultscnmops_cnvs_annotated.vcf \
    --format vcf \
    --cache \
    --dir_cache /caminho/para/cache_vep \
    --species homo_sapiens \
    --offline \
    --vcf

### ONLINE MODE
vep -i resultscnmops_cnvs.vcf \
    -o resultscnmops_cnvs_annotated.vcf \
    --format vcf \
    --species homo_sapiens \
    --vcf