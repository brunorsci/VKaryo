#### REFERENCES ####

## hg38 genome reference do UCSC
cd reference
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

#criar os index
bwa index hg38.fa
samtools faidx hg38.fa
picard CreateSequenceDictionary R=hg38.fa O=hg38.dict

## DBSNP (known variants)
# Criar o arquivo rename.txt com os ids dos cromossomos em formato character:
cat <<EOF > rename.txt
NC_000001.11 chr1
NC_000002.12 chr2
NC_000003.12 chr3
NC_000004.12 chr4
NC_000005.10 chr5
NC_000006.12 chr6
NC_000007.14 chr7
NC_000008.11 chr8
NC_000009.12 chr9
NC_000010.11 chr10
NC_000011.10 chr11
NC_000012.12 chr12
NC_000013.11 chr13
NC_000014.9 chr14
NC_000015.10 chr15
NC_000016.10 chr16
NC_000017.11 chr17
NC_000018.10 chr18
NC_000019.10 chr19
NC_000020.11 chr20
NC_000021.9 chr21
NC_000022.11 chr22
NC_000023.11 chrX
NC_000024.10 chrY
NC_012920.1 chrMT
EOF


# Fazer o download do dbsnp 156
wget https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz
bcftools index ./reference/GCF_000001405.40.gz


# Usar bcftools para renomear. Aqui ja modifico o nome do arquivo.
bcftools reheader -f rename.txt -o dbsnp_156_chr.vcf.gz ./reference/GCF_000001405.40.gz
rm -r GCF_000001405.40.gz

## gnomad
wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/sv/gnomad.v4.1.sv.sites.vcf.gz
#gunzip gnomad.v4.1.sv.sites.vcf.gz

## GTF gene annotation 
# GENCODE (recomendado)
# O GENCODE fornece anotações abrangentes e de alta qualidade para o genoma humano. A versão mais recente disponível é a Release 44, que corresponde ao GRCh38.p13
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz

# ENSEMBL
wget ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
# UCSC Genome Browser
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz

### ARRIBA
mkdir -p arriba_data
cd arriba_data


