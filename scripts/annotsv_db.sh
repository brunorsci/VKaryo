#!/bin/bash

# Definir o diretório onde os bancos de dados serão armazenados
DB_DIR="./annotSV_db"
mkdir -p $DB_DIR

# Baixar os bancos de dados necessários
echo "Baixando os bancos de dados..."

# Banco de dados dbNSFP
#wget -P $DB_DIR https://jhi.sdsu.edu/dbNSFP/dbNSFPv4.0a.zip

# Banco de dados COSMIC
wget -P $DB_DIR https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic_v93.zip

# Banco de dados Decipher
wget -P $DB_DIR https://decipher.sanger.ac.uk/ftp/decipher-current.txt

# Banco de dados ClinVar (opcional, se necessário)
wget -P $DB_DIR ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20230404.vcf.gz

# Extrair os arquivos .zip
echo "Extraindo arquivos..."

unzip $DB_DIR/dbNSFPv4.0a.zip -d $DB_DIR
unzip $DB_DIR/cosmic_v93.zip -d $DB_DIR

# Caso haja outros formatos a serem extraídos, inclua aqui.
# Exemplo de extração de VCF (se necessário)
gunzip $DB_DIR/clinvar_20230404.vcf.gz

# Organizar os arquivos no formato esperado pelo AnnotSV
echo "Organizando os arquivos para o AnnotSV..."

# Criar diretórios adicionais se necessário
mkdir -p $DB_DIR/dbNSFP
mkdir -p $DB_DIR/cosmic

# Mover os arquivos extraídos para as pastas de destino
mv $DB_DIR/dbNSFPv4.0a.txt $DB_DIR/dbNSFP/
mv $DB_DIR/cosmic_v93.txt $DB_DIR/cosmic/

echo "Bancos de dados prontos para o AnnotSV."
