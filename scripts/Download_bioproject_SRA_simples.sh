#!/bin/bash

# Projeto: PRJNA1309308
PROJECT_ID="PRJNA1309308"

# Pega a lista de amostras no ENA (com links para fastq.gz)
wget -q "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT_ID}&result=read_run&fields=run_accession,fastq_ftp&download=true&limit=0" -O ${PROJECT_ID}_runs.txt

# Baixa cada fastq.gz e gzipa se necessário
tail -n +2 ${PROJECT_ID}_runs.txt | while IFS=$'\t' read -r RUN LINKS
do
    echo "==> Baixando $RUN ..."
    for LINK in $(echo $LINKS | tr ";" "\n")
    do
        FILE=$(basename $LINK)
        wget -c ftp://$LINK -O ${FILE}
        
        # Gzipa se não vier comprimido (mas geralmente já vem .fastq.gz)
        if [[ ${FILE} != *.gz ]]; then
            echo "Comprimindo ${FILE} ..."
            gzip ${FILE}
        fi
    done
done
