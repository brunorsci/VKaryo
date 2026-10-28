#!/bin/bash
projects=("PRJNA1309308" "PRJNA1269993" "PRJNA1282768")
#for projeto in $(cat bioproject.txt); do
for projeto in "${projects[@]}"; do
    mkdir -p "$projeto"
    cd "$projeto" || exit 1
    
    # Obtém os IDs dos runs do projeto
    esearch -db sra -query "$projeto"'[bioproject]' | esummary > access_token.txt
    grep -oP 'Run acc="\K[^"]+' access_token.txt > sra_runs.txt
    
    # Verifica quais amostras já foram baixadas
    amostras_baixadas=($(ls *.fastq.gz *.fastq.gz.md5 *.sra 2>/dev/null | sed -E 's/(.fastq.gz|.sra)//g'))
    
    # Itera sobre cada ID de run
    while read -r id; do
        # Verifica se a amostra já foi baixada
        if [[ " ${amostras_baixadas[@]} " =~ " $id " ]]; then
            echo "Amostra $id já baixada. Pulando..."
        else
            # Baixa e comprime os arquivos FASTQ
            fasterq-dump --progress --threads 56 --mem 50000 "$id"
            pigz "$id.fastq"
            pigz "$id"_1.fastq
            pigz "$id"_2.fastq
        fi
    done < sra_runs.txt
    
    cd .. || exit 1
done

