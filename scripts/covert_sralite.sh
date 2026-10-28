#!/bin/bash

# Arquivo com os IDs dos runs
SRA_FILE="sra_runs.txt"

# Verifica se o arquivo existe
if [[ ! -f "$SRA_FILE" ]]; then
    echo "Erro: Arquivo $SRA_FILE não encontrado!"
    exit 1
fi

# Função para converter o arquivo .sralite para fastq.gz
convert_sralite() {
    local id=$1
    local dir="/home/patgen/hd3/cariotipagem/painel100LL/rnaseq/$id"

    # Verifica se a pasta do ID existe
    if [[ -d "$dir" ]]; then
        echo "Entrando na pasta $dir..."
        
        # Verifica se o arquivo .sralite existe
        if [[ -f "$dir/$id.sralite" ]]; then
            echo "Arquivo $id.sralite encontrado. Iniciando a conversão para .fastq.gz..."
            
            # Executa a conversão e compressão
            fasterq-dump --progress --threads 56 --mem 50000 --split-3 --force "$dir/$id.sralite" && \
            pigz "$dir/$id.fastq" && echo "Conversão e compressão concluídas para $id!"
        else
            echo "Erro: Arquivo $id.sralite não encontrado em $dir."
        fi
    else
        echo "Erro: A pasta $dir não existe."
    fi
}

# Itera sobre cada ID no arquivo
while read -r id; do
    echo "Processando a amostra $id..."
    convert_sralite "$id"
done < "$SRA_FILE"

echo "Processo concluído!"
