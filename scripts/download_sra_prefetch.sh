#!/bin/bash

# Arquivo com os IDs dos runs
SRA_FILE="sra_runs.txt"

# Verifica se o arquivo existe
if [[ ! -f "$SRA_FILE" ]]; then
    echo "Erro: Arquivo $SRA_FILE não encontrado!"
    exit 1
fi

# Obtém a lista de amostras já baixadas (remove extensões corretamente)
amostras_baixadas=()
for file in *.fastq.gz *.sra 2>/dev/null; do
    [[ -f "$file" ]] && amostras_baixadas+=("$(basename "$file" | sed -E 's/(\.fastq\.gz|\.sra)$//')")
done

# Converte a lista para uma string separada por espaços para facilitar a comparação
baixadas_str=" ${amostras_baixadas[*]} "

# Itera sobre cada ID no arquivo
while read -r id; do
    # Verifica se o ID já está na lista de amostras baixadas
    if [[ "$baixadas_str" =~ " $id " ]]; then
        echo "Amostra $id já baixada. Pulando..."
    else
        echo "Baixando amostra $id..."
        
        # Baixa o arquivo SRA primeiro com prefetch
        prefetch "$id"
        
        # Verifica se o arquivo foi baixado corretamente
        if [[ -f "$HOME/ncbi/public/sra/$id.sra" ]]; then
            # Converte para FASTQ
            fasterq-dump --threads 56 --mem 50000 --split-3 --force "$HOME/ncbi/public/sra/$id.sra"
            
            # Comprime os arquivos gerados, verificando se existem antes de compactar
            [[ -f "$id.fastq" ]] && pigz "$id.fastq"
            [[ -f "$id"_1.fastq ]] && pigz "$id"_1.fastq
            [[ -f "$id"_2.fastq ]] && pigz "$id"_2.fastq
            
            echo "Amostra $id processada com sucesso."
        else
            echo "Falha ao baixar a amostra $id: Verifique a conexão ou a disponibilidade do arquivo."
        fi
    fi
done < "$SRA_FILE"

echo "Processo concluído!"
