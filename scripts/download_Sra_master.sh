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

# Função para tentar baixar com fasterq-dump
download_fasterq() {
    local id=$1
    echo "Tentando baixar com fasterq-dump a amostra $id..."
    
    # Tenta baixar diretamente o arquivo com fasterq-dump
    fasterq-dump --progress --threads 56 --mem 50000 --split-3 --force "$id" || return 1
}

# Função para tentar baixar com prefetch
download_prefetch() {
    local id=$1
    echo "Tentando baixar com prefetch a amostra $id..."
    # Tenta baixar o arquivo, e verifica se já existe arquivo de lock
    prefetch "$id" || {
        echo "Erro ao baixar com prefetch. Tentando resolver bloqueio e tentar novamente."
        rm -f "/home/patgen/hd3/cariotipagem/painel100LL/rnaseq/$id/$id.sralite.lock" 2>/dev/null
        prefetch "$id" || return 1
    }
}

# Itera sobre cada ID no arquivo
while read -r id; do
    # Verifica se o ID já está na lista de amostras baixadas
    if [[ "$baixadas_str" =~ " $id " ]]; then
        echo "Amostra $id já baixada. Pulando..."
    else
        echo "Baixando amostra $id..."

        # Primeiro tenta baixar com fasterq-dump
        if download_fasterq "$id"; then
            echo "Amostra $id baixada com sucesso via fasterq-dump!"
        else
            echo "Erro ao baixar com fasterq-dump, tentando com prefetch..."

            # Caso o fasterq-dump falhe, tenta baixar com prefetch
            if download_prefetch "$id"; then
                echo "Amostra $id baixada com sucesso via prefetch!"

                # Após o download com prefetch, tenta rodar o fasterq-dump no arquivo .sra
                download_fasterq "$id" && echo "Conversão realizada com sucesso para $id." || echo "Falha na conversão de $id com fasterq-dump."
            else
                echo "Falha ao baixar a amostra $id. Verifique a conexão ou a disponibilidade do arquivo."
            fi
        fi

        # Após o download e a conversão, verifica se os arquivos .fastq existem
        if [[ -f "$id.fastq" ]]; then
            pigz "$id.fastq"
        fi
        if [[ -f "$id"_1.fastq ]]; then
            pigz "$id"_1.fastq
        fi
        if [[ -f "$id"_2.fastq ]]; then
            pigz "$id"_2.fastq
        fi
    fi
done < "$SRA_FILE"

echo "Processo concluído!"
