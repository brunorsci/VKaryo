#!/bin/bash

# Variáveis para armazenar os valores das flags
input_vcf=""
output_vcf=""

zcat $input_vcf > input_vcf


while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) input_vcf="$2"; shift ;; 
        -o|--output) output_vcf="$2"; shift ;;
        *) echo "Flag desconhecida: $1" ; exit 1 ;;
    esac
    shift
done

if [[ -z "$input_vcf" || -z "$output_vcf" ]]; then
    echo "Uso: $0 -i <input_vcf> -o <output_vcf>"
    exit 1
fi

# Descompactar o arquivo .vcf.gz
if [[ "$input_vcf" == *.vcf.gz ]]; then
    temp_vcf="${input_vcf%.gz}"  # Remove a extensão .gz para criar um nome temporário
    zcat "$input_vcf" > "$temp_vcf"
    input_vcf="$temp_vcf"  # Atualiza a variável input_vcf para o arquivo descompactado
fi

# Cabeçalho VCF
echo -e "##fileformat=VCFv4.2" >> "$output_vcf"
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> "$output_vcf"

# Processamento do arquivo
while IFS=$'\t' read -r chrom pos id ref alt qual filter info format sample; do
    # Ajuste dos campos INFO
    info_keys=$(echo "$info" | tr ';' '\n' | cut -d '=' -f 1 | paste -sd ':' -)
    info_values=$(echo "$info" | tr ';' '\n' | cut -d '=' -f 2 | paste -sd ':' -)
    
    # Corrige a linha para BND
    if [[ $id == *"BND"* ]]; then
        format="SVTYPE:MATEID:CIPOS:HOMLEN:HOMSEQ:BND_PAIR_COUNT:PAIR_COUNT"
        sample="BND:$id:$pos:2:AG:2:2"
    fi

    # Gera a linha final
    echo -e "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info_keys\t$format\t$sample" >> "$output_vcf"
done < <(grep -v '^#' "$input_vcf")
echo "File: $input_vcf was processed and $output_vcf was generated!"
