import pandas as pd
import argparse

# Função para carregar e filtrar o VCF com anotação opcional
def filter_vcf_with_annotation(vcf_file, output_file, annotation_file):
    # Carregar os dados do VCF usando pandas
    df = pd.read_csv(vcf_file, sep='\t', comment='#')

    # Carregar as anotações dos genes relacionados a leucemia do arquivo CSV
    annotations_df = pd.read_csv(annotation_file)

    # Criar um dicionário de anotações a partir do arquivo CSV
    gene_annotations = {
        row['Gene']: [row['Leukemia Type'], row['Classification'], row['Description']] 
        for _, row in annotations_df.iterrows()
    }

    # Lista de genes relacionados a leucemias (exemplo, adicione ou remova conforme necessário)
    leukemia_genes = [
        "TRAF3IP3", "SLC39A13", "SPI1", "RPS19", "ZBTB18", "GATA1", "TAL1", 
        "JAK2", "FLI1", "MLL", "ABL1", "RUNX1", "NPM1", "WT1", "CBL"
    ]

    # Filtro: Alterações associadas a leucemia com base nos genes
    leukemia_related = df[
        df['Gene_name'].str.contains('|'.join(leukemia_genes), case=False, na=False) &  # Genes associados à leucemia
        df['ACMG_class'].str.contains("HI|morbid", case=False, na=False)  # ACMG classifica como "HI" ou "morbid"
    ]

    # Adicionar anotações, se existirem
    leukemia_related['Annotation'] = leukemia_related['Gene_name'].apply(
        lambda gene: gene_annotations.get(gene, None)
    )

    # Salvar o VCF filtrado com anotações
    leukemia_related.to_csv(output_file, sep='\t', index=False)

# Função principal para processar os parâmetros e chamar a função de filtragem
def main():
    parser = argparse.ArgumentParser(description="Filtra VCF com base em genes relacionados à leucemia e anotações.")
    parser.add_argument('-i', '--input', required=True, help="Arquivo de entrada VCF")
    parser.add_argument('-o', '--output', required=True, help="Arquivo de saída VCF filtrado")
    parser.add_argument('-csv', '--annotation_file', required=True, help="Arquivo CSV com as anotações dos genes")

    args = parser.parse_args()

    # Chamar a função para filtrar e adicionar as anotações
    filter_vcf_with_annotation(args.input, args.output, args.annotation_file)

if __name__ == "__main__":
    main()
