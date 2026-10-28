library(data.table)
library(dplyr)
library(janitor)

# Diretório base onde estão os arquivos *_fusions.tsv
base_dir <- "/PATH/TO/results/output_vkaryo/arriba/"
output_dir <- "/PATH/TO/results/output_vkaryo/etl_arriba/"

# Diretório para salvar os ETLs individuais (pode ser o mesmo ou outro)
#output_dir <- file.path(base_dir, "etl_output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Listar os arquivos no diretório base
arquivos <- list.files(path = base_dir, pattern = "_fusions.tsv$", full.names = TRUE)

for (arq in arquivos) {
  
  # Ler arquivo
  df <- fread(arq)
  
  # Nome da amostra (removendo o sufixo do arquivo e o caminho)
  sample_id <- sub("_fusions.tsv$", "", basename(arq))
  
  # ETL aplicado ao arquivo individual
  df_etl <- df %>%
    # Remover '#' de colunas
    { setNames(., gsub("^#", "", names(.))) } %>%
    # Padronizar nomes
    clean_names() %>%
    # Criar coluna fusion
    mutate(
      fusion = paste0(toupper(gene1), "--", toupper(gene2)),
      samples_id = sample_id
    ) %>%
    # Reordenar para ter as colunas desejadas primeiro
    relocate(all_of( c("samples_id", "fusion","type", "confidence", "gene1", "gene2", 
                       "reading_frame")), .before = 1) %>%
    # Remover duplicatas
    distinct() %>%
    # Filtrar baixa confiança
    filter(confidence != "low")

  
  # Salvar ETL individual para essa amostra no diretório de saída
  fwrite(df_etl, file.path(output_dir, paste0(sample_id, "_fusion_etl.csv")))
  
}
