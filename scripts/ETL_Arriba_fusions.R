#### ETL ARRIBA FUSIONS ####
### ARRIBA FUSIONS (ok)___####
library(dplyr)
library(janitor)
library(tidyverse)

# Lê o argumento da linha de comando para o diretório base
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Erro: Nenhum diretório de entrada foi fornecido como argumento.")
}
base_dir = args[1]
base_dir = gsub("/$", "", base_dir)  

out_dir = args[2]# diretorio de output
out_dir= "/PATH/TO/results/output_vkaryo/etl_arriba/"

## DIRETORIO
#setwd("~/cariotipagem/results/output_vkaryo/arriba/")

# Lista os arquivos .tsv.gz
tsv_files <- list.files(path=base_dir,pattern = "_fusions.tsv$")

# Loop para ler, adicionar coluna Samples_ID e forçar tudo para character
for (i in seq_along(tsv_files)) {
  sample_name <- gsub("_fusions.tsv$", "", tsv_files[i])
  
  assign(
    x = sample_name,  
    value = fread(tsv_files[i]) %>%
      mutate(Samples_ID = sample_name) %>%
      mutate(across(everything(), as.character))
  )
}

rm(i, tsv_files, sample_name)

Arriba_fusions = eval(parse(text = 
                              paste0("bind_rows(",paste0(ls(pattern = "^ERR"), collapse = ","),")")
))
rm(list = ls(pattern = "^ERR"))

#
library(data.table)

# Remover '#' do início de qualquer nome de coluna
setnames(Arriba_fusions, gsub("^#", "", names(Arriba_fusions)))

# Transformação
Arriba_clean <- Arriba_fusions %>%
  # Remover '#' de colunas como "#gene1"
  { setNames(., gsub("^#", "", names(.))) } %>%
  # Padronizar nomes
  clean_names() %>%
  # Criar coluna fusion
  mutate(fusion = paste0(toupper(gene1), "--", toupper(gene2))) %>%
  # Reordenar para ter as colunas desejadas primeiro
  relocate(all_of( c("samples_id", "fusion","type", "confidence", "gene1", "gene2", 
                      "reading_frame")), .before = 1) %>%
  # Remover duplicatas
  distinct() %>%
  # Filtrar baixa confiança
  filter(confidence != "low")

Arriba_etl <- Arriba_clean %>%
  dplyr::select(samples_id, fusion,type, confidence, 
                retained_protein_domains, gene_id1, gene_id2, reading_frame,
                breakpoint1,breakpoint2, site1, site2, coverage1, coverage2,
                direction1, direction2,fusion_transcript,peptide_sequence)

# Salvar
#write.csv(Arriba_clean, "arriba_fusions_clean.csv", row.names = FALSE)

# Exportar
fwrite(Arriba_clean, paste0(out_dir,"arriba_fusions_all.csv"))
fwrite(Arriba_etl, paste0(out_dir,"arriba_etlclean_all.csv"))
