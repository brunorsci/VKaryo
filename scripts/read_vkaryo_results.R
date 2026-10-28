
library(tidyverse)
library(data.table)
library('vcfppR')


### ARRIBA FUSIONS (ok)___####
setwd("~/cariotipagem/results/output_vkaryo/arriba/")

# Lista os arquivos .tsv.gz
tsv_files <- list.files(pattern = "_fusions.tsv$")

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

arriba_fusions = eval(parse(text = 
                              paste0("bind_rows(",paste0(ls(pattern = "^ERR"), collapse = ","),")")
))
rm(list = ls(pattern = "^ERR"))

# Criar a coluna 'fusao' unindo '#gene1' e 'gene2'
arriba_fusions <- arriba_fusions %>%
  mutate(fusao = paste(`#gene1`, gene2, sep = "-"))

#### ANNOTSV(delly and manta) ---- ####
setwd("~/cariotipagem/results/output_vkaryo/annotsv/delly/")

# Lista os arquivos .tsv.gz
tsv_files <- list.files(pattern = "_annotsv_delly.tsv$")
# Loop para ler, adicionar coluna Samples_ID e forçar tudo para character
for (i in seq_along(tsv_files)) {
  sample_name <- gsub("_annotsv_delly.tsv$", "", tsv_files[i])
  
  assign(
    x = sample_name,  # nome do objeto salvo no R
    value = fread(tsv_files[i]) %>%
      mutate(Samples_ID = sample_name) %>%
      mutate(across(everything(), as.character))
    
  )
}

rm(i)

Annotsv_delly = eval(parse(text = paste0("bind_rows(",paste0(ls(pattern = "^ERR"), collapse = ","),")")
))
rm(list = ls(pattern = "^ERR"))

#Annotsv_delly <- Annotsv_delly %>%
#  mutate(SV_ALT = paste(SV_chrom, ALT, sep = "_"))


#Annotsv_delly=Annotsv_delly %>% 
#  dplyr::select(Samples_ID,SV_ALT, SV_chrom,ALT,CytoBand, 
#         Location, SV_type, Gene_name,GenCC_disease, OMIM_phenotype,P_gain_phen)

## MANTA
setwd("~/cariotipagem/results/output_vkaryo/annotsv/manta/")

# Lista os arquivos .tsv.gz
tsv_files <- list.files(pattern = "_annotsv_manta.tsv$")
# Loop para ler, adicionar coluna Samples_ID e forçar tudo para character
for (i in seq_along(tsv_files)) {
  sample_name <- gsub("_annotsv_manta.tsv$", "", tsv_files[i])
  
  assign(
    x = sample_name,  # nome do objeto salvo no R
    value = fread(tsv_files[i]) %>%
      mutate(Samples_ID = sample_name) %>%
      mutate(across(everything(), as.character))
    
  )
}

rm(i)

Annotsv_manta = eval(parse(text = paste0("bind_rows(",paste0(ls(pattern = "^ERR"), collapse = ","),")")
))
rm(list = ls(pattern = "^ERR"))


#### ANALYSING ######################
SV_leukemia_delly = Annotsv_delly %>%
  janitor::clean_names() %>% 
  filter(str_detect(gen_cc_disease, "leukemia") == TRUE) %>%
  select(annot_sv_id, samples_id, id,  sv_type,
         sv_chrom, sv_start, sv_end, sv_length,
         cyto_band, gene_name, gene_count,
         frameshift, location, location2,
         gnom_ad_p_li,
         gen_cc_disease, gen_cc_classification, gen_cc_pmid,
         omim_id, omim_phenotype, omim_inheritance, 
         acmg, acmg_class,
         starts_with("p_"),
         info, format)

SV_leukemia_manta = Annotsv_manta %>%
  janitor::clean_names() %>% 
  filter(str_detect(gen_cc_disease, "leukemia") == TRUE) %>%
  select(annot_sv_id, samples_id, id, sv_type,
         sv_chrom, sv_start, sv_end, sv_length,
         cyto_band, gene_name, gene_count,
         frameshift, location, location2,
         gnom_ad_p_li,
         gen_cc_disease, gen_cc_classification, gen_cc_pmid,
         omim_id, omim_phenotype, omim_inheritance, 
         acmg, acmg_class,
         starts_with("p_"),
         info, format)
  

# SV COMUNS ENTRE DELLY E MANTA
SV_leukemia_communs <- inner_join(SV_leukemia_delly, SV_leukemia_manta,
                        by = c("gene_name", "samples_id", "sv_type"),
                        suffix = c("_delly", "_manta"))

key_cols <- c("gene_name", "samples_id", "sv_type")

SV_leukemia_communs <- inner_join(SV_leukemia_delly, SV_leukemia_manta, by = key_cols, suffix = c(".delly", ".manta")) %>%
  select(all_of(key_cols), ends_with(".delly")) %>%
  rename_with(~ sub("\\.delly$", "", .x), ends_with(".delly"))

# SALVAR
save.image("~/cariotipagem/results/output_vkaryo/resultados.RData")

#annotsv %>%
#  filter(str_detect(gen_cc_disease, "cholestasis") == TRUE) %>% 
#  write.table(., file = "anotsv_pass_leukemia.tsv", sep = "\t", 
#              row.names = FALSE, quote = FALSE, col.names = TRUE)


