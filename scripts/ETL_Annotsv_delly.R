library(dplyr)
library(janitor)
library(stringr)
library(data.table)

Annotsv_clean <- Annotsv_delly[,c(1:120)] %>%
  # Padronizar nomes de colunas
  clean_names() %>%
  # Criar coluna com nomenclatura HGVS
  mutate(
    hgvs_notation = case_when(
      sv_type == "DEL" ~ paste0("chr", sv_chrom, ":g.", sv_start, "_", sv_end, "del"),
      sv_type == "DUP" ~ paste0("chr", sv_chrom, ":g.", sv_start, "_", sv_end, "dup"),
      sv_type == "INS" ~ paste0("chr", sv_chrom, ":g.", sv_start, "_", sv_end, "ins"),
      sv_type == "INV" ~ paste0("chr", sv_chrom, ":g.", sv_start, "_", sv_end, "inv"),
      sv_type %in% c("BND", "TRA") ~ paste0("chr", sv_chrom, ":g.", sv_start, "::chr", sv_end, "bnd"),
      TRUE ~ paste0("chr", sv_chrom, ":g.", sv_start, "_", sv_end, "_", tolower(sv_type))
    )
  ) %>%
  # Reordenar colunas
  relocate(samples_id, hgvs_notation, .before = 1) %>%
  # Filtrar apenas variantes com filtro PASS
  filter(filter == "PASS") %>%
  # Filtrar se QUALQUER uma das colunas contém os termos de interesse
  filter(
    str_detect(tolower(gen_cc_disease), "leukemia|myeloid|lymphoid") |
      str_detect(tolower(p_gain_phen), "leukemia|myeloid|lymphoid") |
      str_detect(tolower(p_loss_phen), "leukemia|myeloid|lymphoid") |
      str_detect(tolower(p_ins_phen), "leukemia|myeloid|lymphoid") |
      str_detect(tolower(p_snvindel_phen), "leukemia|myeloid|lymphoid") |
      str_detect(tolower(omim_phenotype), "leukemia|myeloid|lymphoid")
  ) %>%
  dplyr::select(id, sample=samples_id, hgvs_notation,svnote=annot_sv_id, crhom=sv_chrom,
                start_pos=sv_start, end_pos=sv_end, length=sv_length, type=sv_type,cyto_band, gene_name,
                gene_count,frameshift,exon_count,gen_cc_disease,gen_cc_classification,
                omim_phenotype,omim_inheritance, p_gain_phen,p_loss_phen,
                p_snvindel_phen)


