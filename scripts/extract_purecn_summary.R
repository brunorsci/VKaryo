#!/usr/bin/env Rscript
# Extrai sumário de pureza e ploidy do PureCN

library(PureCN)

args <- commandArgs(trailingOnly = TRUE)
rds_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Carrega resultado PureCN
ret <- readRDS(rds_file)

# Extrai informações
purity <- ret$results[[1]]$purity
ploidy <- ret$results[[1]]$ploidy
sex <- ret$input$sex
contamination <- ret$results[[1]]$contamination

# LOH events
loh_events <- sum(ret$results[[1]]$CNV$type == "LOH", na.rm = TRUE)

# CNV summary
cnv <- ret$results[[1]]$CNV
n_amp <- sum(cnv$C >= 3, na.rm = TRUE)  # Amplificações
n_del <- sum(cnv$C <= 1, na.rm = TRUE)  # Deleções
n_loh <- sum(cnv$type == "LOH", na.rm = TRUE)

# Escreve sumário
sink(output_file)
cat("=== PureCN Summary ===\n\n")
cat(sprintf("Purity: %.2f%%\n", purity * 100))
cat(sprintf("Ploidy: %.2f\n", ploidy))
cat(sprintf("Sex: %s\n", sex))
cat(sprintf("Contamination: %.4f\n", contamination))
cat("\n")
cat("=== CNV Summary ===\n")
cat(sprintf("Amplifications (C >= 3): %d\n", n_amp))
cat(sprintf("Deletions (C <= 1): %d\n", n_del))
cat(sprintf("LOH events: %d\n", n_loh))
cat("\n")
cat("=== Quality Metrics ===\n")
cat(sprintf("Solution: %d\n", 1))
cat(sprintf("Log-likelihood: %.2f\n", ret$results[[1]]$log.likelihood))
sink()

cat(sprintf("✓ Summary saved to %s\n", output_file))
