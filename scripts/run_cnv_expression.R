#!/usr/bin/env Rscript
# Roda o pipeline completo do pacote CNVExpression (versão instalada neste
# servidor, testada pelo usuário -- github.com/brunorsci/CNVExpression_v1,
# v0.99.0) sobre os inputs concatenados de todas as amostras DNA+RNA
# pareadas (ver scripts/prepare_cnv_expression_input.py).
#
# NOTA IMPORTANTE (12/08/2026): create_gene_granges(id_type="symbol") do
# pacote tem um bug no caminho biomaRt -- SEMPRE nomeia o GRanges resultante
# pelo Ensembl gene_id, mesmo quando id_type="symbol" é pedido (o caminho
# offline via TxDb faz certo, o caminho biomaRt não, porque getBM() nem
# devolve o símbolo de volta quando filtra por external_gene_name -- só
# ensembl_gene_id). Resultado: create_expression_se() nunca encontra genes
# em comum entre a matriz de expressão (símbolos, ex. TP53) e o GRanges
# (nomeado por ENSG...), e o pipeline falha com "No common genes found".
# Por instrução explícita do usuário, NÃO modificamos o pacote instalado
# (é a versão testada e em uso) -- em vez disso, reimplementamos só essa
# etapa aqui (consulta biomaRt direta, nomeando por símbolo) e chamamos as
# demais funções do pacote (todas as outras, sem bug) na mesma ordem que
# cnv_expression_pipeline() usaria internamente.
#
# Uso:
#   Rscript run_cnv_expression.R <cnv.tsv> <expression.tsv> <out_filtered.tsv> <out_merged.tsv> <out_raw.tsv> [logfc_cutoff] [pvalue_cutoff]

suppressMessages(library(CNVExpression))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Uso: run_cnv_expression.R <cnv.tsv> <expression.tsv> <out_filtered.tsv> <out_merged.tsv> <out_raw.tsv> [logfc_cutoff] [pvalue_cutoff]")
}
cnv_file <- args[1]
expr_file <- args[2]
out_filtered <- args[3]
out_merged <- args[4]
out_raw <- args[5]
logfc_cutoff <- if (length(args) >= 6) as.numeric(args[6]) else 1
pvalue_cutoff <- if (length(args) >= 7) as.numeric(args[7]) else 0.05

# ---- workaround pro bug de create_gene_granges(id_type="symbol") ----
create_gene_granges_by_symbol <- function(gene_symbols, genome = "GRCh38") {
    standard_chroms <- c(as.character(1:22), "X", "Y")
    if (requireNamespace("biomaRt", quietly = TRUE)) {
        message("  (consultando biomaRt por símbolo -- workaround do bug id_type=symbol)")
        ensembl <- if (genome == "GRCh37") {
            biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
        } else {
            biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
        }
        gene_info <- biomaRt::getBM(
            attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
            filters = "hgnc_symbol", values = gene_symbols, mart = ensembl
        )
        gene_info <- gene_info[gene_info$chromosome_name %in% standard_chroms & gene_info$hgnc_symbol != "", ]
        gene_info <- gene_info[!duplicated(gene_info$hgnc_symbol), ]
        genes_gr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", gene_info$chromosome_name),
            ranges = IRanges::IRanges(start = gene_info$start_position, end = gene_info$end_position),
            strand = "*", gene_id = gene_info$hgnc_symbol
        )
        names(genes_gr) <- genes_gr$gene_id
        return(genes_gr)
    }
    # fallback offline (TxDb) -- esse caminho do pacote já respeita id_type
    # corretamente, então pode ser usado direto sem workaround.
    message("  biomaRt indisponível, usando create_gene_granges() offline (TxDb) do pacote")
    CNVExpression::create_gene_granges(gene_symbols, id_type = "symbol", genome = genome)
}

message("Carregando CNV: ", cnv_file)
cnv_data <- load_cnv_data(cnv_file)

message("Carregando expressão: ", expr_file)
# Colunas do arquivo se chamam "Ensembl_ID" (rótulo padrão esperado por
# load_expression_data()), mas os valores são SÍMBOLOS de gene -- o GTF
# usado pelo STAR neste pipeline foi indexado com gene_name, não gene_id.
expr_data <- load_expression_data(expr_file, gene_col = "Ensembl_ID",
                                   remove_par_y = FALSE, remove_version = FALSE)

message("Amostras: ", ncol(expr_data), " | Genes: ", nrow(expr_data),
        " | Segmentos de CNV: ", nrow(cnv_data))

result <- tryCatch({
    message("Etapa 1/6: preparando dados de CNV...")
    cnv_processed <- prepare_cnv_data(cnv_data, "segment_mean")
    cnv_processed$sample <- standardize_sample_ids(cnv_processed$sample)
    colnames(expr_data) <- standardize_sample_ids(colnames(expr_data))

    message("Etapa 2/6: criando GRangesList de CNV...")
    cnv_grl <- create_cnv_grangeslist(cnv_processed, sample_col = "sample")

    message("Etapa 3/6: calculando regiões de CNV populacionais...")
    cnvrs <- calculate_population_ranges(cnv_grl, verbose = TRUE)

    message("Etapa 4/6: criando GRanges de genes (por símbolo -- ver nota do workaround acima)...")
    genes_gr <- create_gene_granges_by_symbol(rownames(expr_data), genome = "GRCh38")
    message("  ", length(genes_gr), " genes com coordenadas resolvidas")

    message("Etapa 5/6: criando SummarizedExperiment...")
    rse <- create_expression_se(expr_data, genes_gr)

    message("Etapa 6/6: rodando análise CNV-expressão (CNVRanger::cnvEQTL)...")
    results_raw <- run_cnv_expression_analysis(cnvrs, cnv_grl, rse, verbose = TRUE)

    results_filtered <- filter_cnv_expression_results(
        results_raw, logfc_cutoff = logfc_cutoff, pvalue_cutoff = pvalue_cutoff
    )
    results_merged <- process_integration_results(results_filtered, cnv_processed, sample_col = "sample")

    # add_gene_symbols() do pacote exige uma coluna gene_id (Ensembl/Entrez)
    # que os resultados de cnvEQTL não têm -- em vez disso, casa as
    # coordenadas (seqnames/start/end) de cada resultado de volta contra
    # genes_gr (que já está nomeado por símbolo, graças ao workaround acima)
    # e anexa o símbolo do gene em cada linha, pra facilitar a leitura sem
    # precisar decorar coordenadas genômicas.
    annotate_gene_symbol <- function(df) {
        if (!all(c("seqnames", "start", "end") %in% colnames(df)) || nrow(df) == 0) return(df)
        result_gr <- GenomicRanges::GRanges(
            seqnames = df$seqnames,
            ranges = IRanges::IRanges(start = df$start, end = df$end)
        )
        hits <- GenomicRanges::findOverlaps(result_gr, genes_gr, select = "first")
        df$gene_symbol <- ifelse(is.na(hits), NA_character_, names(genes_gr)[hits])
        df
    }
    results_raw <- annotate_gene_symbol(results_raw)
    results_filtered <- annotate_gene_symbol(results_filtered)
    results_merged <- annotate_gene_symbol(results_merged)

    list(results_raw = results_raw, results_filtered = results_filtered, results_merged = results_merged)
}, error = function(e) {
    message("ERRO no pipeline: ", conditionMessage(e))
    quit(status = 1)
})

message("Salvando resultados...")
write.table(result$results_raw, out_raw, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result$results_filtered, out_filtered, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result$results_merged, out_merged, sep = "\t", quote = FALSE, row.names = FALSE)

message("Concluído: ", nrow(result$results_raw), " associações testadas, ",
        nrow(result$results_filtered), " significativas (logFC >= ", logfc_cutoff,
        ", padj <= ", pvalue_cutoff, ")")
