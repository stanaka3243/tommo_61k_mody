#!/usr/bin/env Rscript
# VEP REST の JSON (chunk_*.json) を 08_merge_vep.R 用の TSV に変換
#   1 変異 x 1 遺伝子 = 1 行 (MANE Select > canonical > 最初の転写産物 の順で選ぶ)
# 使い方: Rscript 07c_parse_vep_rest.R <JSON_DIR> <OUT_TSV>
a <- commandArgs(trailingOnly = TRUE); jdir <- a[1]; out <- a[2]
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  lib <- Sys.getenv("R_LIBS_USER"); dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  install.packages("jsonlite", lib = lib, repos = "https://cloud.r-project.org"); .libPaths(c(lib, .libPaths()))
}
suppressPackageStartupMessages({ library(jsonlite); library(data.table) })

`%||%` <- function(a, b) if (is.null(a)) b else a
g  <- function(x, ...) { for (k in c(...)) { v <- x[[k]]; if (!is.null(v)) return(v) }; NULL }
s1 <- function(x) if (is.null(x) || length(x) == 0) NA_character_ else paste(unlist(x), collapse = "&")
n1 <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else suppressWarnings(as.numeric(unlist(x)[1]))

rows <- list()
for (f in sort(list.files(jdir, "^chunk_[0-9]+\\.json$", full.names = TRUE))) {
  for (v in fromJSON(f, simplifyVector = FALSE)) {
    alt <- sub("^.*/", "", g(v, "allele_string") %||% "")
    ## colocated variants: rsID, ClinVar, gnomAD 頻度
    ev <- c(); cs <- c(); fq <- list()
    for (cv in g(v, "colocated_variants")) {
      ev <- c(ev, g(cv, "id")); cs <- c(cs, unlist(g(cv, "clin_sig")))
      fr <- g(cv, "frequencies"); if (!is.null(fr) && !is.null(fr[[alt]])) fq <- fr[[alt]]
    }
    fqn <- unlist(fq); fqn <- fqn[grepl("^gnomad", names(fqn))]
    maxaf <- if (length(fqn)) max(as.numeric(fqn)) else NA_real_
    maxpop <- if (length(fqn)) paste(names(fqn)[as.numeric(fqn) == maxaf], collapse = "&") else NA_character_
    for (tc in g(v, "transcript_consequences")) {
      am <- g(tc, "alphamissense"); sp <- g(tc, "spliceai")
      rows[[length(rows) + 1]] <- data.table(
        id = g(v, "id"), SYMBOL = s1(g(tc, "gene_symbol")), Gene = s1(g(tc, "gene_id")),
        Feature = s1(g(tc, "transcript_id")),
        MANE_SELECT = s1(g(tc, "mane_select")), CANONICAL = s1(g(tc, "canonical")),
        BIOTYPE = s1(g(tc, "biotype")), Consequence = s1(g(tc, "consequence_terms")),
        IMPACT = s1(g(tc, "impact")), HGVSc = s1(g(tc, "hgvsc")), HGVSp = s1(g(tc, "hgvsp")),
        EXON = s1(g(tc, "exon")), INTRON = s1(g(tc, "intron")),
        Protein_position = s1(g(tc, "protein_start")), Amino_acids = s1(g(tc, "amino_acids")),
        SIFT = s1(paste0(g(tc, "sift_prediction") %||% "", "(", g(tc, "sift_score") %||% "", ")")),
        PolyPhen = s1(paste0(g(tc, "polyphen_prediction") %||% "", "(", g(tc, "polyphen_score") %||% "", ")")),
        gnomADe_AF = n1(fq[["gnomade"]]), gnomADe_EAS_AF = n1(fq[["gnomade_eas"]]),
        gnomADg_AF = n1(fq[["gnomadg"]]), gnomADg_EAS_AF = n1(fq[["gnomadg_eas"]]),
        MAX_AF = maxaf, MAX_AF_POPS = maxpop,
        CLIN_SIG = s1(unique(cs)), Existing_variation = s1(unique(ev)),
        am_class = s1(g(tc, "am_class") %||% g(am, "am_class")),
        am_pathogenicity = n1(g(tc, "am_pathogenicity") %||% g(am, "am_pathogenicity")),
        REVEL = n1(g(tc, "revel")), CADD_PHRED = n1(g(tc, "cadd_phred")),
        SpliceAI_pred_DS_AG = n1(g(sp, "DS_AG")), SpliceAI_pred_DS_AL = n1(g(sp, "DS_AL")),
        SpliceAI_pred_DS_DG = n1(g(sp, "DS_DG")), SpliceAI_pred_DS_DL = n1(g(sp, "DS_DL")),
        LoF = s1(g(tc, "lof")), LoF_filter = s1(g(tc, "lof_filter")), LoF_flags = s1(g(tc, "lof_flags")))
    }
  }
}
d <- rbindlist(rows, fill = TRUE)
d[SIFT == "()", SIFT := NA][PolyPhen == "()", PolyPhen := NA]
## 1 変異 x 1 遺伝子: MANE Select > canonical > 先頭
d[, pri := fifelse(!is.na(MANE_SELECT), 1L, fifelse(CANONICAL %in% c("1", "YES", "TRUE"), 2L, 3L))]
setorder(d, id, SYMBOL, pri)
d <- d[, .SD[1], by = .(id, SYMBOL)][, pri := NULL]
setnames(d, "id", "Uploaded_variation")
fwrite(d, out, sep = "\t", na = "-")
cat("[07c] wrote", out, ":", nrow(d), "rows,", uniqueN(d$Uploaded_variation), "variants\n")
cat("[07c] filled: AlphaMissense", d[!is.na(am_pathogenicity), .N], "| REVEL", d[!is.na(REVEL), .N],
    "| CADD", d[!is.na(CADD_PHRED), .N], "| SpliceAI", d[!is.na(SpliceAI_pred_DS_AG), .N],
    "| gnomADe", d[!is.na(gnomADe_AF), .N], "\n")
