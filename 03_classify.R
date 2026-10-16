#!/usr/bin/env Rscript
# =====================================================================
# 03_classify.R
#  02 で作った panel.tsv / off.tsv を読み, MANE Select 転写産物の結果だけを使って
#  バリアントを分類し, 遺伝子ごとの保因者数と頻度 (1 in X, 95%CI) を出す。
#  使い方: Rscript 03_classify.R <OUTDIR> <REFDIR> <GENES.tsv> [MAX_AC_MISSENSE] [N_PANEL]
# =====================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3) stop("usage: Rscript 03_classify.R OUTDIR REFDIR GENES.tsv [MAX_AC] [N_PANEL]")
outdir  <- a[1]; refdir <- a[2]; genes_f <- a[3]
max_ac  <- if (length(a) >= 4) as.integer(a[4]) else 2L
n_panel <- if (length(a) >= 5) as.integer(a[5]) else 61424L

genes <- fread(genes_f)$gene

## ---- MANE Select (Ensembl transcript, バージョン番号を外す) ----
mane <- fread(cmd = paste("zcat", shQuote(file.path(refdir, "mane_summary.txt.gz"))),
              check.names = TRUE)
mane <- mane[MANE_status == "MANE Select" & symbol %in% genes,
             .(gene = symbol, tx = sub("\\..*$", "", Ensembl_nuc), refseq = RefSeq_nuc)]
missing <- setdiff(genes, mane$gene)
if (length(missing)) warning("MANE Select が見つからない遺伝子: ", paste(missing, collapse = ", "))

## ---- 02 の TSV を読む ----
cols <- c("chrom","pos","ref","alt","AC","AN","nhomalt","cv_sig","cv_rev","cv_dn","bcsq")
rd <- function(f, src) {
  d <- fread(file.path(outdir, f), header = FALSE, col.names = cols, sep = "\t",
             colClasses = "character", na.strings = ".")
  d[, source := src]
}
v <- rbind(rd("panel.tsv", "61KJPN"), rd("off.tsv", "off_panel_singleton"))
num <- function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x }
v[, `:=`(AC = num(AC), AN = num(AN), nhomalt = num(nhomalt))]
v <- v[AC > 0]
v[, carriers := AC - nhomalt]           # ヘテロ + ホモ の人数

## ---- BCSQ を展開: consequence|gene|transcript|biotype|strand|aa|dna ----
v[, vid := .I]
b <- v[!is.na(bcsq), .(csq = unlist(strsplit(bcsq, ","))), by = vid]
b <- b[!startsWith(csq, "@")]                      # 他レコード参照 (@pos) は除外
b[, f := strsplit(csq, "|", fixed = TRUE)]
b[, `:=`(consequence = sapply(f, `[`, 1), gene = sapply(f, `[`, 2),
         tx = sapply(f, `[`, 3), aa = sapply(f, `[`, 6), dna = sapply(f, `[`, 7))]
b[, f := NULL]
b[, tx := sub("\\..*$", "", sub("^transcript:", "", tx))]   # バージョン/接頭辞を外す
b <- merge(b, mane, by = c("gene", "tx"))          # MANE Select の行だけ残す
b <- b[, .SD[1], by = vid]                         # 1 変異 1 行
x <- merge(v, b[, .(vid, gene, tx, refseq, consequence, aa, dna)], by = "vid")

## ---- 分類 ----
ptv_re <- "stop_gained|frameshift|splice_acceptor|splice_donor|start_lost"
x[, clinvar_plp := !is.na(cv_sig) & grepl("Pathogenic|Likely_pathogenic", cv_sig) &
                   !grepl("Conflicting", cv_sig)]
x[, is_ptv := grepl(ptv_re, consequence)]
x[, is_mis := grepl("missense", consequence)]
x[, class := fifelse(clinvar_plp, "ClinVar_PLP",
             fifelse(is_ptv, "PTV",
             fifelse(is_mis & AC <= max_ac, "ultra_rare_missense", NA_character_)))]
cand <- x[!is.na(class)][order(gene, as.numeric(pos))]
cand[, low_AN := AN < 0.9 * max(AN), by = source]  # 呼び出し率が低い部位の目印

keep <- c("source","gene","chrom","pos","ref","alt","tx","refseq","consequence","dna","aa",
          "AC","AN","nhomalt","carriers","cv_sig","cv_rev","cv_dn","class","low_AN")
fwrite(cand[, ..keep], file.path(outdir, "candidates.tsv"), sep = "\t")
fwrite(cand[class == "ultra_rare_missense", ..keep],
       file.path(outdir, "review_missense.tsv"), sep = "\t")   # ACMG で手作業分類する

## ---- 集計: 61KJPN 本体 (分母 = N_PANEL) ----
pan <- cand[source == "61KJPN"]
summ <- function(d, label) {
  s <- d[, .(variants = .N, carriers = sum(carriers)), by = gene]
  s <- merge(data.table(gene = genes), s, by = "gene", all.x = TRUE)
  s[is.na(variants), `:=`(variants = 0L, carriers = 0)]
  tot <- data.table(gene = "TOTAL", variants = sum(s$variants), carriers = sum(s$carriers))
  s <- rbind(s, tot)
  ci <- t(sapply(s$carriers, function(k) binom.test(k, n_panel)$conf.int))
  s[, `:=`(set = label, N = n_panel, prev = carriers / n_panel,
           lo95 = ci[, 1], hi95 = ci[, 2],
           one_in = fifelse(carriers > 0, round(n_panel / carriers), NA_real_))]
  s[]
}
res <- rbind(
  summ(pan[class == "ClinVar_PLP"],                        "A: ClinVar P/LP only"),
  summ(pan[class %in% c("ClinVar_PLP", "PTV")],            "B: ClinVar P/LP + PTV"),
  summ(pan,                                                "C: B + all ultra-rare missense (upper bound)"))
fwrite(res, file.path(outdir, "summary_61KJPN.tsv"), sep = "\t")

## off-panel singleton は別集団 (100KJPN, AN ~ 20万) なので別に数だけ出す
off <- cand[source == "off_panel_singleton", .(variants = .N, carriers = sum(carriers)), by = .(gene, class)]
fwrite(off, file.path(outdir, "summary_off_panel_singleton.tsv"), sep = "\t")

## ---- 画面表示 ----
cat("\n== 61KJPN: carriers by gene (N =", n_panel, ") ==\n")
print(dcast(res, gene ~ set, value.var = "carriers"))
cat("\nSharp et al. (UKB): 1 in 1,052 overall; GCK 1 in 2,787  ->  expected in N:",
    round(n_panel / 1052), "overall,", round(n_panel / 2787), "GCK\n")
cat("\nOutputs:", file.path(outdir, c("candidates.tsv", "review_missense.tsv",
    "summary_61KJPN.tsv", "summary_off_panel_singleton.tsv")), sep = "\n  ")
