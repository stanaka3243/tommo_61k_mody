#!/usr/bin/env Rscript
# =====================================================================
# 03_classify.R
#  02 で作った panel.tsv / off.tsv を読み, 全バリアントに判定列を付けて書き出す
#  (all_variants.tsv: 絞り込みなし)。そのうえで, 候補の条件
#  (既定: FILTER=PASS かつ REGION=CDS_pm{PAD}) を満たすものを従来どおり分類し,
#  遺伝子ごとの保因者数と頻度 (1 in X, 95%CI) を出す。
#  使い方: Rscript 03_classify.R <OUTDIR> <REFDIR> <GENES.tsv> [MAX_AC_MISSENSE] [N_PANEL]
#                                [PAD] [CAND_PASS_ONLY] [CAND_REGIONS]
#    CAND_PASS_ONLY : 1 = 候補を FILTER=PASS に限る, 0 = 限らない
#    CAND_REGIONS   : 候補に使う REGION (カンマ区切り; 既定 CDS_pm{PAD})
# =====================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 3) stop("usage: Rscript 03_classify.R OUTDIR REFDIR GENES.tsv [MAX_AC] [N_PANEL] [PAD] [CAND_PASS_ONLY] [CAND_REGIONS]")
outdir  <- a[1]; refdir <- a[2]; genes_f <- a[3]
max_ac  <- if (length(a) >= 4) as.integer(a[4]) else 2L
n_panel <- if (length(a) >= 5) as.integer(a[5]) else 61424L
pad     <- if (length(a) >= 6) as.integer(a[6]) else 10L
pass_only    <- if (length(a) >= 7) a[7] == "1" else TRUE
cand_regions <- if (length(a) >= 8 && nzchar(a[8])) strsplit(a[8], ",")[[1]] else paste0("CDS_pm", pad)
cds_label    <- paste0("CDS_pm", pad)

genes <- fread(genes_f)$gene

## ---- MANE Select 転写産物と抽出領域 (01 で作成) ----
tx <- fread(file.path(refdir, "mane_tx.tsv"))
tx <- tx[gene %in% genes]
missing <- setdiff(genes, tx$gene)
if (length(missing)) warning("MANE Select が見つからない遺伝子: ", paste(missing, collapse = ", "))
ft <- fread(file.path(refdir, "features.tsv"), header = FALSE,
            col.names = c("gene", "type", "start", "end"))[gene %in% genes]

## ---- 02 の TSV を読む ----
cols <- c("chrom","pos","ref","alt","AC","AN","nhomalt","cv_sig","cv_rev","cv_dn","bcsq",
          "rsid","filter","AF_vcf","platform_bias")
rd <- function(f, src) {
  d <- fread(file.path(outdir, f), header = FALSE, col.names = cols, sep = "\t",
             colClasses = "character", na.strings = ".")
  d[, source := src]
}
v <- rbind(rd("panel.tsv", "61KJPN"), rd("off.tsv", "off_panel_singleton"))
num <- function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x }
v[, `:=`(AC = num(AC), AN = num(AN), nhomalt = num(nhomalt))]
v[, carriers := AC - nhomalt]           # ヘテロ + ホモ の人数
v[, AF := fifelse(AN > 0, signif(AC / AN, 6), NA_real_)]
v[, `:=`(pos_i = as.integer(pos), vend = as.integer(pos) + nchar(ref) - 1L, vid = .I)]

## ---- 遺伝子の割り当て (抽出領域との重なり) と REGION ----
x <- tx[v, on = .(chrom, ext_start <= vend, ext_end >= pos_i), nomatch = 0L,
        .(vid = i.vid, gene = x.gene, tx = x.tx, refseq = x.refseq,
          tx_start = x.tx_start, tx_end = x.tx_end, strand = x.strand)]
x <- merge(x, v, by = "vid")
hit <- function(ftype, pad_bp) {
  f <- ft[type == ftype, .(gene, s = start - pad_bp, e = end + pad_bp)]
  h <- f[x, on = .(gene, s <= vend, e >= pos_i), nomatch = 0L, .(vid = i.vid, gene = i.gene)]
  paste(x$vid, x$gene) %in% paste(h$vid, h$gene)
}
x[, in_cds  := hit("CDS", pad)]
x[, in_exon := hit("exon", 0L)]
x[, region := fcase(in_cds, cds_label,
                    in_exon, "UTR_exon",
                    vend >= tx_start & pos_i <= tx_end, "intron",
                    (vend < tx_start) == (strand == "+"), "upstream",
                    default = "downstream")]

## ---- BCSQ: MANE Select 転写産物の consequence (consequence|gene|transcript|biotype|strand|aa|dna) ----
b <- v[!is.na(bcsq), .(csq = unlist(strsplit(bcsq, ","))), by = vid]
b <- b[!startsWith(csq, "@")]                      # 他レコード参照 (@pos) は除外
b[, f := strsplit(csq, "|", fixed = TRUE)]
b[, `:=`(consequence = sapply(f, `[`, 1), tx = sapply(f, `[`, 3),
         aa = sapply(f, `[`, 6), dna = sapply(f, `[`, 7))]
b[, tx := sub("\\..*$", "", sub("^transcript:", "", tx))]   # バージョン/接頭辞を外す
b <- b[!is.na(tx) & tx != "", .SD[1], by = .(vid, tx)]
x <- merge(x, b[, .(vid, tx, consequence, aa, dna)], by = c("vid", "tx"), all.x = TRUE)
x[, consequence := sub("^\\*", "", consequence)]

## ---- 判定列 (絞り込みはしない) ----
ptv_re <- "stop_gained|frameshift|splice_acceptor|splice_donor|start_lost"
x[, pass := !is.na(filter) & filter == "PASS"]
x[, clinvar_plp := !is.na(cv_sig) & grepl("Pathogenic|Likely_pathogenic", cv_sig) &
                   !grepl("Conflicting", cv_sig)]
x[, is_ptv := !is.na(consequence) & grepl(ptv_re, consequence)]
x[, is_mis := !is.na(consequence) & grepl("missense", consequence)]
x[, class := fifelse(clinvar_plp, "ClinVar_PLP",
             fifelse(is_ptv, "PTV",
             fifelse(is_mis & AC <= max_ac, "ultra_rare_missense", NA_character_)))]
x[, low_AN := AN < 0.9 * max(AN), by = source]   # 呼び出し率が低い部位の目印
setorder(x, source, gene, pos_i, ref, alt)

all_cols <- c("source","gene","chrom","pos","rsid","ref","alt","filter","pass","platform_bias",
              "region","tx","refseq","consequence","dna","aa",
              "AC","AN","AF","AF_vcf","nhomalt","carriers",
              "cv_sig","cv_rev","cv_dn","clinvar_plp","class","low_AN")
fwrite(x[, ..all_cols], file.path(outdir, "all_variants.tsv"), sep = "\t", na = ".", quote = FALSE)

ov <- x[, .N, by = .(source, gene, region, FILTER = fifelse(pass, "PASS", "nonPASS"))]
setorder(ov, source, gene, region, FILTER)
fwrite(ov, file.path(outdir, "overview_counts.tsv"), sep = "\t")

## ---- 候補 (従来の解析): 条件を満たすものを分類 ----
cand <- x[!is.na(class) & AC > 0 & region %in% cand_regions & (pass | !pass_only)]
keep <- c("source","gene","chrom","pos","ref","alt","tx","refseq","consequence","dna","aa",
          "AC","AN","nhomalt","carriers","cv_sig","cv_rev","cv_dn","class","low_AN",
          "rsid","filter","region","platform_bias")      # 07b は先頭 6 列の位置を使い, 08 は id 列を自前で作るので, 追加列は末尾に rsid の名前で置く
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
res[, candidate_filter := paste0(if (pass_only) "PASS; " else "any FILTER; ",
                                 "REGION in ", paste(cand_regions, collapse = ","))]
fwrite(res, file.path(outdir, "summary_61KJPN.tsv"), sep = "\t")

## off-panel singleton は別集団 (100KJPN, AN ~ 20万) なので別に数だけ出す
off <- cand[source == "off_panel_singleton", .(variants = .N, carriers = sum(carriers)), by = .(gene, class)]
fwrite(off, file.path(outdir, "summary_off_panel_singleton.tsv"), sep = "\t")

## ---- 画面表示 ----
cat("\n== all variants (no filtering): source x gene x REGION x FILTER ==\n")
print(dcast(ov, source + gene + region ~ FILTER, value.var = "N", fill = 0))
cat("\n== candidates: ", if (pass_only) "PASS" else "any FILTER", ", REGION in ",
    paste(cand_regions, collapse = ","), " ==\n", sep = "")
cat("== 61KJPN: carriers by gene (N =", n_panel, ") ==\n")
print(dcast(res, gene ~ set, value.var = "carriers"))
cat("\nSharp et al. (UKB): 1 in 1,052 overall; GCK 1 in 2,787  ->  expected in N:",
    round(n_panel / 1052), "overall,", round(n_panel / 2787), "GCK\n")
cat("\nOutputs:", file.path(outdir, c("all_variants.tsv", "overview_counts.tsv", "candidates.tsv",
    "review_missense.tsv", "summary_61KJPN.tsv", "summary_off_panel_singleton.tsv")), sep = "\n  ")
cat("\n")
