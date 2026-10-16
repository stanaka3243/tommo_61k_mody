#!/usr/bin/env Rscript
# =====================================================================
# 08_merge_vep.R
#  03 の candidates.tsv に 07 の VEP 結果 (out/vep/{panel,off}.vep.tsv) を結合し,
#  ACMG 手作業分類用の表を作る。
#  使い方: Rscript 08_merge_vep.R <OUTDIR>
#  出力:
#    candidates_vep.tsv        全候補 + VEP 列 (HGVS, gnomAD, SIFT/PolyPhen, AlphaMissense など)
#    review_missense_vep.tsv   ultra-rare missense のみ (AlphaMissense / REVEL の高い順)
#    plp_core4_vep.tsv         MODY1/2/3/5 (HNF4A, GCK, HNF1A, HNF1B) の ClinVar P/LP + PTV
# =====================================================================
suppressPackageStartupMessages(library(data.table))
outdir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(outdir)) stop("usage: Rscript 08_merge_vep.R OUTDIR")

## ---- 03 の候補 ----
cand <- fread(file.path(outdir, "candidates.tsv"), colClasses = list(character = c("chrom","pos","ref","alt")))
cand[, id := paste(chrom, pos, ref, alt, sep = "_")]
cand[, AF := fifelse(AN > 0, AC / AN, NA_real_)]
cand[, MAF := pmin(AF, 1 - AF)]

## ---- VEP の TSV を読む ("##" 行を飛ばし, "-" を NA に) ----
read_vep <- function(f) {
  if (!file.exists(f)) { warning("not found: ", f); return(NULL) }
  d <- fread(cmd = paste("grep -v '^##'", shQuote(f)), sep = "\t", na.strings = c("-", ""),
             colClasses = "character")
  setnames(d, 1, "id")
  d
}
vep <- rbindlist(list(read_vep(file.path(outdir, "vep", "panel.vep.tsv")),
                      read_vep(file.path(outdir, "vep", "off.vep.tsv"))), fill = TRUE)
vep <- unique(vep)

## --per_gene でも 1 変異が複数遺伝子に重なることがあるので, 候補の遺伝子と SYMBOL を合わせる
v <- merge(cand, vep, by.x = c("id", "gene"), by.y = c("id", "SYMBOL"), all.x = TRUE)
v <- v[, .SD[1], by = .(id, gene, source)]

## 数値化
numc <- intersect(c("gnomADe_AF","gnomADe_EAS_AF","gnomADg_AF","gnomADg_EAS_AF","MAX_AF",
                    "am_pathogenicity","REVEL","CADD_PHRED",
                    "SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG","SpliceAI_pred_DS_DL"),
                  names(v))
for (cc in numc) set(v, j = cc, value = suppressWarnings(as.numeric(v[[cc]])))
if (all(c("SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG","SpliceAI_pred_DS_DL") %in% names(v)))
  v[, SpliceAI_max := pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL, na.rm = TRUE)]

## gnomAD に無い / 極めて稀 (ACMG PM2_Supporting の目安)
v[, gnomAD_absent := (is.na(gnomADe_AF) | gnomADe_AF == 0) & (is.na(gnomADg_AF) | gnomADg_AF == 0)]

## ---- 列の並び ----
front <- c("source","gene","id","HGVSc","HGVSp","Consequence","IMPACT","class",
           "AC","AN","AF","MAF","carriers","low_AN",
           "gnomADe_AF","gnomADe_EAS_AF","gnomADg_AF","gnomADg_EAS_AF","MAX_AF","MAX_AF_POPS","gnomAD_absent",
           "am_class","am_pathogenicity","REVEL","CADD_PHRED","SpliceAI_max","SIFT","PolyPhen",
           "LoF","LoF_filter","LoF_flags",
           "cv_sig","cv_rev","cv_dn","CLIN_SIG","Existing_variation",
           "MANE_SELECT","Feature","EXON","INTRON","Protein_position","Amino_acids")
front <- intersect(front, names(v))
setcolorder(v, c(front, setdiff(names(v), front)))
setorder(v, gene, source, id)
fwrite(v, file.path(outdir, "candidates_vep.tsv"), sep = "\t")

## ultra-rare missense: AlphaMissense (なければ REVEL) の高い順
mis <- v[class == "ultra_rare_missense"]
ord <- intersect(c("am_pathogenicity","REVEL"), names(mis))
if (length(ord)) setorderv(mis, ord, order = rep(-1L, length(ord)), na.last = TRUE)
fwrite(mis, file.path(outdir, "review_missense_vep.tsv"), sep = "\t")

## MODY1/2/3/5 の ClinVar P/LP + PTV
core <- v[gene %in% c("HNF4A","GCK","HNF1A","HNF1B") & class %in% c("ClinVar_PLP","PTV")]
fwrite(core, file.path(outdir, "plp_core4_vep.tsv"), sep = "\t")

## ---- 画面表示 ----
cat("VEP matched:", v[!is.na(Consequence), .N], "/", nrow(v), "candidates\n")
if ("am_class" %in% names(mis)) {
  cat("\nultra-rare missense x AlphaMissense class (61KJPN):\n")
  print(mis[source == "61KJPN", .N, by = .(gene, am_class)][order(gene, am_class)])
}
cat("\nOutputs:\n ", paste(file.path(outdir, c("candidates_vep.tsv","review_missense_vep.tsv","plp_core4_vep.tsv")), collapse = "\n  "), "\n")
