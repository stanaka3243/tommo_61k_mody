#!/usr/bin/env Rscript
# 04 の出力から bin ごとの CN0/CN1/CN3+ 人数を計算し, HNF1B 遺伝子上の欠失 (CN1) 頻度を推定
# 使い方: Rscript 05_cnv_17q12.R <OUTDIR>
suppressPackageStartupMessages(library(data.table))
outdir <- commandArgs(trailingOnly = TRUE)[1]
d <- fread(file.path(outdir, "cnv_17q12.tsv"), header = FALSE,
           col.names = c("chrom","start","end","alt","SN","SC"), colClasses = "character")
cnt <- function(alt, sc, cn) {
  al <- gsub("[<>CN]", "", strsplit(alt, ",")[[1]]); n <- as.numeric(strsplit(sc, ",")[[1]])
  sum(n[as.numeric(al) %in% cn])
}
d[, `:=`(start = as.numeric(start), end = as.numeric(end), SN = as.numeric(SN))]
d[, `:=`(CN0 = mapply(cnt, alt, SC, MoreArgs = list(cn = 0)),
         CN1 = mapply(cnt, alt, SC, MoreArgs = list(cn = 1)),
         CN3p = mapply(cnt, alt, SC, MoreArgs = list(cn = 3:20)))]
fwrite(d[, .(chrom, start, end, SN, CN0, CN1, CN3p)], file.path(outdir, "cnv_17q12_bins.tsv"), sep = "\t")

h <- d[start >= 37686431 & end <= 37745105]          # HNF1B 遺伝子上の bin
cat("HNF1B bins:", nrow(h), "\n")
cat("CN1 (heterozygous loss) per bin: median =", median(h$CN1), " range =", range(h$CN1), "\n")
cat("-> estimated 1 in", round(median(h$SN) / max(median(h$CN1), 1)), "(N =", median(h$SN), ")\n")

png(file.path(outdir, "cnv_17q12.png"), width = 1400, height = 600, res = 120)
plot(d$start / 1e6, d$CN1, type = "h", col = "steelblue", xlab = "chr17 (Mb, GRCh38)",
     ylab = "Number of people with CN1", main = "JCNVv1: copy-number loss around 17q12")
rect(37.686431, 0, 37.745105, max(d$CN1, 1), border = "red")
text(37.716, max(d$CN1, 1), "HNF1B", col = "red", pos = 3, xpd = TRUE)
invisible(dev.off())
cat("plot:", file.path(outdir, "cnv_17q12.png"), "\n")
