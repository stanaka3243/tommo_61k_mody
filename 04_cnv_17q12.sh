#!/usr/bin/env bash
# JCNVv1 (48,874 人, 1 kb bin) から 17q12 周辺の bin ごとのコピー数分布を取り出す
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$REFDIR" "$OUTDIR"; cd "$REFDIR"
[ -s jcnvv1.vcf.gz ]     || wget -q -O jcnvv1.vcf.gz     "$CNV_URL"
[ -s jcnvv1.vcf.gz.tbi ] || wget -q -O jcnvv1.vcf.gz.tbi "$CNV_URL.tbi"
# 17q12 周辺 (表示用の広め領域; HNF1B は chr17:37,686,431-37,745,105)
bcftools query -r chr17:35500000-38500000 \
  -f '%CHROM\t%POS\t%INFO/END\t%ALT\t%INFO/SN\t%INFO/SC\n' jcnvv1.vcf.gz > "$OUTDIR/cnv_17q12.tsv"
echo "[04] done: $OUTDIR/cnv_17q12.tsv"
