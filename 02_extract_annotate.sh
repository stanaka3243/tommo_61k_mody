#!/usr/bin/env bash
# jMorp 61KJPN (本体) と off-panel singleton から MODY_BED の領域だけを取得し,
# 正規化 -> ClinVar 付与 -> bcftools csq で機能予測 -> TSV に書き出す
# (データ全体はダウンロードせず, HTTPS + tabix index で領域だけ読む)
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$OUTDIR/tmp"; cd "$OUTDIR/tmp"
rm -f panel.*.vcf.gz* off.vcf.gz* off.*.vcf.gz*

REF="$REFDIR/ref_sub.fa"
GFF="$REFDIR/ensembl.chr.gff3.gz"
CV="$REFDIR/clinvar.mody.chr.vcf.gz"
BED="$REFDIR/targets.chr.bed"          # 01 で整形した BED (chr 付き, PAD 済み)

## 1) 61KJPN 本体: 染色体ごとのファイルなので, 染色体ごとに BED を分けて取得 (FILTER=PASS)
for chr in $(cut -f1 "$BED" | sort -u); do
  awk -v c="$chr" '$1==c' "$BED" > "bed.$chr.bed"
  url="${AF_URL/\{CHR\}/$chr}"
  echo "[02] 61KJPN $chr ($(wc -l < "bed.$chr.bed") intervals)"
  bcftools view -f PASS -R "bed.$chr.bed" "$url" -Oz -o "panel.$chr.vcf.gz" </dev/null
  bcftools index -t -f "panel.$chr.vcf.gz"
done

## 2) off-panel singleton: 常染色体が 1 ファイルなので BED 全体で取得
echo "[02] off-panel singleton"
bcftools view -f PASS -R "$BED" "$OFF_URL" -Oz -o off.vcf.gz </dev/null
bcftools index -t -f off.vcf.gz

## 3) 結合 -> 多アレル分割/正規化 -> ClinVar -> csq
annotate () {   # $1 = 出力名, 残り = 入力 VCF
  local name=$1; shift
  # (a) 結合 -> ソート -> 多アレル分割/正規化  (annotate -a には索引付きファイルが必要なので一旦書き出す)
  bcftools concat -a -D "$@" -Ou \
  | bcftools sort -Ou \
  | bcftools norm -m -any -f "$REF" -Oz -o "$name.norm.vcf.gz"
  bcftools index -t -f "$name.norm.vcf.gz"
  # (b) ClinVar を CV_* として付与
  bcftools annotate -a "$CV" \
      -c 'INFO/CV_SIG:=INFO/CLNSIG,INFO/CV_REV:=INFO/CLNREVSTAT,INFO/CV_DN:=INFO/CLNDN,INFO/CV_ALLELEID:=INFO/ALLELEID' \
      "$name.norm.vcf.gz" -Oz -o "$name.cv.vcf.gz"
  bcftools index -t -f "$name.cv.vcf.gz"
  # (c) bcftools csq で機能予測 (INFO/BCSQ)
  bcftools csq -f "$REF" -g "$GFF" -l "$name.cv.vcf.gz" -Oz -o "$OUTDIR/$name.annot.vcf.gz"
  bcftools index -t -f "$OUTDIR/$name.annot.vcf.gz"
}
annotate panel panel.chr*.vcf.gz
annotate off   off.vcf.gz

## 4) TSV に書き出し (R で分類する)
#   panel: 61KJPN 本体 (血縁除外 61,424 人)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/nhomalt\t%INFO/CV_SIG\t%INFO/CV_REV\t%INFO/CV_DN\t%INFO/BCSQ\n' \
  "$OUTDIR/panel.annot.vcf.gz" > "$OUTDIR/panel.tsv"
#   off: 100KJPN で AC=1 だが 61KJPN では AC=0 の変異 (nhomalt は 0 とみなす)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t0\t%INFO/CV_SIG\t%INFO/CV_REV\t%INFO/CV_DN\t%INFO/BCSQ\n' \
  "$OUTDIR/off.annot.vcf.gz" > "$OUTDIR/off.tsv"
echo "[02] done: $OUTDIR/panel.tsv ($(wc -l < "$OUTDIR/panel.tsv") variants), $OUTDIR/off.tsv ($(wc -l < "$OUTDIR/off.tsv") variants)"
