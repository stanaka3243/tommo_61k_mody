#!/usr/bin/env bash
# jMorp 61KJPN (本体) と off-panel singleton から抽出領域 (01 の targets.chr.bed) を取得し,
# 正規化 -> ClinVar 付与 -> bcftools csq で機能予測 -> TSV に書き出す
# 抽出時の絞り込みは行わない (FILTER も列として残し, 絞り込みは 03 以降で行う)
#   DATA_SITE=slade    : 61KJPN はローカルの VCF (AF_LOCAL) を読む
#   DATA_SITE=shirokane: 61KJPN は jMorp から HTTPS + tabix index で領域だけ読む
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$OUTDIR/tmp"; cd "$OUTDIR/tmp"
rm -f panel.*.vcf.gz* off.vcf.gz* off.*.vcf.gz*

REF="$REFDIR/ref_sub.fa"
GFF="$REFDIR/ensembl.chr.gff3.gz"
CV="$REFDIR/clinvar.mody.chr.vcf.gz"
BED="$REFDIR/targets.chr.bed"          # 01 で作った抽出領域 (chr 付き)

case "$DATA_SITE" in
  slade)     AF_SRC="$AF_LOCAL" ;;
  shirokane) AF_SRC="$AF_URL" ;;
  *) echo "[02] ERROR: DATA_SITE must be slade or shirokane" >&2; exit 1 ;;
esac
OFF_SRC="${OFF_LOCAL:-$OFF_URL}"

## 1) 61KJPN 本体: 染色体ごとのファイルなので, 染色体ごとに BED を分けて取得 (FILTER では絞らない)
for chr in $(cut -f1 "$BED" | sort -u); do
  awk -v c="$chr" '$1==c' "$BED" > "bed.$chr.bed"
  src="${AF_SRC/\{CHR\}/$chr}"
  if [ "$DATA_SITE" = slade ]; then
    [ -s "$src" ] || { echo "[02] ERROR: $src not found" >&2; exit 1; }
    [ -e "$src.tbi" ] || [ -e "$src.csi" ] || bcftools index -t "$src"
  fi
  echo "[02] 61KJPN $chr ($(wc -l < "bed.$chr.bed") regions) <- $src"
  bcftools view -R "bed.$chr.bed" "$src" -Oz -o "panel.$chr.vcf.gz" </dev/null
  bcftools index -t -f "panel.$chr.vcf.gz"
done

## 2) off-panel singleton: 常染色体が 1 ファイルなので BED 全体で取得
echo "[02] off-panel singleton <- $OFF_SRC"
bcftools view -R "$BED" "$OFF_SRC" -Oz -o off.vcf.gz </dev/null
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
##    列: chrom pos ref alt AC AN nhomalt cv_sig cv_rev cv_dn bcsq id filter AF_vcf platform_bias
##    (ヘッダーに無い INFO 項目は "." を出す; off-panel の nhomalt は 0 とみなす)
qfmt () {   # $1 = VCF, $2 = nhomalt が無いときの値
  local h; h=$(bcftools view -h "$1")
  has () { grep -q "^##INFO=<ID=$1," <<< "$h"; }
  local nh="$2" af='.' pb='.'
  has nhomalt && nh='%INFO/nhomalt'
  has AF && af='%INFO/AF'
  has TOMMO_POSSIBLE_PLATFORM_BIAS && pb='%INFO/TOMMO_POSSIBLE_PLATFORM_BIAS'
  printf '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%INFO/AC\\t%%INFO/AN\\t%s\\t%%INFO/CV_SIG\\t%%INFO/CV_REV\\t%%INFO/CV_DN\\t%%INFO/BCSQ\\t%%ID\\t%%FILTER\\t%s\\t%s\\n' "$nh" "$af" "$pb"
}
#   panel: 61KJPN 本体 (血縁除外 61,424 人)
bcftools query -f "$(qfmt "$OUTDIR/panel.annot.vcf.gz" .)" "$OUTDIR/panel.annot.vcf.gz" > "$OUTDIR/panel.tsv"
#   off: 100KJPN で AC=1 だが 61KJPN では AC=0 の変異
bcftools query -f "$(qfmt "$OUTDIR/off.annot.vcf.gz" 0)" "$OUTDIR/off.annot.vcf.gz" > "$OUTDIR/off.tsv"
echo "[02] done: $OUTDIR/panel.tsv ($(wc -l < "$OUTDIR/panel.tsv") variants), $OUTDIR/off.tsv ($(wc -l < "$OUTDIR/off.tsv") variants)"
