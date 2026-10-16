#!/usr/bin/env bash
# 参照データの準備: 対象領域 BED の整形, 参照配列 (対象染色体のみ), Ensembl GFF3, MANE, ClinVar
# インターネット接続のあるノードで実行する
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$REFDIR"; cd "$REFDIR"

## 0) 対象領域 BED を整形: chr 付き / 前後 PAD bp / ソート / マージ
##    (元の BED の染色体名が "7" でも "chr7" でも動く)
grep -v -E '^(#|track|browser)' "$MODY_BED" \
 | awk -v p="$PAD" 'BEGIN{OFS="\t"} {c=$1; if (c !~ /^chr/) c="chr"c; s=$2-p; if (s<0) s=0; print c, s, $3+p}' \
 | sort -k1,1 -k2,2n \
 | awk 'BEGIN{OFS="\t"} $1==c && $2<=e {if ($3>e) e=$3; next} {if (c!="") print c,s,e; c=$1; s=$2; e=$3} END{if (c!="") print c,s,e}' \
 > targets.chr.bed
sed 's/^chr//; s/^M\t/MT\t/' targets.chr.bed > targets.nochr.bed       # ClinVar 用 (chr なし)
CHRS=$(cut -f1 targets.chr.bed | sort -u)
echo "[01] targets: $(wc -l < targets.chr.bed) intervals, $(awk '{s+=$3-$2} END{print s}' targets.chr.bed) bp on: $CHRS"

## 1) 参照配列 (UCSC hg38, 必要な染色体だけ)
for c in $CHRS; do
  [ -s "$c.fa.gz" ] || wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$c.fa.gz"
done
if [ ! -s ref_sub.fa.fai ]; then
  for c in $CHRS; do zcat "$c.fa.gz"; done > ref_sub.fa
  samtools faidx ref_sub.fa
fi

## 2) Ensembl GFF3 (染色体名に chr を付ける)
GFF="Homo_sapiens.GRCh38.${ENSEMBL_REL}.gff3.gz"
[ -s "$GFF" ] || wget -q "https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gff3/homo_sapiens/$GFF"
if [ ! -s ensembl.chr.gff3.gz ]; then
  zcat "$GFF" | sed -E '/^#/!s/^/chr/' | bgzip > ensembl.chr.gff3.gz
fi

## 3) MANE Select リスト (NCBI)
MANE="MANE.GRCh38.v${MANE_VER}.summary.txt.gz"
[ -s "$MANE" ] || wget -q "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANE_VER}/$MANE"
ln -sf "$MANE" mane_summary.txt.gz

## 4) ClinVar (GRCh38) -> 対象領域だけ切り出し, chr 付きに変換
[ -s clinvar.vcf.gz ]     || wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
[ -s clinvar.vcf.gz.tbi ] || wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
{ for c in $(seq 1 22) X Y; do printf "%s\tchr%s\n" "$c" "$c"; done; printf "MT\tchrM\n"; } > chrmap.txt
bcftools view -R targets.nochr.bed clinvar.vcf.gz -Ou \
 | bcftools annotate --rename-chrs chrmap.txt -Ou \
 | bcftools sort -Oz -o clinvar.mody.chr.vcf.gz
bcftools index -t -f clinvar.mody.chr.vcf.gz
bcftools view -h clinvar.vcf.gz | grep -m1 '^##fileDate' > clinvar_version.txt || true
echo "[01] done: $REFDIR"
