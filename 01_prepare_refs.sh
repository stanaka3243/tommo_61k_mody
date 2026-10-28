#!/usr/bin/env bash
# 参照データの準備: Ensembl GFF3, MANE, 抽出領域 (MANE Select 転写産物 + 上流/下流),
#                   参照配列 (対象染色体のみ), ClinVar
# インターネット接続のあるノードで実行する
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$REFDIR"; cd "$REFDIR"

## 1) Ensembl GFF3 (染色体名に chr を付ける)
GFF="Homo_sapiens.GRCh38.${ENSEMBL_REL}.gff3.gz"
[ -s "$GFF" ] || wget -q "https://ftp.ensembl.org/pub/release-${ENSEMBL_REL}/gff3/homo_sapiens/$GFF"
if [ ! -s ensembl.chr.gff3.gz ]; then
  zcat "$GFF" | sed -E '/^#/!s/^/chr/' | bgzip > ensembl.chr.gff3.gz
fi

## 2) MANE Select リスト (NCBI)
MANE="MANE.GRCh38.v${MANE_VER}.summary.txt.gz"
[ -s "$MANE" ] || wget -q "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANE_VER}/$MANE"
ln -sf "$MANE" mane_summary.txt.gz

## 3) 抽出領域: MANE Select 転写産物の全長 + 上流 UPSTREAM bp + 下流 DOWNSTREAM bp (鎖の向きを考慮)
##    (a) 遺伝子 -> MANE Select の Ensembl 転写産物 ID
GENE_LIST=$(awk -F'\t' 'NR>1{print $1}' "$GENES" | paste -sd, -)
zcat mane_summary.txt.gz | awk -F'\t' -v OFS='\t' -v gl="$GENE_LIST" '
  BEGIN { n = split(gl, a, ","); for (i = 1; i <= n; i++) want[a[i]] = 1 }
  NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); col[h] = i }; next }
  $col["MANE_status"] == "MANE Select" && ($col["symbol"] in want) {
    tx = $col["Ensembl_nuc"]; sub(/\..*$/, "", tx); print $col["symbol"], tx, $col["RefSeq_nuc"] }
' > mane_select.tsv
for g in ${GENE_LIST//,/ }; do
  grep -q "^$g"$'\t' mane_select.tsv || { echo "[01] ERROR: MANE Select not found for $g" >&2; exit 1; }
done

##    (b) GFF3 から転写産物の座標と, その CDS / exon を取り出す (1-based)
##        mane_tx.tsv : gene tx refseq chrom tx_start tx_end strand ext_start ext_end
##        features.tsv: gene type(CDS/exon) start end
zcat ensembl.chr.gff3.gz | awk -F'\t' -v OFS='\t' -v up="$UPSTREAM" -v dn="$DOWNSTREAM" '
  FNR == NR { gene[$2] = $1; refseq[$2] = $3; next }
  /^#/ || NF < 9 { next }
  $3 == "mRNA" && match($9, /ID=transcript:[^;]+/) {
    t = substr($9, RSTART + 14, RLENGTH - 14)
    if (t in gene) { chr[t] = $1; s[t] = $4; e[t] = $5; str[t] = $7 }
    next
  }
  ($3 == "CDS" || $3 == "exon") && match($9, /Parent=transcript:[^;]+/) {
    t = substr($9, RSTART + 18, RLENGTH - 18)
    if (t in gene) print gene[t], $3, $4, $5 > "features.tsv"
  }
  END {
    print "gene", "tx", "refseq", "chrom", "tx_start", "tx_end", "strand", "ext_start", "ext_end" > "mane_tx.tsv"
    for (t in gene) {
      if (!(t in chr)) { print "[01] ERROR: " t " (" gene[t] ") not in GFF3" > "/dev/stderr"; bad = 1; continue }
      if (str[t] == "+") { a = s[t] - up; b = e[t] + dn } else { a = s[t] - dn; b = e[t] + up }
      if (a < 1) a = 1
      print gene[t], t, refseq[t], chr[t], s[t], e[t], str[t], a, b > "mane_tx.tsv"
    }
    exit bad
  }
' mane_select.tsv -
sort -k1,1 -k2,2 -k3,3n features.tsv -o features.tsv
awk -F'\t' -v OFS='\t' 'NR > 1 { print $4, $8 - 1, $9, $1 }' mane_tx.tsv | sort -k1,1 -k2,2n > targets.chr.bed
sed 's/^chr//; s/^M\t/MT\t/' targets.chr.bed > targets.nochr.bed       # ClinVar 用 (chr なし)
CHRS=$(cut -f1 targets.chr.bed | sort -u)
echo "[01] targets (MANE span, up ${UPSTREAM} / down ${DOWNSTREAM} bp):"
column -t mane_tx.tsv 2>/dev/null || cat mane_tx.tsv

## 4) 参照配列 (UCSC hg38, 必要な染色体だけ)
for c in $CHRS; do
  [ -s "$c.fa.gz" ] || wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/$c.fa.gz"
done
if [ ! -s ref_sub.fa.fai ] || [ "$(cut -f1 ref_sub.fa.fai | sort | paste -sd' ')" != "$(echo $CHRS | tr ' ' '\n' | sort | paste -sd' ')" ]; then
  for c in $CHRS; do zcat "$c.fa.gz"; done > ref_sub.fa
  samtools faidx ref_sub.fa
fi

## 5) ClinVar (GRCh38) -> 対象領域だけ切り出し, chr 付きに変換
[ -s clinvar.vcf.gz ]     || wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
[ -s clinvar.vcf.gz.tbi ] || wget -q https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
{ for c in $(seq 1 22) X Y; do printf "%s\tchr%s\n" "$c" "$c"; done; printf "MT\tchrM\n"; } > chrmap.txt
bcftools view -R targets.nochr.bed clinvar.vcf.gz -Ou \
 | bcftools annotate --rename-chrs chrmap.txt -Ou \
 | bcftools sort -Oz -o clinvar.mody.chr.vcf.gz
bcftools index -t -f clinvar.mody.chr.vcf.gz
bcftools view -h clinvar.vcf.gz | grep -m1 '^##fileDate' > clinvar_version.txt || true
echo "[01] done: $REFDIR"
