#!/usr/bin/env bash
# VEP の準備 (初回のみ): Singularity イメージ, オフラインキャッシュ, FASTA, プラグイン, AlphaMissense
# インターネット接続のあるノードで実行。キャッシュ展開に ~30 GB 必要。
set -euo pipefail
source "$(dirname "$0")/config.sh"
set +u; module load apptainer 2>/dev/null || true; set -u
SING=$(command -v singularity || command -v apptainer) || { echo "singularity/apptainer が見つかりません"; exit 1; }
mkdir -p "$VEP_DATA" "$VEP_CACHE" "$VEP_PLUGINS"; cd "$VEP_DATA"

## 1) VEP コンテナ
#   ログインノードではスレッド数制限で mksquashfs が失敗することがあるため:
#   (a) mksquashfs を 1 スレッドで SIF 作成 -> (b) 失敗したら sandbox (ディレクトリ形式, 圧縮なし) で作成
export APPTAINER_CACHEDIR="$VEP_DATA/.apptainer_cache"  APPTAINER_TMPDIR="$VEP_DATA/.apptainer_tmp"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"     SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
IMG="docker://ensemblorg/ensembl-vep:release_${VEP_REL}.0"
if [ ! -s "$VEP_SIF" ] && [ ! -d "$VEP_SANDBOX" ]; then
  if ! "$SING" build --mksquashfs-args "-processors 1" "$VEP_SIF" "$IMG"; then
    echo "[06] SIF の作成に失敗 -> sandbox 形式で作成します"
    rm -f "$VEP_SIF"
    "$SING" build --sandbox "$VEP_SANDBOX" "$IMG"
  fi
fi

## 2) オフラインキャッシュ (GRCh38, gnomAD v4 頻度・ClinVar 等を含む)
CT="homo_sapiens_vep_${VEP_REL}_GRCh38.tar.gz"
if [ ! -d "$VEP_CACHE/homo_sapiens/${VEP_REL}_GRCh38" ]; then
  [ -s "$CT" ] || wget -c "https://ftp.ensembl.org/pub/release-${VEP_REL}/variation/indexed_vep_cache/$CT"
  tar -xzf "$CT" -C "$VEP_CACHE" && rm -f "$CT"
fi

## 3) FASTA (HGVS 表記に必要) -> bgzip + faidx
if [ ! -s "$VEP_FASTA.fai" ]; then
  wget -c -O tmp.fa.gz "https://ftp.ensembl.org/pub/release-${VEP_REL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
  zcat tmp.fa.gz | bgzip -@4 > "$VEP_FASTA" && rm -f tmp.fa.gz
  samtools faidx "$VEP_FASTA"
fi

## 4) VEP プラグイン (公式リポジトリ)
[ -d "$VEP_PLUGINS/.git" ] || git clone -b "release/${VEP_REL}" --depth 1 https://github.com/Ensembl/VEP_plugins.git "$VEP_PLUGINS"

## 5) AlphaMissense (自動)
if [ ! -s "$AM_FILE.tbi" ]; then
  wget -c -O "$AM_FILE" https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
  tabix -s 1 -b 2 -e 2 -f -S 1 "$AM_FILE"
fi

## 6) 手動で用意するもの (規約/ログインが必要)
cat <<MSG
[06] 自動準備は完了。以下は必要に応じて手動で置く (無ければ 07 で自動的にスキップ):
  REVEL    : https://sites.google.com/site/revelgenomics/ -> $REVEL_FILE (+ .tbi)
             (zip を展開し VEP REVEL プラグインの説明どおり GRCh38 用に整形・bgzip・tabix)
  SpliceAI : Illumina BaseSpace の precomputed scores -> $SPLICEAI_SNV / $SPLICEAI_INDEL (+ .tbi)
  CADD     : https://cadd.gs.washington.edu/download -> $CADD_SNV / $CADD_INDEL (+ .tbi)
  LOFTEE   : git clone -b grch38 https://github.com/konradjk/loftee.git $LOFTEE_DIR
             + human_ancestor.fa.gz / GERP bigwig / loftee.sql を $(dirname "$LOFTEE_ANC") に
MSG
