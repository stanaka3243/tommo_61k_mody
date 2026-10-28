#!/usr/bin/env bash
# ===== 実行場所 =====
#   slade     : Exeter Slade。61KJPN はローカルの VCF を読む (AF_LOCAL)
#   shirokane : SHIROKANE。61KJPN は jMorp から HTTPS で領域だけ読む (AF_URL)
DATA_SITE="${DATA_SITE:-slade}"

# ===== SHIROKANE: environment modules (最初に読み込む) =====
if [ "$DATA_SITE" = shirokane ]; then
  if ! command -v module >/dev/null 2>&1; then          # 非対話シェルでは module 関数が未定義のことがある
    for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do
      [ -f "$f" ] && source "$f" && break
    done
  fi
  set +u                                               # module スクリプトは未定義変数を参照することがある
  module use /usr/local/package/modulefiles
  module load samtools
  set -u
fi
# ===== 共通設定 (各スクリプトから source される) =====
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$SCRIPTDIR}"          # 作業ディレクトリ (変更可)
REFDIR="$WORKDIR/ref"
OUTDIR="$WORKDIR/out"
GENES="${GENES:-$SCRIPTDIR/genes_core4.tsv}"   # 対象遺伝子 (MODY1/2/3/5)。10 遺伝子なら genes.tsv

# --- 抽出領域 (01 で MANE Select 転写産物から作る) ---
UPSTREAM=5000          # 転写産物の 5' 側に足す bp (HNF4A P1 プロモーターなどを含める)
DOWNSTREAM=1000        # 転写産物の 3' 側に足す bp
PAD=10                 # CDS の前後 PAD bp を REGION=CDS_pm{PAD} とする (splice donor/acceptor ±2 を含む)

# --- 参照データのバージョン (必要に応じて更新) ---
ENSEMBL_REL=113                            # Ensembl GFF3 (bcftools csq 用)
MANE_VER=1.4                               # NCBI MANE summary

# --- jMorp のファイル ---
AF_LOCAL_DEFAULT='/slade/home/st1041/tommo/61kjpn/tommo-61kjpn-20250616-GRCh38-snvindel-af-autosome-{CHR}.vcf.gz'
AF_LOCAL="${AF_LOCAL:-$AF_LOCAL_DEFAULT}"   # Slade 上の 61KJPN ({CHR} は chr7 などに置換される)
OFF_LOCAL="${OFF_LOCAL:-}"                 # Slade 上の off-panel singleton (空なら OFF_URL を HTTPS で読む)
JMORP="https://jmorp.megabank.tohoku.ac.jp/datasets"
AF_URL="$JMORP/tommo-61kjpn-20250616-af_snvindelall/files/tommo-61kjpn-20250616-GRCh38-snvindel-af-autosome-{CHR}.vcf.gz"
OFF_URL="$JMORP/tommo-61kjpn-20260831-off_panel_singleton/files/tommo-61kjpn-20260831-GRCh38-snvindel-af-autosome-chr1_22-off_panel_singleton.vcf.gz"
CNV_URL="$JMORP/tommo-jcnvv1-20230828/files/tommo-jcnvv1-20230828-GRCh38.vcf.gz"

# --- 解析パラメータ ---
MAX_AC_MISSENSE=2      # ultra-rare missense の上限 (Sharp et al. の gnomAD AC<=2 に準拠)
N_PANEL=61424          # 61KJPN の人数 (分母)
# 03 の候補 (従来の集計) に使う条件。all_variants.tsv はこの条件に関係なく全件出る
CAND_PASS_ONLY=1               # 1 = FILTER=PASS のみ, 0 = FILTER を問わない
CAND_REGIONS="CDS_pm${PAD}"    # カンマ区切り: CDS_pm10, UTR_exon, intron, upstream, downstream

# ===== VEP (06_setup_vep.sh / 07_run_vep.sh) =====
VEP_REL=113                                        # Ensembl/VEP リリース (GFF3 と揃える)
VEP_DATA="${VEP_DATA:-$HOME/vep_data}"             # キャッシュ・プラグイン・データの置き場 (数十 GB)
VEP_SIF="$VEP_DATA/ensembl-vep_${VEP_REL}.sif"     # Apptainer イメージ (SIF)
VEP_SANDBOX="$VEP_DATA/ensembl-vep_${VEP_REL}_sandbox"   # SIF が作れない場合の sandbox
VEP_CACHE="$VEP_DATA/cache"
VEP_PLUGINS="$VEP_DATA/plugins"
VEP_FASTA="$VEP_DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"   # HGVS 用 (bgzip)
VEP_FORK=4
# --- プラグイン用データ (ファイルがあれば自動で使う / 無ければスキップ) ---
AM_FILE="$VEP_DATA/AlphaMissense_hg38.tsv.gz"                      # 06 で自動ダウンロード
REVEL_FILE="$VEP_DATA/revel/new_tabbed_revel_grch38.tsv.gz"        # 手動 (学術利用規約あり)
SPLICEAI_SNV="$VEP_DATA/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz"     # 手動 (Illumina BaseSpace)
SPLICEAI_INDEL="$VEP_DATA/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz"
CADD_SNV="$VEP_DATA/cadd/whole_genome_SNVs.tsv.gz"                 # 任意 (約 80 GB)
CADD_INDEL="$VEP_DATA/cadd/gnomad.genomes.r4.0.indel.tsv.gz"
LOFTEE_DIR="$VEP_DATA/loftee"                                      # 任意 (GRCh38 ブランチ)
LOFTEE_ANC="$VEP_DATA/loftee_data/human_ancestor.fa.gz"
LOFTEE_GERP="$VEP_DATA/loftee_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
LOFTEE_SQL="$VEP_DATA/loftee_data/loftee.sql"

# ===== VEP REST API (07b_vep_rest.sh) =====
VEP_REST_SERVER="https://rest.ensembl.org"
# 使うオプション (未対応のものがあると 400 エラー -> ここから外す)
VEP_REST_PARAMS="hgvs=1&mane=1&canonical=1&numbers=1&variant_class=1&af_gnomade=1&af_gnomadg=1&AlphaMissense=1&REVEL=1&SpliceAI=1&CADD=1&LoF=1"
