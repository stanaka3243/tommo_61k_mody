#!/usr/bin/env bash
# ===== 共通設定 (各スクリプトから source される) =====
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$SCRIPTDIR}"          # 作業ディレクトリ (変更可)
REFDIR="$WORKDIR/ref"
OUTDIR="$WORKDIR/out"
GENES="$SCRIPTDIR/genes.tsv"              # 遺伝子名リスト (MANE の選択と集計に使う)
MODY_BED="/home/tanaka-s/MODY/mody_exons_merged.bed"   # 抽出領域 (GRCh38, BED 0-based)
PAD=10                                     # BED の各区間の前後に足す bp (splice donor/acceptor ±2 を含めるため)

# --- 参照データのバージョン (必要に応じて更新) ---
ENSEMBL_REL=113                            # Ensembl GFF3 (bcftools csq 用)
MANE_VER=1.4                               # NCBI MANE summary

# --- jMorp のファイル (2026-10 時点の URL) ---
JMORP="https://jmorp.megabank.tohoku.ac.jp/datasets"
AF_URL="$JMORP/tommo-61kjpn-20250616-af_snvindelall/files/tommo-61kjpn-20250616-GRCh38-snvindel-af-autosome-{CHR}.vcf.gz"
OFF_URL="$JMORP/tommo-61kjpn-20260831-off_panel_singleton/files/tommo-61kjpn-20260831-GRCh38-snvindel-af-autosome-chr1_22-off_panel_singleton.vcf.gz"
CNV_URL="$JMORP/tommo-jcnvv1-20230828/files/tommo-jcnvv1-20230828-GRCh38.vcf.gz"

# --- 解析パラメータ ---
MAX_AC_MISSENSE=2      # ultra-rare missense の上限 (Sharp et al. の gnomAD AC<=2 に準拠)
N_PANEL=61424          # 61KJPN の人数 (分母)
