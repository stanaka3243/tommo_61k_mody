#!/usr/bin/env bash
# 02 で作った panel / off の VCF を VEP でアノテーション
#   - MANE Select を優先して 1 遺伝子 1 行 (--per_gene)
#   - HGVS (c./p.), gnomAD v4 (exome/genome, EAS), SIFT/PolyPhen, ClinVar(キャッシュ)
#   - プラグイン: AlphaMissense / REVEL / SpliceAI / CADD / LOFTEE (データがあるものだけ)
# 出力: out/vep/{panel,off}.vep.vcf.gz と {panel,off}.vep.tsv (ID = CHROM_POS_REF_ALT で 03 の結果と結合可)
set -euo pipefail
source "$(dirname "$0")/config.sh"
set +u; module load apptainer 2>/dev/null || true; set -u
SING=$(command -v singularity || command -v apptainer) || { echo "singularity/apptainer が見つかりません"; exit 1; }
mkdir -p "$OUTDIR/vep"; cd "$OUTDIR/vep"

## プラグイン引数 (データがあるものだけ)
PLUG=(); FIELDS_EXTRA=()
if [ -s "$AM_FILE" ]; then
  PLUG+=(--plugin "AlphaMissense,file=$AM_FILE"); FIELDS_EXTRA+=(am_class am_pathogenicity)
fi
if [ -s "$REVEL_FILE" ]; then
  PLUG+=(--plugin "REVEL,file=$REVEL_FILE"); FIELDS_EXTRA+=(REVEL)
fi
if [ -s "$SPLICEAI_SNV" ] && [ -s "$SPLICEAI_INDEL" ]; then
  PLUG+=(--plugin "SpliceAI,snv=$SPLICEAI_SNV,indel=$SPLICEAI_INDEL")
  FIELDS_EXTRA+=(SpliceAI_pred_DS_AG SpliceAI_pred_DS_AL SpliceAI_pred_DS_DG SpliceAI_pred_DS_DL)
fi
if [ -s "$CADD_SNV" ] && [ -s "$CADD_INDEL" ]; then
  PLUG+=(--plugin "CADD,snv=$CADD_SNV,indels=$CADD_INDEL"); FIELDS_EXTRA+=(CADD_PHRED)
fi
if [ -d "$LOFTEE_DIR" ] && [ -s "$LOFTEE_ANC" ]; then
  # VEP の --dir_plugins は 1 か所だけなので, 他のプラグインを LOFTEE ディレクトリにリンクする
  for pm in AlphaMissense REVEL SpliceAI CADD; do ln -sf "$VEP_PLUGINS/$pm.pm" "$LOFTEE_DIR/$pm.pm"; done
  PLUG+=(--dir_plugins "$LOFTEE_DIR"
         --plugin "LoF,loftee_path:$LOFTEE_DIR,human_ancestor_fa:$LOFTEE_ANC,conservation_file:$LOFTEE_SQL,gerp_bigwig:$LOFTEE_GERP")
  FIELDS_EXTRA+=(LoF LoF_filter LoF_flags)
else
  PLUG+=(--dir_plugins /plugins)
fi
echo "[07] plugins: ${FIELDS_EXTRA[*]:-(none)}"

FIELDS="Uploaded_variation,Location,Allele,SYMBOL,Gene,Feature,MANE_SELECT,CANONICAL,BIOTYPE,Consequence,IMPACT,HGVSc,HGVSp,EXON,INTRON,Protein_position,Amino_acids,Codons,SIFT,PolyPhen,gnomADe_AF,gnomADe_EAS_AF,gnomADg_AF,gnomADg_EAS_AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,Existing_variation"
[ ${#FIELDS_EXTRA[@]} -gt 0 ] && FIELDS="$FIELDS,$(IFS=,; echo "${FIELDS_EXTRA[*]}")"

BIND="$VEP_DATA,$OUTDIR,$WORKDIR,$VEP_PLUGINS:/plugins"
[ -d "$LOFTEE_DIR" ] && BIND="$BIND,$LOFTEE_DIR"

for name in panel off; do
  IN="$OUTDIR/$name.annot.vcf.gz"
  [ -s "$IN" ] || { echo "[07] skip $name (no $IN)"; continue; }
  echo "[07] VEP: $name ($(bcftools view -H "$IN" | wc -l) variants)"
  # ID を CHROM_POS_REF_ALT にし, INFO は軽くして VEP に渡す
  bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -x INFO "$IN" -Oz -o "$name.in.vcf.gz"
  bcftools index -t -f "$name.in.vcf.gz"

  common=(--offline --cache --dir_cache "$VEP_CACHE" --cache_version "$VEP_REL" --assembly GRCh38
          --fasta "$VEP_FASTA" --fork "$VEP_FORK" --force_overwrite --no_stats
          --hgvs --symbol --canonical --mane --biotype --numbers --sift b --polyphen b
          --af_gnomade --af_gnomadg --max_af --check_existing
          --per_gene --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length
          "${PLUG[@]}")

  "$SING" exec -B "$BIND" "$VEP_SIF" vep -i "$OUTDIR/vep/$name.in.vcf.gz" \
      "${common[@]}" --vcf --compress_output bgzip -o "$OUTDIR/vep/$name.vep.vcf.gz"
  "$SING" exec -B "$BIND" "$VEP_SIF" vep -i "$OUTDIR/vep/$name.in.vcf.gz" \
      "${common[@]}" --tab --fields "$FIELDS" -o "$OUTDIR/vep/$name.vep.tsv"
  bcftools index -t -f "$name.vep.vcf.gz"
done
echo "[07] done: $OUTDIR/vep/"
