#!/usr/bin/env bash
# Ensembl VEP REST API で candidates.tsv のバリアントをアノテーション (キャッシュ不要)
#   - 200 変異ずつ POST (/vep/homo_sapiens/region), 失敗・混雑時は待って再試行
#   - JSON は out/vep/rest/ に保存し, 07c_parse_vep_rest.R で TSV に変換
# 出力: out/vep/candidates.vep.tsv  (08_merge_vep.R がそのまま読む)
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$OUTDIR/vep/rest"; cd "$OUTDIR/vep/rest"

## 1) 入力 (VCF 風: CHROM POS ID REF ALT . . .), ID = chrN_POS_REF_ALT, CHROM は chr なし
awk -F'\t' 'NR>1 { c=$3; sub(/^chr/,"",c); id=$3"_"$4"_"$5"_"$6
                   if (!(id in s)) { s[id]=1; print c, $4, id, $5, $6, ".", ".", "." } }' \
    "$OUTDIR/candidates.tsv" > input.txt
echo "[07b] variants: $(wc -l < input.txt)"
rm -f chunk_???; split -l 200 -d -a 3 input.txt chunk_

## 2) POST (済んだ chunk は飛ばす)
for f in chunk_???; do
  out="$f.json"; [ -s "$out" ] && { echo "[07b] $f done (skip)"; continue; }
  body=$(awk 'BEGIN{printf "{\"variants\":["} {printf "%s\"%s\"", (NR>1?",":""), $0} END{printf "]}"}' "$f")
  for i in 1 2 3 4 5; do
    code=$(curl -sS -o "$out.tmp" -D hdr.txt -w '%{http_code}' -X POST \
           "$VEP_REST_SERVER/vep/homo_sapiens/region?$VEP_REST_PARAMS" \
           -H 'Content-Type: application/json' -H 'Accept: application/json' -d "$body" || echo 000)
    if [ "$code" = 200 ]; then mv "$out.tmp" "$out"; echo "[07b] $f OK"; break; fi
    if [ "$code" = 400 ]; then echo "[07b] 400 Bad Request:"; head -c 500 "$out.tmp"; echo
         echo "-> config.sh の VEP_REST_PARAMS から未対応のオプションを外して再実行"; exit 1; fi
    wait=$(awk 'tolower($1)=="retry-after:"{print int($2)+1}' hdr.txt); wait=${wait:-$((i*15))}
    echo "[07b] $f HTTP $code -> retry $i in ${wait}s"; sleep "$wait"
  done
  [ -s "$out" ] || { echo "[07b] ERROR: $f failed"; exit 1; }
  sleep 1
done

## 3) JSON -> TSV
Rscript "$SCRIPTDIR/07c_parse_vep_rest.R" "$OUTDIR/vep/rest" "$OUTDIR/vep/candidates.vep.tsv"
