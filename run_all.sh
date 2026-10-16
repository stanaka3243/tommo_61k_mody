#!/usr/bin/env bash
# まとめて実行:  bash run_all.sh   (インターネット接続のあるノードで)
set -euo pipefail
source "$(dirname "$0")/config.sh"
bash "$SCRIPTDIR/01_prepare_refs.sh"
bash "$SCRIPTDIR/02_extract_annotate.sh"
Rscript "$SCRIPTDIR/03_classify.R" "$OUTDIR" "$REFDIR" "$GENES" "$MAX_AC_MISSENSE" "$N_PANEL"
bash "$SCRIPTDIR/04_cnv_17q12.sh"
Rscript "$SCRIPTDIR/05_cnv_17q12.R" "$OUTDIR"
