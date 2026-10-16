# jMorp 61KJPN を使った MODY 保因者頻度の探索的解析 (SHIROKANE 用, bash + R)

## 何をするか
Sharp et al. (JCEM 2026) の UK Biobank 解析に近い方法で、jMorp の公開データ
(61KJPN-SNV/INDEL, off-panel singleton, JCNVv1) から MODY 10 遺伝子の保因者頻度を推定する。
データ全体 (約100 GB) はダウンロードせず、HTTPS + tabix index で対象領域だけを読む。

## 準備 (初回のみ)
    micromamba create -f env.yml      # または conda env create -f env.yml
    micromamba activate jmorp_mody

## 実行 (インターネットに出られるノードで)
    bash run_all.sh
個別に動かす場合は 01 → 02 → 03 → 04 → 05 の順。設定は config.sh。抽出領域は MODY_BED (/home/tanaka-s/MODY/mody_exons_merged.bed) を PAD (既定 10 bp) だけ広げて使う。
genes.tsv は遺伝子名リスト (MANE Select の選択と集計に使用; BED の遺伝子と一致させること)。

| スクリプト | 内容 | 主な出力 |
|---|---|---|
| 01_prepare_refs.sh | BED 整形 (chr 付与/PAD/マージ), 参照配列(対象染色体), Ensembl GFF3, MANE, ClinVar | ref/targets.chr.bed ほか |
| 02_extract_annotate.sh | 61KJPN と off-panel singleton から領域取得 → 正規化 → ClinVar → bcftools csq | out/panel.tsv, out/off.tsv |
| 03_classify.R | MANE Select で分類し保因者数・頻度 (1 in X, 95%CI) を集計 | out/candidates.tsv, review_missense.tsv, summary_61KJPN.tsv |
| 04_cnv_17q12.sh | JCNVv1 から 17q12 周辺の bin を抽出 | out/cnv_17q12.tsv |
| 05_cnv_17q12.R | HNF1B 上の CN1 (欠失) 人数と図 | out/cnv_17q12_bins.tsv, cnv_17q12.png |

## VEP アノテーション (任意, 06 → 07)
| スクリプト | 内容 | 主な出力 |
|---|---|---|
| 06_setup_vep.sh | 初回のみ: VEP コンテナ, キャッシュ(GRCh38), FASTA, プラグイン, AlphaMissense | $VEP_DATA |
| 07_run_vep.sh | panel/off の VCF を VEP で注釈 (MANE 優先, HGVS, gnomAD v4 EAS, SIFT/PolyPhen, AlphaMissense; REVEL/SpliceAI/CADD/LOFTEE はデータがあれば) | out/vep/*.vep.tsv |
| 08_merge_vep.R | candidates.tsv と VEP 結果を結合 (Rscript 08_merge_vep.R out) | out/candidates_vep.tsv, review_missense_vep.tsv, plp_core4_vep.tsv |
- 結合キーは CHROM_POS_REF_ALT と遺伝子名。AF/MAF 列 (AC/AN), gnomAD_absent (PM2 の目安) も追加される。
- REVEL, SpliceAI, CADD, LOFTEE のデータは規約・ログインの都合で手動準備 (06 の最後に案内が出る)。

## 分類の定義 (03_classify.R)
- ClinVar_PLP : ClinVar で Pathogenic / Likely pathogenic (Conflicting は除外)
- PTV         : stop_gained, frameshift, splice_acceptor/donor, start_lost (MANE Select 上)
- ultra_rare_missense : missense かつ AC <= 2 → review_missense.tsv を ACMG/ClinGen MDEP で手作業分類
- summary_61KJPN.tsv は A (ClinVar のみ), B (+PTV), C (+全 ultra-rare missense = 上限) の3段階

## 注意点
- 分母: 61KJPN 本体は 61,424 人。off-panel singleton は 100KJPN (AN ~ 20万) 由来なので別集計。
- jMorp の VCF の CLNSIG 欄は空のことがあるため、ClinVar は自分で付与している (CV_* 列)。
- 必ず MANE Select 転写産物で判定する (GCK の肝型エクソン、HNF1B の別アイソフォームで偽陽性が出る)。
- AN が低い部位 (low_AN = TRUE) は呼び出しの質に注意。
- 頻度の合計は近似 (同一人物が複数変異を持つ場合を区別できない)。
- jMorp の利用規約 (https://jmorp.megabank.tohoku.ac.jp/help/conditions-of-use) を確認し、論文では引用元を明記する。
- MANE/Ensembl/ClinVar のバージョンは config.sh と ref/clinvar_version.txt で記録される。

## このフォルダの古いファイル
jmorp_mody_pilot.py / candidates.tsv / *.tbi は Python での試行版。新しい解析では使わない (削除可)。
