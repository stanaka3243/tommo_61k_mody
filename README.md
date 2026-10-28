# jMorp 61KJPN を使った MODY 保因者頻度の探索的解析 (Slade / SHIROKANE, bash + R)

## 何をするか
Sharp et al. (JCEM 2026) の UK Biobank 解析に近い方法で、jMorp の公開データ
(61KJPN-SNV/INDEL, off-panel singleton, JCNVv1) から MODY 遺伝子の保因者頻度を推定する。
既定の対象は MODY1/2/3/5 の 4 遺伝子 (genes_core4.tsv; 10 遺伝子なら GENES=genes.tsv)。
抽出の段階ではバリアントを絞り込まず (FILTER・頻度・領域・影響度はすべて列として残す)、
絞り込みは 03 以降で列の条件として行う。
実行場所は config.sh の DATA_SITE で選ぶ: slade (既定) はローカルの 61KJPN VCF を読み、
shirokane はデータ全体 (約100 GB) をダウンロードせず HTTPS + tabix index で対象領域だけを読む。

## 準備 (初回のみ)
    micromamba create -f env.yml      # または conda env create -f env.yml
    micromamba activate jmorp_mody

## 実行 (インターネットに出られるノードで)
    bash run_all.sh
個別に動かす場合は 01 → 02 → 03 → 04 → 05 の順。設定は config.sh。
抽出領域は各遺伝子の MANE Select 転写産物の全長に、上流 UPSTREAM (既定 5 kb) と下流 DOWNSTREAM (既定 1 kb) を足したもの (鎖の向きを考慮)。
Slade では qrun で投入してよい (01 と、off-panel を HTTPS で読む場合の 02 はインターネット接続が必要)。

| スクリプト | 内容 | 主な出力 |
|---|---|---|
| 01_prepare_refs.sh | Ensembl GFF3, MANE, 抽出領域 (MANE 転写産物 + 上流/下流), 参照配列(対象染色体), ClinVar | ref/mane_tx.tsv, ref/features.tsv, ref/targets.chr.bed ほか |
| 02_extract_annotate.sh | 61KJPN と off-panel singleton から領域取得 (FILTER で絞らない) → 正規化 → ClinVar → bcftools csq | out/panel.tsv, out/off.tsv |
| 03_classify.R | 全バリアントに判定列を付与; 候補条件を満たすものを MANE Select で分類し保因者数・頻度 (1 in X, 95%CI) を集計 | out/all_variants.tsv, overview_counts.tsv, candidates.tsv, review_missense.tsv, summary_61KJPN.tsv |
| 04_cnv_17q12.sh | JCNVv1 から 17q12 周辺の bin を抽出 | out/cnv_17q12.tsv |
| 05_cnv_17q12.R | HNF1B 上の CN1 (欠失) 人数と図 | out/cnv_17q12_bins.tsv, cnv_17q12.png |

## VEP アノテーション (任意, 06 → 07)
| スクリプト | 内容 | 主な出力 |
|---|---|---|
| 06_setup_vep.sh | 初回のみ: VEP コンテナ, キャッシュ(GRCh38), FASTA, プラグイン, AlphaMissense | $VEP_DATA |
| 07_run_vep.sh | panel/off の VCF を VEP で注釈 (MANE 優先, HGVS, gnomAD v4 EAS, SIFT/PolyPhen, AlphaMissense; REVEL/SpliceAI/CADD/LOFTEE はデータがあれば) | out/vep/*.vep.tsv |
| 07b_vep_rest.sh | (キャッシュ不要) Ensembl VEP REST API で candidates.tsv を注釈。REVEL/SpliceAI/CADD/AlphaMissense/LOFTEE も取得 | out/vep/candidates.vep.tsv |
| 08_merge_vep.R | candidates.tsv と VEP 結果を結合 (Rscript 08_merge_vep.R out) | out/candidates_vep.tsv, review_missense_vep.tsv, plp_core4_vep.tsv |
- 結合キーは CHROM_POS_REF_ALT と遺伝子名。AF/MAF 列 (AC/AN), gnomAD_absent (PM2 の目安) も追加される。
- REVEL, SpliceAI, CADD, LOFTEE のデータは規約・ログインの都合で手動準備 (06 の最後に案内が出る)。

## 全バリアント表 (03_classify.R → all_variants.tsv)
- 抽出領域内のすべての変異アレル (多アレル座位は 1 ALT 1 行)。AC=0 や FILTER≠PASS も含む
- filter / pass : VCF の FILTER と PASS かどうか
- platform_bias : TOMMO_POSSIBLE_PLATFORM_BIAS (ヘッダーにある場合)
- region        : CDS_pm10 (CDS ±PAD bp), UTR_exon, intron, upstream, downstream (MANE Select 転写産物基準)
- consequence   : bcftools csq の MANE Select 上の結果 (イントロン・上下流は転写産物 ID が付かないので ".")
- clinvar_plp / class : 下記の定義。all_variants.tsv では絞り込みに使わず列として付けるだけ
- overview_counts.tsv : source × gene × region × PASS/nonPASS の件数

## 候補の条件 (config.sh)
- CAND_PASS_ONLY=1, CAND_REGIONS=CDS_pm10 が既定 (以前のエクソン BED ±10 bp・PASS のみに近い)
- 例: FILTER を問わず UTR も含める → CAND_PASS_ONLY=0, CAND_REGIONS=CDS_pm10,UTR_exon
- summary_61KJPN.tsv の candidate_filter 列に、使った条件が記録される

## 分類の定義 (03_classify.R)
- ClinVar_PLP : ClinVar で Pathogenic / Likely pathogenic (Conflicting は除外)
- PTV         : stop_gained, frameshift, splice_acceptor/donor, start_lost (MANE Select 上)
- ultra_rare_missense : missense かつ AC <= 2 → review_missense.tsv を ACMG/ClinGen MDEP で手作業分類
- summary_61KJPN.tsv は A (ClinVar のみ), B (+PTV), C (+全 ultra-rare missense = 上限) の3段階

## 注意点
- 分母: 61KJPN 本体は 61,424 人。off-panel singleton は 100KJPN (AN ~ 20万) 由来なので別集計。
- jMorp の VCF の CLNSIG 欄は空のことがあるため、ClinVar は自分で付与している (CV_* 列)。
- 必ず MANE Select 転写産物で判定する (GCK の肝型エクソン、HNF1B の別アイソフォームで偽陽性が出る)。
- MANE Select にないエクソン (別アイソフォーム固有) は REGION=intron などになる。
- HNF1B の 17q12 全欠失は SNV/INDEL パネルでは見えない (04/05 の JCNVv1 で別に確認)。
- AN が低い部位 (low_AN = TRUE) は呼び出しの質に注意。
- 頻度の合計は近似 (同一人物が複数変異を持つ場合を区別できない)。
- jMorp の利用規約 (https://jmorp.megabank.tohoku.ac.jp/help/conditions-of-use) を確認し、論文では引用元を明記する。
- MANE/Ensembl/ClinVar のバージョンは config.sh と ref/clinvar_version.txt で記録される。

## このフォルダの古いファイル
jmorp_mody_pilot.py / candidates.tsv / *.tbi は Python での試行版。新しい解析では使わない (削除可)。
