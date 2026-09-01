# v0.3 Progress Report

本轮建立了 designed DBP experimental specificity benchmark v0.3，核心数据源为 GEO GSE237017。已程序化解析 series/GSM metadata，并下载 12 个 usable uPBM samples，覆盖 7 个 designed DBP：DBP1、DBP3、DBP5、DBP6、DBP9、DBP35、DBP48。每个 sample 的 processed 7-mer 文件原始为 8192 行、两个 7-mer 列；解析时显式展开 reverse-complement companion column 后，每个 sample 覆盖 16,384 个 unique 7-mers，缺失数为 0。最终 consensus benchmark 包含 114688 条 protein-7mer experimental measurements。

Replicate QC 显示 DBP1、DBP3、DBP6、DBP9、DBP35 有 replicate；DBP5 和 DBP48 为 single replicate only。E-score replicate Pearson 范围为 0.551-0.765，中位数 0.665；Spearman 范围为 0.469-0.734，中位数 0.591。7 个 DBP 的 protein sequence 和 intended target sequence 均已从官方 supplementary workbook 恢复，confidence 为 high。

Target rank 分析只基于 intended target 的 overlapping 7-mers，不解释为 full-target affinity。多数设计的 best target-derived 7-mer 位于较高 percentile，但 DBP48 的 best target-derived 7-mer percentile 为 0.773，相对较弱。Sequence-only baseline 与 PBM E-score 的 per-protein Spearman 整体较低，中位数接近 0；同时发现 140 个 top 1% E-score 但 target-derived 7-mer similarity 不高的候选，说明仅靠 DNA sequence similarity 难以解释 designed DBP specificity landscape。

当前最大限制是 PBM processed score 是 7-mer 级别，不能直接代表完整 target DNA 的 binding affinity；DBP5/DBP48 缺少 replicate；E-score 也不能跨 protein 当作绝对 affinity。下一步应进入 v0.4 protein-conditioned baseline：以 protein sequence、intended target-derived context 和 candidate 7-mer 为输入，先做非神经或轻量模型的 per-protein ranking baseline，再和 sequence-only baseline 比较。
