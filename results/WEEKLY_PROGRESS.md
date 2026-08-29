# Weekly Progress

本周完成：
本周完成了一个可复现的 DBP-DNA off-target prototype。首先基于 RCSB PDB 整理了 16 条真实 protein-DNA 配对，并保留了结构条目、链 ID、PDB ID、来源 URL 和检索日期。
随后构建了 single mutant 数据集 843 条、double mutant 数据集 23904 条，以及 GC-matched/random negatives 共 32000 条。
整合后的 benchmark_v0.1 共 56763 条，已生成 sequence-only proxy baseline 和 preliminary figures。
genome scan demo 以 GRCh38 chr22 为目标，共返回 200 个候选位点。

下周计划：
下周重点是扩充具有 quantitative specificity ground truth 的数据源，并把 protein-conditioned scoring 模型接口接到现有 benchmark 上。

限制：当前结果仍是 sequence-only proxy，尚不能替代 protein-conditioned binding prediction。