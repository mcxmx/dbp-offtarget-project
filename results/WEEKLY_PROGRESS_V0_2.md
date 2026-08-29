# Weekly Progress v0.2

本周对 DBP off-target benchmark 原型完成了一次 scientific audit 和 v0.2 修正。首先检查了 v0.1 的数据生成逻辑，发现 PDB complex 曾被误标为 specificity ground truth，并且使用“最长蛋白链/最长 DNA 链”作为 benchmark 选择规则。v0.2 将证据拆分为 structural cognate、direct DNA-binding evidence、sequence-specificity evidence 和 quantitative specificity ground truth；当前 16 个 PDB pair 均不再被视为 quantitative ground truth。

随后对 16 个 PDB pair 逐条完成机制分类，保留 8 个 core benchmark pair；guide-dependent、lesion-specific、non-specific 和 transposase/substrate 样本被保留在 curated all table，但不进入 core specificity benchmark。基于 curated core 重新生成 benchmark_v0_2，共 25761 行，其中 single mutants 393 条、double mutants 9360 条、random negatives 16000 条。

分析方面，已重新输出 v0.2 figures，并完成 sequence-only proxy positional bias 检查，确认 k-mer/combined proxy 会产生位置效应，不能解释为生物学 specificity landscape。另建立了 Layer C experimental specificity pilot，从 JASPAR CORE 和 UniProt 获取 5 个 protein、1209 条 PFM-derived k-mer score 记录。下周重点是继续接入 raw PBM/HT-SELEX/CIS-BP 数据，并评估 protein-conditioned scoring interface。
