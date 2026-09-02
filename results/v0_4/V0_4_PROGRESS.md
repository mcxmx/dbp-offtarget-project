# v0.4 进展报告

本轮在 v0.3.1 已通过 validation 的 GSE237017 designed-DBP uPBM benchmark 上，建立了 strong baseline arena 的第一版。v0.3.1 数据保持冻结，v0.4 新结果单独保存在 `data/processed/v0_4/`、`metadata/v0_4/` 和 `results/v0_4/`。

已完成 DeepPBS、NA-MPNN 和 Tier 1 SimpleProteinConditionalBaseline 的可行性审计。DeepPBS 本轮未公平运行，原因是官方预处理依赖 Linux/结构特征工具链和额外图神经网络依赖；SimpleProteinConditionalBaseline 只保留为未训练的 protein-conditioned 接口，因为尚未加入 assay-matched natural PBM/uPBM 训练集。NA-MPNN 使用官方 specificity checkpoint 成功对 DBP35 和 DBP48 产生结构 PPM 诊断预测，但 DBP48/8TAC 在 NA-MPNN split 文件中出现，不能作为 zero-shot 结果。

当前评估覆盖 57344 个 protein-RC-class 单位。最好的 sequence-only baseline 是 `sequence_kmer3`，macro median Spearman 为 0.232。NA-MPNN 诊断结果覆盖 2/7 个 DBP，macro median Spearman 为 -0.041。五个暂未能结构评估的 DBP 是：DBP1, DBP3, DBP5, DBP6, DBP9。v0.3.1 的 1,515 个 sequence-vs-experiment disagreement candidates 中，NA-MPNN 可评估 398 个，按预设 top-10% 规则解析 50 个。

本轮生成了 baseline performance、per-protein heatmap、prediction-vs-experiment、replicate reference、motif-distance 分层和 failure landscape 六张图。最终 gate 为 `CONDITIONAL GO`：sequence-only gap 明显，但 DeepPBS/NA-MPNN 的全覆盖强 baseline 证据仍不足。下一步最应优先补充 assay-matched natural PBM/uPBM training control，并在可复现环境中补跑 DeepPBS 或替代结构-aware baseline。
