"""Render Figure 4 from frozen v1.3/v1.3.2/v1.3.3 result tables."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results/v1_3"
FIGURES = ROOT / "figures/v1_3"


def build_source() -> pd.DataFrame:
    rows = []
    dev = pd.read_csv(RESULTS / "figure2_data.tsv", sep="\t")
    dev = dev[(dev.method == "DeepPBS") & (dev.analysis_status == "PRIMARY")]
    for _, r in dev.iterrows():
        rows.append({"cohort": "Designed DBP development", "status": "MODEL PERFORMANCE", "unit": r.protein_id, "global": r.global_spearman, "position": r.position_spearman, "identity": r.within_position_residual_spearman, "pairwise": r.within_position_pairwise_accuracy})
    hold = pd.read_csv(RESULTS / "holdout_decomposition.tsv", sep="\t")
    hold = hold[hold.row_type == "protein"]
    for _, r in hold.iterrows():
        rows.append({"cohort": "Designed DBP prospective holdout", "status": "PROSPECTIVE HOLDOUT", "unit": r.protein_id, "global": r.global_spearman, "position": r.position_spearman, "identity": r.within_residual_spearman, "pairwise": r.within_pairwise_accuracy})
    samba = pd.read_csv(RESULTS / "external_samba_decomposition_per_tf.tsv", sep="\t")
    samba = samba[samba.perturbation_type == "canonical_watson_crick"]
    for _, r in samba.iterrows():
        rows.append({"cohort": "SaMBA technical repeatability", "status": "TECHNICAL REPEATABILITY", "unit": r.tf_id, "global": r.global_spearman, "position": r.position_spearman, "identity": r.identity_residual_spearman, "pairwise": r.identity_pairwise_accuracy})
    natural = pd.read_csv(RESULTS / "natural_tf_decomposition_per_tf.tsv", sep="\t")
    for _, r in natural.iterrows():
        rows.append({"cohort": "Natural TF DeepPBS", "status": "MODEL PERFORMANCE / NOT EVALUABLE", "unit": r.site_id, "global": r.global_spearman, "position": r.position_spearman, "identity": r.identity_residual_spearman, "pairwise": r.identity_pairwise_accuracy})
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "figure4_natural_vs_designed_identity.tsv", sep="\t", index=False, na_rep="NA")
    return out


def render(source: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    metrics = [("global", "Global"), ("position", "Position"), ("identity", "Identity residual"), ("pairwise", "Identity pairwise")]
    order = ["SaMBA technical repeatability", "Natural TF DeepPBS", "Designed DBP development", "Designed DBP prospective holdout"]
    colors = {"SaMBA technical repeatability": "#2F7D8C", "Natural TF DeepPBS": "#B8C0CC", "Designed DBP development": "#667085", "Designed DBP prospective holdout": "#C43D3D"}
    fig, axes = plt.subplots(1, 4, figsize=(13.8, 4.5), sharey=False)
    rng = np.random.default_rng(1301)
    for ax, (metric, title) in zip(axes, metrics):
        for x, cohort in enumerate(order):
            g = source[source.cohort == cohort]
            values = pd.to_numeric(g[metric], errors="coerce").dropna().to_numpy(float)
            if len(values):
                jitter = rng.uniform(-0.12, 0.12, len(values))
                ax.scatter(np.full(len(values), x) + jitter, values, s=20, alpha=0.65, color=colors[cohort], edgecolor="white", linewidth=0.3)
                ax.plot([x - 0.18, x + 0.18], [np.median(values)] * 2, color=colors[cohort], lw=3)
            else:
                ax.scatter([x], [0.05 if metric != "pairwise" else 0.5], marker="x", s=45, color="#98A2B3", linewidth=1.5)
                ax.text(x, 0.10 if metric != "pairwise" else 0.57, "NA", ha="center", va="bottom", fontsize=8, color="#667085")
        ax.set_title(title, fontsize=10)
        ax.set_xticks(range(len(order)), ["SaMBA\ntechnical", "Natural TF\nDeepPBS", "Designed\ndevelopment", "Designed\nholdout"], rotation=30, ha="right", fontsize=8)
        ax.grid(axis="y", color="#E4E7EC", lw=0.7)
        if metric == "pairwise":
            ax.axhline(0.5, color="black", lw=0.8, ls="--")
            ax.set_ylim(0.25, 1.0)
        else:
            ax.axhline(0, color="black", lw=0.8)
            ax.set_ylim(-0.25, 1.05)
    axes[0].set_ylabel("Correlation / pairwise accuracy")
    fig.suptitle("Natural versus designed identity-level challenge", y=1.03, fontsize=13)
    fig.text(0.5, -0.02, "SaMBA = technical within-assay repeatability; DeepPBS points = model performance; natural-TF primary challenge is not evaluable after frozen preprocessing failure", ha="center", fontsize=8.5, color="#475467")
    fig.tight_layout()
    fig.savefig(FIGURES / "figure4_natural_vs_designed_identity.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "figure4_natural_vs_designed_identity.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    render(build_source())
