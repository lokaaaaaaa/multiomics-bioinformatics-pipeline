# ============================================================
# Proteomics Workflow: Label-Free Quantification (LFQ) Analysis
# Tools: pandas, scipy, statsmodels, matplotlib, seaborn
# Compatible with MaxQuant proteinGroups.txt output
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import warnings
import os

warnings.filterwarnings("ignore")
os.makedirs("results", exist_ok=True)
os.makedirs("plots",   exist_ok=True)

# ── 1. Load / Simulate Data ────────────────────────────────────
# In real use: df = pd.read_csv("data/proteinGroups.txt", sep="\t")
# and filter/select LFQ intensity columns

np.random.seed(42)
n_proteins = 1500
n_ctrl     = 3
n_treat    = 3

ctrl_cols  = [f"LFQ_Control_{i}"   for i in range(1, n_ctrl  + 1)]
treat_cols = [f"LFQ_Treatment_{i}" for i in range(1, n_treat + 1)]
all_cols   = ctrl_cols + treat_cols

# Simulate LFQ intensities (log-normal distribution ~ real proteomics)
base_intensity = np.random.lognormal(mean=20, sigma=2, size=n_proteins)
ctrl_data  = np.column_stack([base_intensity * np.random.lognormal(0, 0.2, n_proteins)
                               for _ in range(n_ctrl)])
treat_data = np.column_stack([base_intensity * np.random.lognormal(0.5, 0.2, n_proteins)
                               for _ in range(n_treat)])

# Introduce ~10% truly up-regulated proteins
up_idx   = np.random.choice(n_proteins, size=150, replace=False)
treat_data[up_idx] *= np.random.uniform(1.5, 4, size=(150, n_treat))

df = pd.DataFrame(
    np.hstack([ctrl_data, treat_data]),
    columns=all_cols
)
df.insert(0, "Protein_ID", [f"PROT_{i:04d}" for i in range(n_proteins)])
df.insert(1, "Gene_names", [f"GENE{i}"      for i in range(n_proteins)])

# Introduce ~15% missing values (MCAR)
mask = np.random.rand(*df[all_cols].shape) < 0.15
df[all_cols] = df[all_cols].where(~mask, other=np.nan)

print(f"Loaded {len(df)} protein groups, {df[all_cols].isna().sum().sum()} missing values")

# ── 2. Pre-processing ──────────────────────────────────────────
## 2a. Log2 transform
df_log = df.copy()
df_log[all_cols] = np.log2(df[all_cols].replace(0, np.nan))

## 2b. Filter: keep proteins with ≥ 2 valid values per group
valid_ctrl  = df_log[ctrl_cols].notna().sum(axis=1)  >= 2
valid_treat = df_log[treat_cols].notna().sum(axis=1) >= 2
df_log = df_log[valid_ctrl & valid_treat].reset_index(drop=True)
print(f"After filtering: {len(df_log)} proteins retained")

## 2c. Median normalisation
medians = df_log[all_cols].median(skipna=True)
global_med = medians.mean()
df_log[all_cols] = df_log[all_cols].subtract(medians - global_med)

## 2d. KNN-style imputation (column-wise median for simplicity)
imputer = SimpleImputer(strategy="median")
df_log[all_cols] = imputer.fit_transform(df_log[all_cols])

# ── 3. Differential Abundance Analysis ────────────────────────
results = []
for _, row in df_log.iterrows():
    ctrl_vals  = row[ctrl_cols].values.astype(float)
    treat_vals = row[treat_cols].values.astype(float)
    fc         = np.mean(treat_vals) - np.mean(ctrl_vals)   # log2 FC
    t_stat, p  = stats.ttest_ind(treat_vals, ctrl_vals, equal_var=False)
    results.append({
        "Protein_ID":       row["Protein_ID"],
        "Gene_names":       row["Gene_names"],
        "log2FC":           fc,
        "p_value":          p,
        "mean_ctrl":        np.mean(ctrl_vals),
        "mean_treat":       np.mean(treat_vals)
    })

res_df = pd.DataFrame(results)

# Multiple-testing correction (BH/FDR)
reject, p_adj, _, _ = multipletests(res_df["p_value"].fillna(1), method="fdr_bh")
res_df["padj"]       = p_adj
res_df["significant"] = (res_df["padj"] < 0.05) & (res_df["log2FC"].abs() > 1)

sig = res_df[res_df["significant"]].sort_values("padj")
print(f"Significantly changed proteins: {len(sig)}")
print(f"  Up-regulated:   {(sig['log2FC'] > 1).sum()}")
print(f"  Down-regulated: {(sig['log2FC'] < -1).sum()}")

res_df.to_csv("results/proteomics_DA_results.csv", index=False)
sig.to_csv("results/proteomics_significant.csv",   index=False)

# ── 4. Visualisations ─────────────────────────────────────────
## 4a. Volcano plot
fig, ax = plt.subplots(figsize=(9, 7))
colors = np.where(
    (res_df["log2FC"] > 1)  & (res_df["padj"] < 0.05), "#e63946",
    np.where(
    (res_df["log2FC"] < -1) & (res_df["padj"] < 0.05), "#457b9d",
    "#adb5bd"))
ax.scatter(res_df["log2FC"], -np.log10(res_df["padj"] + 1e-300),
           c=colors, s=12, alpha=0.7, linewidths=0)
ax.axhline(-np.log10(0.05), color="black", lw=1, ls="--", alpha=0.6)
ax.axvline(1,  color="black", lw=1, ls="--", alpha=0.6)
ax.axvline(-1, color="black", lw=1, ls="--", alpha=0.6)
ax.set_xlabel("log₂ Fold Change (Treatment / Control)", fontsize=12)
ax.set_ylabel("-log₁₀ Adjusted p-value",               fontsize=12)
ax.set_title("Proteomics Volcano Plot",                 fontsize=14, fontweight="bold")
patches = [mpatches.Patch(color="#e63946", label=f"Up ({(sig['log2FC']>1).sum()})"),
           mpatches.Patch(color="#457b9d", label=f"Down ({(sig['log2FC']<-1).sum()})"),
           mpatches.Patch(color="#adb5bd", label="NS")]
ax.legend(handles=patches, framealpha=0.8)
plt.tight_layout()
plt.savefig("plots/proteomics_volcano.png", dpi=150)
plt.close()

## 4b. PCA of samples
scaler  = StandardScaler()
X_scaled = scaler.fit_transform(df_log[all_cols].T)      # samples × proteins
pca     = PCA(n_components=2)
pcs     = pca.fit_transform(X_scaled)
colors_pca = ["#e63946"] * n_ctrl + ["#2a9d8f"] * n_treat
labels_pca = ctrl_cols + treat_cols

fig, ax = plt.subplots(figsize=(7, 6))
for i, (x, y, lab, col) in enumerate(zip(pcs[:, 0], pcs[:, 1], labels_pca, colors_pca)):
    ax.scatter(x, y, color=col, s=80, zorder=3)
    ax.annotate(lab, (x, y), textcoords="offset points", xytext=(6, 4), fontsize=8)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("Sample PCA — Proteomics", fontweight="bold")
ctrl_patch  = mpatches.Patch(color="#e63946", label="Control")
treat_patch = mpatches.Patch(color="#2a9d8f", label="Treatment")
ax.legend(handles=[ctrl_patch, treat_patch])
plt.tight_layout()
plt.savefig("plots/proteomics_PCA.png", dpi=150)
plt.close()

## 4c. Heatmap of top 50 DA proteins
top50_ids = sig.head(50)["Protein_ID"].tolist()
heat_data = df_log[df_log["Protein_ID"].isin(top50_ids)].set_index("Gene_names")[all_cols]
fig, ax   = plt.subplots(figsize=(10, 12))
sns.heatmap(heat_data, cmap="RdBu_r", center=0, ax=ax,
            xticklabels=True, yticklabels=True, linewidths=0.2,
            cbar_kws={"label": "log₂ Intensity (normalised)"})
ax.set_title("Top 50 Differentially Abundant Proteins", fontweight="bold", fontsize=13)
plt.tight_layout()
plt.savefig("plots/proteomics_heatmap.png", dpi=150)
plt.close()

print("\n✅ Proteomics pipeline complete — results/ and plots/ updated.")
