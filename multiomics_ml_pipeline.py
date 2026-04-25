# ============================================================
# AI/ML Multi-Omics Integration Pipeline
# Methods: MOFA+, Random Forest, Autoencoder, SHAP
# Tools: scikit-learn, torch, shap, mofapy2
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import warnings
import os

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (roc_auc_score, roc_curve, confusion_matrix,
                              classification_report, ConfusionMatrixDisplay)
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectFromModel
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

warnings.filterwarnings("ignore")
os.makedirs("results", exist_ok=True)
os.makedirs("plots",   exist_ok=True)

SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

# ── 1. Simulate Multi-Omics Data ──────────────────────────────
n_samples   = 120
n_ctrl      = 60
labels      = np.array([0] * n_ctrl + [1] * (n_samples - n_ctrl))

# Transcriptomics (500 genes)
rna = np.random.randn(n_samples, 500)
rna[n_ctrl:, :50] += 2.0      # up-regulated in cases

# Proteomics (300 proteins)
prot = np.random.randn(n_samples, 300)
prot[n_ctrl:, :30] += 1.5

# Methylation (200 CpG beta-values → M-values)
meth = np.random.randn(n_samples, 200)
meth[n_ctrl:, :20] -= 1.8

sample_ids = [f"S{i:03d}" for i in range(n_samples)]
print(f"Simulated data: {n_samples} samples | {n_ctrl} controls, {n_samples-n_ctrl} cases")
print(f"Omics layers: RNA ({rna.shape[1]}) | Protein ({prot.shape[1]}) | Methylation ({meth.shape[1]})")

# ── 2. Concatenated Feature Matrix ────────────────────────────
X_concat = np.hstack([rna, prot, meth])
feature_names = (
    [f"RNA_{i}"  for i in range(rna.shape[1])] +
    [f"PROT_{i}" for i in range(prot.shape[1])] +
    [f"METH_{i}" for i in range(meth.shape[1])]
)
X_df = pd.DataFrame(X_concat, columns=feature_names, index=sample_ids)
y    = pd.Series(labels, name="label", index=sample_ids)

# ── 3. Autoencoder for Latent Omics Representation ────────────
class OmicsAutoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim=32):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 256), nn.BatchNorm1d(256), nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 64),        nn.BatchNorm1d(64),  nn.ReLU(),
            nn.Linear(64, latent_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 64), nn.ReLU(),
            nn.Linear(64, 256),        nn.ReLU(),
            nn.Linear(256, input_dim)
        )
    def forward(self, x):
        z = self.encoder(x)
        return self.decoder(z), z

scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X_concat).astype(np.float32)
X_tensor = torch.tensor(X_scaled)
dataset  = TensorDataset(X_tensor, X_tensor)
loader   = DataLoader(dataset, batch_size=16, shuffle=True)

ae     = OmicsAutoencoder(input_dim=X_scaled.shape[1], latent_dim=32)
opt    = optim.Adam(ae.parameters(), lr=1e-3, weight_decay=1e-4)
criterion = nn.MSELoss()

print("\nTraining Autoencoder...")
losses = []
for epoch in range(60):
    ae.train()
    epoch_loss = 0
    for xb, _ in loader:
        opt.zero_grad()
        recon, _ = ae(xb)
        loss = criterion(recon, xb)
        loss.backward()
        opt.step()
        epoch_loss += loss.item()
    losses.append(epoch_loss / len(loader))
    if (epoch + 1) % 10 == 0:
        print(f"  Epoch {epoch+1:3d}/60 | Loss: {losses[-1]:.4f}")

# Extract latent embeddings
ae.eval()
with torch.no_grad():
    _, Z = ae(X_tensor)
Z_np = Z.numpy()

# ── 4. Random Forest Classifier + Cross-Validation ────────────
print("\nTraining Random Forest on concatenated omics features...")
rf_pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("selector", SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=SEED), threshold="mean")),
    ("clf", RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=-1))
])

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
auc_scores = cross_val_score(rf_pipe, X_concat, labels, cv=cv, scoring="roc_auc")
print(f"  5-fold CV AUC: {auc_scores.mean():.3f} ± {auc_scores.std():.3f}")

# Full-data fit for feature importance & SHAP
rf_pipe.fit(X_concat, labels)
rf_model      = rf_pipe.named_steps["clf"]
selector      = rf_pipe.named_steps["selector"]
X_selected    = rf_pipe[:-1].transform(X_concat)
selected_feats = np.array(feature_names)[selector.get_support()]

# ── 5. SHAP Explainability ────────────────────────────────────
print("Computing SHAP values (this may take ~30 seconds)...")
explainer   = shap.TreeExplainer(rf_model)
shap_values = explainer.shap_values(X_selected)
# shap_values[1] = class 1 (case)
sv_class1 = shap_values[1] if isinstance(shap_values, list) else shap_values

fig, ax = plt.subplots(figsize=(10, 7))
shap.summary_plot(sv_class1, X_selected,
                  feature_names=selected_feats,
                  show=False, max_display=20)
plt.title("SHAP Summary — Top 20 Features (Case vs Control)")
plt.tight_layout()
plt.savefig("plots/shap_summary.png", dpi=150)
plt.close()

# ── 6. ROC Curve (train/test split for display) ───────────────
from sklearn.model_selection import train_test_split
X_tr, X_te, y_tr, y_te = train_test_split(X_concat, labels, test_size=0.25, stratify=labels, random_state=SEED)
rf_pipe.fit(X_tr, y_tr)
probs = rf_pipe.predict_proba(X_te)[:, 1]
fpr, tpr, _ = roc_curve(y_te, probs)
auc = roc_auc_score(y_te, probs)

fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(fpr, tpr, color="#e63946", lw=2, label=f"RF (AUC = {auc:.3f})")
ax.plot([0, 1], [0, 1], color="grey", lw=1, ls="--")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ROC Curve — Multi-Omics Classifier")
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig("plots/roc_curve.png", dpi=150)
plt.close()

# ── 7. Latent Space Visualisation (Autoencoder embeddings) ────
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA as skPCA

pca2   = skPCA(n_components=2)
Z_pca  = pca2.fit_transform(Z_np)

tsne   = TSNE(n_components=2, perplexity=30, random_state=SEED)
Z_tsne = tsne.fit_transform(Z_np)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
cols = ["#457b9d" if l == 0 else "#e63946" for l in labels]
for ax, coords, title in zip(axes,
                              [Z_pca, Z_tsne],
                              ["PCA of Autoencoder Latent Space",
                               "t-SNE of Autoencoder Latent Space"]):
    ax.scatter(coords[:, 0], coords[:, 1], c=cols, s=30, alpha=0.8)
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")
import matplotlib.patches as mpatches
p1 = mpatches.Patch(color="#457b9d", label="Control")
p2 = mpatches.Patch(color="#e63946", label="Case")
axes[1].legend(handles=[p1, p2])
plt.suptitle("Autoencoder Latent Representations", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("plots/latent_space.png", dpi=150)
plt.close()

# ── 8. Autoencoder Loss Curve ─────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(losses, color="#2a9d8f", lw=2)
ax.set_xlabel("Epoch")
ax.set_ylabel("Reconstruction Loss (MSE)")
ax.set_title("Autoencoder Training Loss")
plt.tight_layout()
plt.savefig("plots/autoencoder_loss.png", dpi=150)
plt.close()

# ── 9. Save Results ───────────────────────────────────────────
pd.DataFrame({"CV_Fold": range(1, 6), "AUC": auc_scores}).to_csv(
    "results/cross_validation_AUC.csv", index=False)
pd.DataFrame(Z_np, index=sample_ids,
             columns=[f"Latent_{i}" for i in range(Z_np.shape[1])]).to_csv(
    "results/latent_embeddings.csv")
pd.DataFrame({"Feature": selected_feats,
              "MeanSHAP": np.abs(sv_class1).mean(axis=0)
             }).sort_values("MeanSHAP", ascending=False).to_csv(
    "results/feature_importance_SHAP.csv", index=False)

print("\n✅ Multi-omics ML pipeline complete.")
print(f"   5-fold CV AUC: {auc_scores.mean():.3f} | Hold-out AUC: {auc:.3f}")
print("   Outputs saved to results/ and plots/")
