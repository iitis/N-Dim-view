import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler

df = pd.read_csv("winequality-red.csv", delimiter=";")

avatar_cols = ['fixed_acidity', 'total_sulfur_dioxide', 'pH', 'alcohol', 'quality' ]
spatial_cols = []
pca_cols = ['residual_sugar', 'volatile_acidity', 'citric_acid', 'chlorides', 'free_sulfur_dioxide', 'density', 'sulphates']
short_labels = ["SUG", "VAC", "CAC", "CLO", "FSD", "DEN", "SUL"]


print("Avatar features:", avatar_cols)
print("Spatial features:", spatial_cols)
print("PCA (anonymous) features:", pca_cols)

if pca_cols:
    X_pca = df[pca_cols].values
    scaler = StandardScaler()
    X_pca_scaled = scaler.fit_transform(X_pca)

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X_pca_scaled)
    loadings = pca.components_.T
    
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={'width_ratios':[1,1.2]})

    sns.heatmap(loadings, annot=True, fmt=".2f", center=0, cmap='vlag',
            yticklabels=short_labels, xticklabels=['PC1', 'PC2'], ax=axes[0])
    
    axes[0].set_title("PCA loadings")

    sc = axes[1].scatter(pcs[:,0], pcs[:,1], c=df["quality"], cmap='viridis', s=70)
    axes[1].set_xlabel("PC1")
    axes[1].set_ylabel("PC2")
    axes[1].set_title("PCA scatter (colored by quality)")
    cb = plt.colorbar(sc, ax=axes[1], fraction=0.045)
    cb.set_label("Quality")

    plt.tight_layout()
    plt.savefig("politicians_pca_combo.png", dpi=300)
    plt.show()


else:
    print("No anonymous features selected for PCA.")

print("Sample avatar feature values:\n", df[avatar_cols].head())
print("Sample spatial feature values:\n", df[spatial_cols].head())
