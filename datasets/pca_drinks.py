import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler

df = pd.read_csv("synthetic_drinks.csv")

avatar_cols = ["sweetness", "fizziness", "overall_rating"]   # mapped to avatar

spatial_cols = []
pca_cols = ["price", "caffeine", 'sourness', 'sugar', 'citric_acid', 'co2_pressure', 'preservatives', 'color_intensity']
short_labels = [ "PRI", "CAF", "SOU", "SUG", "CIT", "CO2", "PRE", "COL" ]

# spatial_cols = ["price", "caffeine"]
# pca_cols = ['sourness', 'sugar', 'citric_acid', 'co2_pressure', 'preservatives', 'color_intensity']
# short_labels = [ "SOU", "SUG", "CIT", "CO2", "PRE", "COL" ]

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

    # Heatmap
    sns.heatmap(loadings, annot=True, fmt=".2f", center=0, cmap='vlag',
            yticklabels=short_labels, xticklabels=['PC1', 'PC2'], ax=axes[0])
    
    axes[0].set_title("PCA loadings")

    # Scatter plot PCA
    sc = axes[1].scatter(pcs[:,0], pcs[:,1], c=df["overall_rating"], cmap='viridis', s=70)
    for i, n in enumerate(df["name"]):
        axes[1].text(pcs[i,0]+0.05, pcs[i,1], n, fontsize=8)
    axes[1].set_xlabel("PC1")
    axes[1].set_ylabel("PC2")
    axes[1].set_title("PCA scatter (colored by overall_rating)")
    cb = plt.colorbar(sc, ax=axes[1], fraction=0.045)
    cb.set_label("Overall_Rating")

    plt.tight_layout()
    plt.savefig("politicians_pca_combo.png", dpi=300)
    plt.show()


else:
    print("No anonymous features selected for PCA.")

print("Sample avatar feature values:\n", df[avatar_cols].head())
print("Sample spatial feature values:\n", df[spatial_cols].head())

