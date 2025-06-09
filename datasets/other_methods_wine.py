import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

df = pd.read_csv('VinhoVerde/winequality-red.csv', delimiter=";")

feature_cols = [
    "fixed_acidity", "volatile_acidity", "citric_acid", "residual_sugar",
    "chlorides", "free_sulfur_dioxide", "total_sulfur_dioxide",
    "density", "pH", "sulphates", "alcohol"
]
label_col = 'quality'

if label_col:
    uniq_q = sorted(df[label_col].unique())
    color_map = plt.colormaps['plasma'].resampled(len(uniq_q))

    df['color'] = df[label_col].map(lambda q: color_map(uniq_q.index(q)))
else:
    df['color'] = 'gray'

fig, axs = plt.subplots(1, 2, figsize=(14,6))

# PCA scatter
X = df[feature_cols].values
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)
axs[0].scatter(
    X_pca[:,0], X_pca[:,1],
    c=df['color'].tolist(),
    s=18, alpha=0.6, edgecolor='none'
)
axs[0].set_xlabel('PC1')
axs[0].set_ylabel('PC2')
axs[0].set_title('PCA scatter (colored by quality)')

# t-SNE scatter
tsne = TSNE(n_components=2, random_state=42, perplexity=40)
X_tsne = tsne.fit_transform(X)
axs[1].scatter(
    X_tsne[:,0], X_tsne[:,1],
    c=df['color'].tolist(),
    s=18, alpha=0.6, edgecolor='none'
)
axs[1].set_title('t-SNE projection (colored by quality)')

plt.tight_layout()
fig.suptitle('Classical projection methods for the wine dataset', fontsize=15)
plt.subplots_adjust(top=0.89)
plt.show()
