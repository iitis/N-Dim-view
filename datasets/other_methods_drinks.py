import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from pandas.plotting import parallel_coordinates
import matplotlib.patches as patches

df = pd.read_csv('synthetic_drinks.csv')

feature_cols = [
    'sweetness', 'sourness', 'fizziness', 'color_intensity', 'overall_rating',
    'caffeine', 'preservatives', 'price'
]
label_col = 'group_numeric'
group_names = {0: "Classic", 1: "Premium", 2: "Extreme"}
group_colors = ['tab:blue', 'tab:orange', 'tab:green']

fig, axs = plt.subplots(2, 2, figsize=(16,12))

# PCA scatter
X = df[feature_cols].values
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)
for grp in sorted(df[label_col].unique()):
    axs[0,0].scatter(
        X_pca[df[label_col]==grp,0],
        X_pca[df[label_col]==grp,1],
        label=group_names.get(grp, f"Group {grp}"),
        alpha=0.7, s=36, color=group_colors[grp % len(group_colors)]
    )
axs[0,0].set_xlabel('PC1')
axs[0,0].set_ylabel('PC2')
axs[0,0].set_title('PCA scatter')
axs[0,0].legend(fontsize=9)

# t-SNE scatter
tsne = TSNE(n_components=2, random_state=42, perplexity=20)
X_tsne = tsne.fit_transform(X)
for grp in sorted(df[label_col].unique()):
    axs[0,1].scatter(
        X_tsne[df[label_col]==grp,0],
        X_tsne[df[label_col]==grp,1],
        label=group_names.get(grp, f"Group {grp}"),
        alpha=0.7, s=36, color=group_colors[grp % len(group_colors)]
    )
axs[0,1].set_title('t-SNE projection')

# Parallel coordinates
pc_df = df.copy()
pc_df[label_col] = pc_df[label_col].astype(str)
parallel_coordinates(pc_df[[label_col]+feature_cols], class_column=label_col, alpha=0.18, ax=axs[1,0], color=group_colors)
axs[1,0].set_title('Parallel coordinates')
axs[1,0].tick_params(axis='x', rotation=30, labelsize=8)

# Chernoff faces
def plot_face(ax, features, cmap=plt.cm.viridis):
    features = np.clip(features, 0, 1)
    face = patches.Ellipse((0,0), 1, 1.2, color=cmap(features[0]), alpha=0.6)
    ax.add_patch(face)
    ax.add_patch(patches.Ellipse((-0.25,0.2), 0.14,0.10, color='black'))
    ax.add_patch(patches.Ellipse((0.25,0.2), 0.14,0.10, color='black'))
    ax.add_patch(patches.Ellipse((0,0), 0.09,0.15*features[1]+0.07, color='brown'))
    mouth_x = np.linspace(-0.2,0.2,100)
    mouth_y = -0.25 + 0.15*(features[2]-0.5)*np.sin(np.pi*mouth_x/0.4)
    ax.plot(mouth_x, mouth_y, color='red', lw=2)
    face.set_edgecolor(cmap(features[3]))
    face.set_linewidth(2+2*features[4])
    ax.set_xlim(-0.6,0.6)
    ax.set_ylim(-0.7,0.7)
    ax.axis('off')

glyph_cols = ['sweetness','sourness','fizziness','color_intensity','overall_rating']
glyph_arr = df[glyph_cols].to_numpy()
glyph_min = np.min(glyph_arr, axis=0)
glyph_max = np.max(glyph_arr, axis=0)
glyph_norm = (glyph_arr - glyph_min) / (glyph_max - glyph_min + 1e-9)

n_faces = 24
n_rows, n_cols = 4, 6
idxs = np.linspace(0, len(df)-1, n_faces, dtype=int)
axs[1,1].set_title("Chernoff faces (examples, 6x4 grid)")
axs[1,1].axis('off')
dx, dy = 1.0 / n_cols, 1.0 / n_rows
for i, idx in enumerate(idxs):
    r, c = divmod(i, n_cols)
    left = c * dx
    bottom = 1.0 - (r+1) * dy
    rect = [left, bottom, dx*0.92, dy*0.88]
    sub_ax = axs[1,1].inset_axes(rect)
    feat = glyph_norm[idx]
    g = df.iloc[idx][label_col]
    plot_face(sub_ax, feat)
    sub_ax.set_title(f"G{int(g)}", fontsize=8)
    sub_ax.axis('off')

plt.tight_layout()
plt.subplots_adjust(top=0.92)
fig.suptitle('Comparison of classical visualization methods for the synthetic soft drinks dataset', fontsize=16)
plt.show()
