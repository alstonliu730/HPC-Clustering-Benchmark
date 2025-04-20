
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import numpy as np

# Load the CSV file
file_path = 'output.csv'  # Update this path as needed
df = pd.read_csv(file_path, header=None, names=["x", "y", "z", "cluster_id"])

# Get unique cluster IDs
unique_clusters = np.unique(df['cluster_id'])
unique_clusters.sort()

# Use a large colormap to support many clusters
cmap = cm.get_cmap('nipy_spectral', len(unique_clusters))
color_dict = {}

for i, cid in enumerate(unique_clusters):
    if cid == 0:
        color_dict[cid] = (0.5, 0.5, 0.5, 1.0)  # gray for noise
    else:
        color_dict[cid] = cmap(i % cmap.N)  # wrap colors if needed

# Map cluster IDs to colors
colors = df['cluster_id'].map(color_dict)

# Plotting
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(df['x'], df['y'], df['z'], c=colors, s=1, marker='o')

ax.set_title("3D Cluster Visualization")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.tight_layout()
plt.show()
