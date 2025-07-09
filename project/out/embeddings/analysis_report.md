
# 🔬 Carbon Cluster Embedding Analysis Report

## 📊 Dataset Overview
- **Total Clusters:** 5,232
- **Features:** 9 columns
- **Numeric Features:** 8

## 📈 Statistical Summary
| Feature | Mean | Std | Min | Max |
|---------|------|-----|-----|-----|
| n_atoms | 239.966 | 2.489 | 60.000 | 240.000 |
| E_cluster | -9.045 | 0.820 | -15.714 | -2.157 |
| BE_per_atom | -7.332 | 0.003 | -7.340 | -7.305 |
| radius_of_gyration | 23.239 | 0.522 | 11.931 | 25.201 |
| sphericity | 0.498 | 0.114 | 0.300 | 0.700 |

## 🧠 Embedding Methods Generated
- **RAW_FEATURES:** 7 dimensions
- **PCA:** 7 dimensions
- **TSNE:** 3 dimensions
- **SVD:** 7 dimensions
- **UMAP:** 3 dimensions

## 💡 Recommendations
- Large dataset: UMAP recommended for preserving global structure
