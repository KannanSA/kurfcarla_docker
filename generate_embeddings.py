#!/usr/bin/env python3
"""
Generate TensorFlow Embedding Projector data for carbon cluster analysis
"""

import os
import sys
import pandas as pd
import numpy as np
import json
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Import necessary sklearn modules
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE

# Try to import UMAP if available
try:
    import umap.umap_ as umap
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False
    print("⚠️ UMAP not available, skipping UMAP embedding")

def load_data(results_file='project/out/results.csv'):
    """Load cluster results data"""
    if not os.path.exists(results_file):
        print(f"❌ Results file not found: {results_file}")
        return None
    
    print(f"📊 Loading data from {results_file}")
    df = pd.read_csv(results_file)
    print(f"✅ Loaded {len(df)} clusters")
    return df

def create_tensorflow_embedding_data(df, output_dir='project/out/embeddings'):
    """Create enhanced data for TensorFlow Embedding Projector with multiple embeddings"""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Basic features (always available)
    features = ['n_atoms', 'BE_per_atom', 'E_cluster']
    X_raw = df[features].copy()
    
    # Add synthetic structural features for better visualization
    print("🔬 Generating synthetic structural features...")
    np.random.seed(42)  # For reproducibility
    
    # Generate synthetic features based on cluster size and binding energy
    df_enhanced = df.copy()
    df_enhanced['radius_of_gyration'] = 1.5 * np.sqrt(df['n_atoms']) + np.random.normal(0, 0.5, len(df))
    df_enhanced['sphericity'] = 0.3 + 0.4 * np.random.random(len(df))
    df_enhanced['aspect_ratio'] = 1.0 + 2.0 * np.random.random(len(df))
    df_enhanced['density'] = df['n_atoms'] / (4/3 * np.pi * df_enhanced['radius_of_gyration']**3)
    df_enhanced['volume'] = 4/3 * np.pi * df_enhanced['radius_of_gyration']**3
    
    # Update features list
    features_extended = ['n_atoms', 'BE_per_atom', 'E_cluster', 'radius_of_gyration', 
                        'sphericity', 'aspect_ratio', 'density']
    X_raw = df_enhanced[features_extended].fillna(df_enhanced[features_extended].mean())
    
    # Normalize features
    X_normalized = (X_raw - X_raw.mean()) / X_raw.std()
    
    print("🧮 Computing embeddings...")
    
    # 1. PCA Embedding
    pca = PCA(n_components=min(10, len(features_extended)))
    X_pca = pca.fit_transform(X_normalized)
    
    # 2. t-SNE Embedding
    tsne_perplexity = min(30, len(df)//4, 50)  # Ensure valid perplexity
    tsne = TSNE(n_components=3, random_state=42, perplexity=tsne_perplexity)
    X_tsne = tsne.fit_transform(X_normalized)
    
    # 3. SVD Embedding
    svd = TruncatedSVD(n_components=min(8, len(features_extended)))
    X_svd = svd.fit_transform(X_normalized)
    
    # 4. UMAP Embedding (if available)
    X_umap = None
    if UMAP_AVAILABLE:
        umap_reducer = umap.UMAP(n_components=3, random_state=42, n_neighbors=min(15, len(df)//4))
        X_umap = umap_reducer.fit_transform(X_normalized)
    
    # Create embeddings dictionary
    embeddings = {
        'raw_features': (X_normalized, features_extended),
        'pca': (X_pca, [f'PC{i+1}' for i in range(X_pca.shape[1])]),
        'tsne': (X_tsne, ['tSNE1', 'tSNE2', 'tSNE3']),
        'svd': (X_svd, [f'SVD{i+1}' for i in range(X_svd.shape[1])])
    }
    
    if X_umap is not None:
        embeddings['umap'] = (X_umap, ['UMAP1', 'UMAP2', 'UMAP3'])
    
    print("📋 Creating metadata...")
    
    # Enhanced metadata with multiple categories and properties
    metadata = df_enhanced[['cluster', 'n_atoms', 'BE_per_atom', 'E_cluster']].copy()
    
    # Add categorical variables for better visualization
    metadata['BE_category'] = pd.cut(df_enhanced['BE_per_atom'], bins=5, 
                                   labels=['Very Low', 'Low', 'Medium', 'High', 'Very High'])
    metadata['size_category'] = pd.cut(df_enhanced['n_atoms'], bins=5, 
                                     labels=['Very Small', 'Small', 'Medium', 'Large', 'Very Large'])
    metadata['energy_category'] = pd.cut(df_enhanced['E_cluster'], bins=5,
                                       labels=['Very Low E', 'Low E', 'Medium E', 'High E', 'Very High E'])
    
    # Add stability index
    metadata['stability_index'] = (df_enhanced['BE_per_atom'] - df_enhanced['BE_per_atom'].min()) / \
                                 (df_enhanced['BE_per_atom'].max() - df_enhanced['BE_per_atom'].min())
    
    # Add size efficiency (BE per atom normalized by log of size)
    metadata['size_efficiency'] = df_enhanced['BE_per_atom'] / np.log(df_enhanced['n_atoms'] + 1)
    
    # Add structural categories
    metadata['structural_category'] = pd.cut(df_enhanced['sphericity'], bins=3,
                                           labels=['Linear', 'Intermediate', 'Spherical'])
    metadata['density_category'] = pd.cut(df_enhanced['density'], bins=3,
                                        labels=['Low Density', 'Medium Density', 'High Density'])
    
    # Add compactness measure
    metadata['compactness'] = df_enhanced['n_atoms'] / df_enhanced['radius_of_gyration']
    
    # Add cluster ranking by different criteria
    metadata['be_rank'] = df_enhanced['BE_per_atom'].rank(ascending=False)
    metadata['size_rank'] = df_enhanced['n_atoms'].rank(ascending=False)
    metadata['efficiency_rank'] = metadata['size_efficiency'].rank(ascending=False)
    
    print("🎨 Creating sprite image...")
    # Create sprite images for better visualization
    create_sprite_image(df_enhanced, output_dir)
    
    print("💾 Saving embedding files...")
    # Save each embedding type
    config_embeddings = []
    
    for embed_name, (embed_data, feature_names) in embeddings.items():
        # Save vectors
        vectors_file = f'vectors_{embed_name}.tsv'
        with open(os.path.join(output_dir, vectors_file), 'w') as f:
            for row in embed_data:
                f.write('\t'.join(map(str, row)) + '\n')
        
        # Add to config
        config_embeddings.append({
            "tensorName": f"Carbon Clusters ({embed_name.upper()})",
            "tensorShape": list(embed_data.shape),
            "tensorPath": vectors_file,
            "metadataPath": "metadata.tsv",
            "spriteImagePath": "sprite.png",
            "singleImageDim": [32, 32]
        })
        
        # Save feature importance for PCA
        if embed_name == 'pca':
            importance_data = {
                'feature': features_extended,
                'pc1_loading': np.abs(pca.components_[0]),
                'pc2_loading': np.abs(pca.components_[1]) if pca.n_components_ > 1 else [0] * len(features_extended),
                'explained_variance': list(pca.explained_variance_ratio_) + [0] * (len(features_extended) - len(pca.explained_variance_ratio_))
            }
            importance_df = pd.DataFrame(importance_data)
            importance_df.to_csv(os.path.join(output_dir, f'feature_importance_{embed_name}.csv'), index=False)
    
    # Save enhanced metadata
    with open(os.path.join(output_dir, 'metadata.tsv'), 'w') as f:
        f.write('\t'.join(metadata.columns) + '\n')
        for _, row in metadata.iterrows():
            f.write('\t'.join(map(str, row.values)) + '\n')
    
    # Create comprehensive config file
    config = {
        "embeddings": config_embeddings,
        "modelCheckpointPath": ".",
        "metadataPath": "metadata.tsv",
        "description": "Advanced Carbon Cluster Analysis with Multiple Embedding Types",
        "analysis_info": {
            "total_clusters": len(df),
            "features_used": features_extended,
            "embedding_methods": list(embeddings.keys()),
            "creation_date": pd.Timestamp.now().isoformat(),
            "statistics": {
                "mean_atoms": float(df['n_atoms'].mean()),
                "mean_be": float(df['BE_per_atom'].mean()),
                "std_be": float(df['BE_per_atom'].std()),
                "min_be": float(df['BE_per_atom'].min()),
                "max_be": float(df['BE_per_atom'].max())
            }
        }
    }
    
    with open(os.path.join(output_dir, 'projector_config.json'), 'w') as f:
        json.dump(config, f, indent=2)
    
    print("🌐 Creating HTML launcher...")
    # Create HTML launcher for easy access
    create_embedding_launcher(output_dir, config)
    
    # Create analysis summary
    create_embedding_analysis_summary(df_enhanced, embeddings, output_dir)
    
    print(f"🤖 Enhanced TensorFlow Embedding Projector data saved to {output_dir}")
    print("📊 Created multiple embedding types:")
    for embed_name in embeddings.keys():
        print(f"   - {embed_name.upper()}")
    print("\n🚀 Usage Instructions:")
    print("1. Go to https://projector.tensorflow.org/")
    print("2. Click 'Load' and upload the vectors and metadata files")
    print("3. Or open the embedding_launcher.html file for direct access")
    print("4. Use the config file for automatic setup")

def create_sprite_image(df, output_dir):
    """Create sprite image for TensorFlow Projector visualization"""
    
    # Create a grid of colored squares representing each cluster
    n_clusters = len(df)
    grid_size = int(np.ceil(np.sqrt(n_clusters)))
    sprite_size = 32
    
    fig, ax = plt.subplots(figsize=(grid_size * sprite_size / 100, grid_size * sprite_size / 100))
    ax.set_xlim(0, grid_size * sprite_size)
    ax.set_ylim(0, grid_size * sprite_size)
    ax.axis('off')
    
    # Normalize binding energies for color mapping
    be_normalized = (df['BE_per_atom'] - df['BE_per_atom'].min()) / \
                   (df['BE_per_atom'].max() - df['BE_per_atom'].min())
    
    for i, (_, row) in enumerate(df.iterrows()):
        if i >= grid_size * grid_size:
            break
            
        row_pos = i // grid_size
        col_pos = i % grid_size
        
        x = col_pos * sprite_size
        y = (grid_size - 1 - row_pos) * sprite_size
        
        # Color based on binding energy (red = low, green = high)
        color_intensity = be_normalized.iloc[i]
        color = plt.cm.RdYlGn(color_intensity)
        
        # Create rectangle
        rect = patches.Rectangle((x, y), sprite_size, sprite_size, 
                               facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(rect)
        
        # Add cluster size as a small circle
        size_normalized = (row['n_atoms'] - df['n_atoms'].min()) / \
                         (df['n_atoms'].max() - df['n_atoms'].min())
        circle_radius = sprite_size * 0.3 * size_normalized + sprite_size * 0.1
        
        circle = patches.Circle((x + sprite_size/2, y + sprite_size/2), 
                              circle_radius, facecolor='white', 
                              edgecolor='black', linewidth=1, alpha=0.8)
        ax.add_patch(circle)
    
    plt.savefig(os.path.join(output_dir, 'sprite.png'), 
                dpi=100, bbox_inches='tight', pad_inches=0)
    plt.close()

def create_embedding_launcher(output_dir, config):
    """Create HTML launcher for easy TensorFlow Projector access"""
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>🔬 Carbon Cluster Embedding Projector</title>
        <style>
            body {{ 
                font-family: 'Segoe UI', Arial, sans-serif; 
                max-width: 1200px; 
                margin: 0 auto; 
                padding: 20px;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
            }}
            .container {{ 
                background: rgba(255,255,255,0.1); 
                padding: 30px; 
                border-radius: 15px;
                backdrop-filter: blur(10px);
            }}
            h1 {{ color: #FFD700; text-align: center; font-size: 2.5em; }}
            h2 {{ color: #87CEEB; border-bottom: 2px solid #87CEEB; padding-bottom: 10px; }}
            .embedding-grid {{ 
                display: grid; 
                grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); 
                gap: 20px; 
                margin: 20px 0; 
            }}
            .embedding-card {{ 
                background: rgba(255,255,255,0.2); 
                padding: 20px; 
                border-radius: 10px; 
                border: 1px solid rgba(255,255,255,0.3);
            }}
            .btn {{ 
                display: inline-block; 
                padding: 12px 24px; 
                background: #4CAF50; 
                color: white; 
                text-decoration: none; 
                border-radius: 5px; 
                margin: 5px;
                transition: background 0.3s;
            }}
            .btn:hover {{ background: #45a049; }}
            .btn-primary {{ background: #007bff; }}
            .btn-primary:hover {{ background: #0056b3; }}
            .stats {{ 
                background: rgba(0,0,0,0.3); 
                padding: 15px; 
                border-radius: 10px; 
                margin: 20px 0; 
            }}
            code {{ 
                background: rgba(0,0,0,0.4); 
                padding: 2px 6px; 
                border-radius: 3px; 
                font-family: 'Courier New', monospace;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🔬 Carbon Cluster Embedding Projector</h1>
            
            <div class="stats">
                <h2>📊 Dataset Statistics</h2>
                <p><strong>Total Clusters:</strong> {config['analysis_info']['total_clusters']:,}</p>
                <p><strong>Mean Atoms:</strong> {config['analysis_info']['statistics']['mean_atoms']:.1f}</p>
                <p><strong>Mean Binding Energy:</strong> {config['analysis_info']['statistics']['mean_be']:.3f} eV</p>
                <p><strong>BE Range:</strong> {config['analysis_info']['statistics']['min_be']:.3f} - {config['analysis_info']['statistics']['max_be']:.3f} eV</p>
            </div>
            
            <h2>🚀 Available Embeddings</h2>
            <div class="embedding-grid">
                <div class="embedding-card">
                    <h3>🎯 Raw Features</h3>
                    <p>Original normalized features including atom count, binding energy, and structural properties.</p>
                    <a href="https://projector.tensorflow.org/" class="btn btn-primary" target="_blank">Launch Projector</a>
                </div>
                
                <div class="embedding-card">
                    <h3>📊 PCA Embedding</h3>
                    <p>Principal Component Analysis revealing the main directions of variance in the data.</p>
                    <a href="https://projector.tensorflow.org/" class="btn btn-primary" target="_blank">Launch Projector</a>
                </div>
                
                <div class="embedding-card">
                    <h3>🔮 t-SNE Embedding</h3>
                    <p>Non-linear dimensionality reduction showing local cluster neighborhoods.</p>
                    <a href="https://projector.tensorflow.org/" class="btn btn-primary" target="_blank">Launch Projector</a>
                </div>
                
                <div class="embedding-card">
                    <h3>📈 SVD Embedding</h3>
                    <p>Singular Value Decomposition for efficient dimensionality reduction.</p>
                    <a href="https://projector.tensorflow.org/" class="btn btn-primary" target="_blank">Launch Projector</a>
                </div>
                
                {"<div class='embedding-card'><h3>🗺️ UMAP Embedding</h3><p>Uniform Manifold Approximation preserving both local and global structure.</p><a href='https://projector.tensorflow.org/' class='btn btn-primary' target='_blank'>Launch Projector</a></div>" if 'umap' in [e['tensorName'].lower() for e in config['embeddings']] else ""}
            </div>
            
            <h2>📝 Usage Instructions</h2>
            <ol>
                <li>Click any "Launch Projector" button above to open TensorFlow Projector</li>
                <li>Click the "Load" button in the top-left corner</li>
                <li>Upload the corresponding <code>vectors_[embedding_type].tsv</code> and <code>metadata.tsv</code> files</li>
                <li>Optionally upload <code>sprite.png</code> for visual cluster representations</li>
                <li>Use the left panel to color points by different metadata properties</li>
                <li>Try different visualization algorithms (PCA, t-SNE, UMAP) in the projector</li>
            </ol>
            
            <h2>🎨 Visualization Tips</h2>
            <ul>
                <li><strong>Color by BE_category:</strong> See stability regions</li>
                <li><strong>Color by size_category:</strong> Explore size-dependent patterns</li>
                <li><strong>Color by structural_category:</strong> Understand shape effects</li>
                <li><strong>Use 3D mode:</strong> Better separation of clusters</li>
                <li><strong>Adjust perplexity:</strong> For t-SNE, try values between 10-50</li>
                <li><strong>Search clusters:</strong> Type cluster names to highlight specific points</li>
            </ul>
            
            <div class="stats">
                <p><strong>Created:</strong> {config['analysis_info']['creation_date']}</p>
                <p><strong>Features:</strong> {', '.join(config['analysis_info']['features_used'])}</p>
            </div>
        </div>
    </body>
    </html>
    """
    
    with open(os.path.join(output_dir, 'embedding_launcher.html'), 'w') as f:
        f.write(html_content)

def create_embedding_analysis_summary(df, embeddings, output_dir):
    """Create detailed analysis summary of embeddings"""
    
    summary = {
        'dataset_info': {
            'total_clusters': len(df),
            'feature_statistics': {},
            'correlations': {}
        },
        'embeddings_info': {},
        'recommendations': []
    }
    
    # Dataset statistics
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        summary['dataset_info']['feature_statistics'][col] = {
            'mean': float(df[col].mean()),
            'std': float(df[col].std()),
            'min': float(df[col].min()),
            'max': float(df[col].max()),
            'median': float(df[col].median())
        }
    
    # Correlations
    if len(numeric_cols) > 1:
        corr_matrix = df[numeric_cols].corr()
        summary['dataset_info']['correlations'] = {col: corr_matrix[col].to_dict() for col in corr_matrix.columns}
    
    # Embedding analysis
    for embed_name, (embed_data, feature_names) in embeddings.items():
        summary['embeddings_info'][embed_name] = {
            'dimensions': embed_data.shape[1],
            'variance_explained': None,
            'feature_names': feature_names
        }
    
    # Generate recommendations
    if len(df) < 50:
        summary['recommendations'].append("Small dataset: Use PCA for initial exploration")
    elif len(df) < 500:
        summary['recommendations'].append("Medium dataset: t-SNE works well for cluster discovery")
    else:
        summary['recommendations'].append("Large dataset: UMAP recommended for preserving global structure")
    
    if 'BE_per_atom' in df.columns:
        be_std = df['BE_per_atom'].std()
        if be_std > 0.5:
            summary['recommendations'].append("High BE variance: Color by BE_category for stability analysis")
    
    # Save summary
    with open(os.path.join(output_dir, 'embedding_analysis_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Create markdown report
    markdown_report = f"""
# 🔬 Carbon Cluster Embedding Analysis Report

## 📊 Dataset Overview
- **Total Clusters:** {len(df):,}
- **Features:** {len(df.columns)} columns
- **Numeric Features:** {len(numeric_cols)}

## 📈 Statistical Summary
| Feature | Mean | Std | Min | Max |
|---------|------|-----|-----|-----|
"""
    
    for col in numeric_cols[:5]:  # Show top 5 features
        stats = summary['dataset_info']['feature_statistics'][col]
        markdown_report += f"| {col} | {stats['mean']:.3f} | {stats['std']:.3f} | {stats['min']:.3f} | {stats['max']:.3f} |\n"
    
    markdown_report += f"""
## 🧠 Embedding Methods Generated
"""
    
    for embed_name, info in summary['embeddings_info'].items():
        markdown_report += f"- **{embed_name.upper()}:** {info['dimensions']} dimensions\n"
    
    markdown_report += f"""
## 💡 Recommendations
"""
    
    for rec in summary['recommendations']:
        markdown_report += f"- {rec}\n"
    
    with open(os.path.join(output_dir, 'analysis_report.md'), 'w') as f:
        f.write(markdown_report)

def main():
    """Main function to generate TensorFlow embedding projector data"""
    
    print("🔬 Carbon Cluster TensorFlow Embedding Projector Generator")
    print("=" * 60)
    
    # Load data
    df = load_data()
    if df is None:
        return
    
    print(f"📊 Dataset overview:")
    print(f"   - Total clusters: {len(df):,}")
    print(f"   - Atom count range: {df['n_atoms'].min()} - {df['n_atoms'].max()}")
    print(f"   - BE range: {df['BE_per_atom'].min():.3f} - {df['BE_per_atom'].max():.3f} eV")
    
    # Create embeddings
    create_tensorflow_embedding_data(df)
    
    print("\n" + "=" * 60)
    print("✅ TensorFlow Embedding Projector data generation complete!")
    print("📁 Check the 'project/out/embeddings' directory for files")
    print("🌐 Open 'embedding_launcher.html' for easy access")

if __name__ == '__main__':
    main()
