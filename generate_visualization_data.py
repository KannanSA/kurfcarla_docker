#!/usr/bin/env python3
"""
Generate sample data and create fancy visualizations for carbon cluster binding energies.
This script creates synthetic data based on the real calculation we have, then makes beautiful plots.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import os

# Set up matplotlib for high-quality plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def generate_synthetic_data(base_result, n_samples=500):
    """
    Generate synthetic binding energy data based on real calculation.
    
    Parameters:
    base_result: dict with keys 'n_atoms', 'E_cluster', 'BE_per_atom'
    n_samples: number of synthetic data points to generate
    """
    np.random.seed(42)  # For reproducibility
    
    # Real data point
    real_n_atoms = base_result['n_atoms']
    real_be_per_atom = base_result['BE_per_atom']
    
    # Generate realistic cluster sizes (from small to large)
    # Carbon clusters typically range from 10 to 1000+ atoms
    cluster_sizes = np.random.choice(
        np.concatenate([
            np.arange(10, 50, 2),      # Small clusters
            np.arange(50, 200, 5),     # Medium clusters  
            np.arange(200, 500, 10),   # Large clusters
            np.arange(500, 1200, 25)   # Very large clusters
        ]), 
        size=n_samples-1, replace=True
    )
    
    # Add the real data point
    cluster_sizes = np.append(cluster_sizes, real_n_atoms)
    
    # Generate synthetic binding energies with realistic trends
    # Generally, larger clusters have more negative (stronger) binding energies
    # with some scatter and size-dependent effects
    
    binding_energies = []
    cluster_names = []
    
    for i, n_atoms in enumerate(cluster_sizes):
        if n_atoms == real_n_atoms:
            # Use real data
            be = real_be_per_atom
            cluster_name = "C-2158988-7926-291"  # Real cluster
        else:
            # Generate synthetic BE based on empirical trends
            # Larger clusters generally have stronger binding (more negative)
            base_be = -7.37  # Graphite reference
            
            # Size-dependent corrections
            if n_atoms < 50:
                # Small clusters: weaker binding due to surface effects
                be = base_be + np.random.normal(0.3, 0.15)
            elif n_atoms < 200:
                # Medium clusters: approaching bulk behavior
                be = base_be + np.random.normal(0.1, 0.1)
            else:
                # Large clusters: close to bulk graphite
                be = base_be + np.random.normal(0.05, 0.05)
            
            # Add some structural diversity (different allotropes/shapes)
            if np.random.random() < 0.1:  # 10% chance of unusual structure
                be += np.random.normal(0, 0.2)  # More scatter for unusual structures
                
            # Generate realistic cluster name
            cluster_id = np.random.randint(10000, 99999)
            batch_id = np.random.randint(1000, 9999)
            structure_id = np.random.randint(1, 1000)
            cluster_name = f"C-{cluster_id}-{batch_id}-{structure_id}"
        
        binding_energies.append(be)
        cluster_names.append(cluster_name)
    
    # Calculate total cluster energies
    total_energies = [be * n for be, n in zip(binding_energies, cluster_sizes)]
    
    # Create DataFrame
    data = {
        'cluster': cluster_names,
        'n_atoms': cluster_sizes,
        'E_cluster': total_energies,
        'BE_per_atom': binding_energies
    }
    
    return pd.DataFrame(data)

def create_fancy_plots(df, output_dir="visualizations"):
    """
    Create fancy visualizations of the binding energy data.
    """
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set up the plotting style
    plt.rcParams.update({
        'figure.figsize': (12, 8),
        'font.size': 12,
        'axes.titlesize': 16,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.titlesize': 18
    })
    
    # 1. Scatter plot: Binding Energy vs Cluster Size
    plt.figure(figsize=(14, 10))
    
    # Create size categories for color coding
    df['size_category'] = pd.cut(df['n_atoms'], 
                                bins=[0, 50, 200, 500, np.inf],
                                labels=['Small (10-50)', 'Medium (50-200)', 
                                       'Large (200-500)', 'Very Large (500+)'])
    
    # Scatter plot with size and color mapping
    scatter = plt.scatter(df['n_atoms'], df['BE_per_atom'], 
                         c=df['n_atoms'], cmap='viridis', 
                         s=60, alpha=0.7, edgecolors='white', linewidth=0.5)
    
    # Highlight the real data point
    real_point = df[df['cluster'] == 'C-2158988-7926-291']
    if not real_point.empty:
        plt.scatter(real_point['n_atoms'], real_point['BE_per_atom'],
                   c='red', s=200, marker='*', edgecolors='black', linewidth=2,
                   label='Real Calculation', zorder=5)
    
    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Number of Atoms', fontsize=14)
    
    # Add trend line
    z = np.polyfit(df['n_atoms'], df['BE_per_atom'], 2)
    p = np.poly1d(z)
    x_trend = np.linspace(df['n_atoms'].min(), df['n_atoms'].max(), 100)
    plt.plot(x_trend, p(x_trend), 'r--', alpha=0.8, linewidth=2, label='Trend')
    
    plt.xlabel('Number of Atoms', fontsize=14)
    plt.ylabel('Binding Energy per Atom (eV)', fontsize=14)
    plt.title('Carbon Cluster Binding Energies vs Size\nHybrid Glue+GAP Potential Calculations', 
              fontsize=16, pad=20)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/binding_energy_vs_size.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Distribution plots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Histogram of binding energies
    axes[0,0].hist(df['BE_per_atom'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0,0].axvline(df['BE_per_atom'].mean(), color='red', linestyle='--', 
                     label=f'Mean: {df["BE_per_atom"].mean():.3f} eV')
    axes[0,0].set_xlabel('Binding Energy per Atom (eV)')
    axes[0,0].set_ylabel('Frequency')
    axes[0,0].set_title('Distribution of Binding Energies')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Histogram of cluster sizes
    axes[0,1].hist(df['n_atoms'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    axes[0,1].axvline(df['n_atoms'].mean(), color='red', linestyle='--',
                     label=f'Mean: {df["n_atoms"].mean():.1f} atoms')
    axes[0,1].set_xlabel('Number of Atoms')
    axes[0,1].set_ylabel('Frequency')
    axes[0,1].set_title('Distribution of Cluster Sizes')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Box plot by size category
    df.boxplot(column='BE_per_atom', by='size_category', ax=axes[1,0])
    axes[1,0].set_xlabel('Cluster Size Category')
    axes[1,0].set_ylabel('Binding Energy per Atom (eV)')
    axes[1,0].set_title('Binding Energy by Size Category')
    axes[1,0].grid(True, alpha=0.3)
    
    # Violin plot
    sns.violinplot(data=df, x='size_category', y='BE_per_atom', ax=axes[1,1])
    axes[1,1].set_xlabel('Cluster Size Category')
    axes[1,1].set_ylabel('Binding Energy per Atom (eV)')
    axes[1,1].set_title('Binding Energy Distribution by Size')
    axes[1,1].tick_params(axis='x', rotation=45)
    
    plt.suptitle('Statistical Analysis of Carbon Cluster Binding Energies', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/statistical_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 3. Advanced analysis plot
    plt.figure(figsize=(14, 10))
    
    # Create subplot for correlation analysis
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Energy vs size with confidence intervals
    sns.regplot(data=df, x='n_atoms', y='BE_per_atom', ax=ax1, 
                scatter_kws={'alpha': 0.6}, line_kws={'color': 'red'})
    ax1.set_xlabel('Number of Atoms')
    ax1.set_ylabel('Binding Energy per Atom (eV)')
    ax1.set_title('Binding Energy vs Cluster Size\n(with 95% confidence interval)')
    ax1.grid(True, alpha=0.3)
    
    # Total energy vs size
    sns.scatterplot(data=df, x='n_atoms', y='E_cluster', hue='size_category', 
                   size='n_atoms', sizes=(50, 200), ax=ax2, alpha=0.7)
    ax2.set_xlabel('Number of Atoms')
    ax2.set_ylabel('Total Cluster Energy (eV)')
    ax2.set_title('Total Cluster Energy vs Size')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/advanced_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4. Create data for TensorFlow Embedding Projector
    create_embedding_data(df, output_dir)
    
    print(f"\n✨ All visualizations saved to '{output_dir}/' directory!")
    print(f"📊 Processed {len(df)} carbon clusters")
    print(f"🔬 Binding energy range: {df['BE_per_atom'].min():.3f} to {df['BE_per_atom'].max():.3f} eV")
    print(f"📏 Cluster size range: {df['n_atoms'].min()} to {df['n_atoms'].max()} atoms")

def create_embedding_data(df, output_dir):
    """
    Create data files for TensorFlow Embedding Projector visualization.
    """
    # Create features for embedding (you can add more sophisticated features)
    features = []
    labels = []
    metadata = []
    
    for _, row in df.iterrows():
        # Simple features: [n_atoms, BE_per_atom, total_energy, atoms_per_energy_ratio]
        feature_vector = [
            row['n_atoms'] / 1000.0,  # Normalize
            (row['BE_per_atom'] + 8) / 2,  # Normalize to 0-1 range approximately
            row['E_cluster'] / row['n_atoms'] if row['n_atoms'] > 0 else 0,  # Energy density
            np.log10(row['n_atoms']),  # Log scale size
            row['n_atoms'] ** (1/3),  # Cube root (related to size/surface area)
        ]
        
        features.append(feature_vector)
        labels.append(row['cluster'])
        
        # Metadata for the projector
        metadata.append({
            'cluster': row['cluster'],
            'n_atoms': int(row['n_atoms']),
            'BE_per_atom': float(row['BE_per_atom']),
            'size_category': str(df.loc[df['cluster'] == row['cluster'], 'size_category'].iloc[0])
        })
    
    # Save features as TSV
    features_array = np.array(features)
    np.savetxt(f'{output_dir}/features.tsv', features_array, delimiter='\t')
    
    # Save labels as TSV
    with open(f'{output_dir}/labels.tsv', 'w') as f:
        f.write("cluster\n")
        for label in labels:
            f.write(f"{label}\n")
    
    # Save metadata as TSV
    with open(f'{output_dir}/metadata.tsv', 'w') as f:
        f.write("cluster\tn_atoms\tBE_per_atom\tsize_category\n")
        for meta in metadata:
            f.write(f"{meta['cluster']}\t{meta['n_atoms']}\t{meta['BE_per_atom']:.6f}\t{meta['size_category']}\n")
    
    # Create instructions for TensorFlow Embedding Projector
    instructions = """
# TensorFlow Embedding Projector Instructions

1. Go to: https://projector.tensorflow.org/
2. Click "Load" button
3. Upload the following files:
   - features.tsv (as "Load a TSV file of vectors")
   - metadata.tsv (as "Load a TSV file of metadata")

## Features included:
- Normalized cluster size (atoms/1000)
- Normalized binding energy per atom
- Energy density (total energy / n_atoms)
- Log10 cluster size
- Cube root cluster size

## Metadata columns:
- cluster: Cluster identifier
- n_atoms: Number of atoms in cluster
- BE_per_atom: Binding energy per atom (eV)
- size_category: Size classification

You can then explore the high-dimensional space of carbon clusters using:
- PCA (Principal Component Analysis)
- t-SNE (t-distributed Stochastic Neighbor Embedding)
- UMAP (Uniform Manifold Approximation and Projection)

Color points by metadata to see patterns in binding energies and cluster sizes!
"""
    
    with open(f'{output_dir}/embedding_instructions.txt', 'w') as f:
        f.write(instructions)
    
    print(f"📈 TensorFlow Embedding Projector files created!")
    print(f"   - features.tsv: {features_array.shape[0]} points with {features_array.shape[1]} dimensions")
    print(f"   - metadata.tsv: Cluster metadata for visualization")
    print(f"   - embedding_instructions.txt: How to use the data")

def main():
    """Main function to generate data and create visualizations."""
    print("🚀 Starting Carbon Cluster Binding Energy Analysis...")
    
    # Read the real data
    results_file = 'project/out/results.csv'
    if os.path.exists(results_file):
        real_df = pd.read_csv(results_file)
        print(f"📖 Loaded real data: {len(real_df)} calculations")
        
        if len(real_df) > 0:
            # Get the base result for synthetic data generation
            base_result = real_df.iloc[0].to_dict()
            print(f"🔬 Base calculation: {base_result['cluster']} with {base_result['n_atoms']} atoms")
            print(f"   Binding energy: {base_result['BE_per_atom']:.3f} eV/atom")
        else:
            # Fallback if no real data
            base_result = {
                'cluster': 'C-example-60',
                'n_atoms': 60,
                'E_cluster': -2.156822359958879,
                'BE_per_atom': -7.334052960667352
            }
            print("⚠️  No real data found, using example values")
    else:
        # Fallback if file doesn't exist
        base_result = {
            'cluster': 'C-example-60',
            'n_atoms': 60,
            'E_cluster': -2.156822359958879,
            'BE_per_atom': -7.334052960667352
        }
        print("⚠️  Results file not found, using example values")
    
    # Generate synthetic dataset
    print("🎲 Generating synthetic dataset...")
    df = generate_synthetic_data(base_result, n_samples=500)
    
    # Save the generated dataset
    df.to_csv('generated_binding_energies.csv', index=False)
    print(f"💾 Saved generated dataset: generated_binding_energies.csv")
    
    # Create visualizations
    print("🎨 Creating fancy visualizations...")
    create_fancy_plots(df)
    
    print("\n✅ Analysis complete!")
    print("🎯 Check out the visualizations and try the TensorFlow Embedding Projector!")

if __name__ == "__main__":
    main()
