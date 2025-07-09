#!/usr/bin/env python3
"""
Docker wrapper script to run LAMMPS calculations on all xyz files.
This script builds the Docker image and runs the calculations inside the container.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path

def build_docker_image():
    """Build the Docker image for LAMMPS calculations"""
    print("🐳 Building Docker image...")
    
    try:
        result = subprocess.run([
            'docker', 'build', '-t', 'kurfcarla_lammps', '.'
        ], check=True, capture_output=True, text=True)
        print("✅ Docker image built successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Docker build failed:")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return False

def run_docker_lammps(xyz_dir, out_dir, results_file, lammps_cmd='lmp'):
    """Run LAMMPS calculations in Docker container"""
    
    # Get absolute paths
    current_dir = os.getcwd()
    xyz_dir_abs = os.path.abspath(xyz_dir)
    out_dir_abs = os.path.abspath(out_dir)
    results_file_abs = os.path.abspath(results_file)
    
    # Ensure output directory exists
    os.makedirs(out_dir_abs, exist_ok=True)
    
    print(f"🚀 Running LAMMPS calculations in Docker...")
    print(f"   XYZ directory: {xyz_dir_abs}")
    print(f"   Output directory: {out_dir_abs}")
    print(f"   Results file: {results_file_abs}")
    
    # Docker run command
    docker_cmd = [
        'docker', 'run', '--rm',
        '-v', f'{current_dir}:/app',  # Mount entire project directory
        '-v', f'{xyz_dir_abs}:/app/project/data/clusters',  # Mount xyz files
        '-v', f'{out_dir_abs}:/app/project/out',  # Mount output directory
        'kurfcarla_lammps',
        'python', 'lammps_runner.py',
        '--dir', '/app/project/data/clusters',
        '--out', '/app/project/out',
        '--results', '/app/project/out/results.csv',
        '--lammps', lammps_cmd
    ]
    
    try:
        print(f"Running command: {' '.join(docker_cmd)}")
        result = subprocess.run(docker_cmd, check=True, text=True)
        print("✅ LAMMPS calculations completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Docker run failed:")
        print(f"Return code: {e.returncode}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Run LAMMPS calculations in Docker for all xyz files')
    parser.add_argument('--xyz-dir', type=str, default='project/data/clusters', 
                       help='Directory containing xyz files')
    parser.add_argument('--out-dir', type=str, default='project/out', 
                       help='Output directory for results')
    parser.add_argument('--results', type=str, default='project/out/results.csv', 
                       help='Results CSV file')
    parser.add_argument('--lammps', type=str, default='lmp', 
                       help='LAMMPS command inside container')
    parser.add_argument('--rebuild', action='store_true', 
                       help='Force rebuild of Docker image')
    
    args = parser.parse_args()
    
    # Check if xyz directory exists
    if not os.path.exists(args.xyz_dir):
        print(f"❌ XYZ directory not found: {args.xyz_dir}")
        sys.exit(1)
    
    # Count xyz files
    xyz_files = list(Path(args.xyz_dir).glob('*.xyz'))
    print(f"📊 Found {len(xyz_files)} xyz files to process")
    
    if len(xyz_files) == 0:
        print("❌ No xyz files found!")
        sys.exit(1)
    
    # Build Docker image (only if it doesn't exist or rebuild is requested)
    if args.rebuild:
        if not build_docker_image():
            sys.exit(1)
    else:
        # Check if image exists
        try:
            result = subprocess.run(['docker', 'image', 'inspect', 'kurfcarla_lammps'], 
                                  capture_output=True, check=True)
            print("✅ Docker image already exists")
        except subprocess.CalledProcessError:
            print("🔨 Docker image not found, building...")
            if not build_docker_image():
                sys.exit(1)
    
    # Run LAMMPS calculations
    if run_docker_lammps(args.xyz_dir, args.out_dir, args.results, args.lammps):
        print(f"\n🎉 All calculations completed!")
        print(f"Results saved to: {os.path.abspath(args.results)}")
        
        # Show results summary if file exists
        if os.path.exists(args.results):
            try:
                import pandas as pd
                df = pd.read_csv(args.results)
                print(f"\n📈 Results Summary:")
                print(f"   Total clusters processed: {len(df)}")
                print(f"   Atom count range: {df['n_atoms'].min()} - {df['n_atoms'].max()}")
                print(f"   BE per atom range: {df['BE_per_atom'].min():.3f} - {df['BE_per_atom'].max():.3f} eV")
                print(f"   Mean BE per atom: {df['BE_per_atom'].mean():.3f} ± {df['BE_per_atom'].std():.3f} eV")
            except Exception as e:
                print(f"Could not read results file: {e}")
    else:
        print("❌ Calculations failed!")
        sys.exit(1)

if __name__ == '__main__':
    main()
