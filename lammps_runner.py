import os
import sys
import pandas as pd
from ase import Atoms
from ase.io import read
from ase.calculators.lammpsrun import LAMMPS
import glob
import subprocess

def detect_potential_path():
    """Detect the correct path for the potential file"""
    possible_paths = [
        '/app/project/results/Carbon_GAP_20.xml',  # Docker path
        'project/results/Carbon_GAP_20.xml',       # Local relative path
        os.path.join(os.getcwd(), 'project/results/Carbon_GAP_20.xml'),  # Local absolute path
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    return None

# Try to detect the correct potential path
POTENTIAL_XML_PATH = detect_potential_path()

# Reference energy per atom (eV) for graphite.
E_REF = -7.37

def run_binding_energy(xyz_file, out_dir, lammps_executable):
    """
    Runs a single-point energy calculation using the hybrid Glue+GAP potential.
    """
    atoms = read(xyz_file)
    n_atoms = len(atoms)
    base = os.path.splitext(os.path.basename(xyz_file))[0]
    os.makedirs(out_dir, exist_ok=True)

    # Check that the main potential file exists at its absolute path.
    if not POTENTIAL_XML_PATH or not os.path.exists(POTENTIAL_XML_PATH):
        raise FileNotFoundError(f"Potential file not found. Checked paths: {['/app/project/results/Carbon_GAP_20.xml', 'project/results/Carbon_GAP_20.xml']}")
    
    # Configure LAMMPS parameters with the working hybrid potential format
    params = {
        'pair_style': 'quip',
        'pair_coeff': [f'* * {POTENTIAL_XML_PATH} "IP Glue" 6'],  # Working format for hybrid Glue+GAP
        'mass': ['1 12.011']
    }

    print(f"Running calculation for {base} with {n_atoms} atoms")
    print(f"Using hybrid Glue+GAP potential: {POTENTIAL_XML_PATH}")
    
    # Configure LAMMPS executable
    actual_lammps_executable = lammps_executable
    try:
        # Check if the LAMMPS executable exists
        try:
            subprocess.run(['which', actual_lammps_executable], check=True, capture_output=True)
        except subprocess.CalledProcessError:
            print(f"LAMMPS executable '{actual_lammps_executable}' not found. Trying common alternatives...")
            # Try common LAMMPS executable names, preferring non-MPI versions
            for alt_cmd in ['lmp_serial', 'lammps', 'lmp', 'lmp_mpi']:
                try:
                    subprocess.run(['which', alt_cmd], check=True, capture_output=True)
                    actual_lammps_executable = alt_cmd
                    print(f"Found LAMMPS executable: {actual_lammps_executable}")
                    break
                except subprocess.CalledProcessError:
                    continue
            else:
                raise RuntimeError("No LAMMPS executable found!")
        
        # Handle MPI versions with proper configuration
        if 'mpi' in actual_lammps_executable.lower():
            print(f"Using MPI version {actual_lammps_executable} with container configuration...")
            # Set environment variables to help with MPI in containers
            os.environ['OMPI_ALLOW_RUN_AS_ROOT'] = '1'
            os.environ['OMPI_ALLOW_RUN_AS_ROOT_CONFIRM'] = '1'
            os.environ['OMPI_MCA_btl_vader_single_copy_mechanism'] = 'none'
        
        calc = LAMMPS(
            **params,
            command=actual_lammps_executable,
            keep_tmp_files=True,  # Keep temporary files for debugging
            no_data_file=False,   # Ensure data file is written
            always_triclinic=False  # Use simple orthorhombic box
        )
        
        atoms.calc = calc
        E_cluster = atoms.get_potential_energy()
        BE_per_atom = (n_atoms * E_REF - E_cluster) / n_atoms
        
        print(f"Calculation completed successfully!")
        print(f"Cluster energy: {E_cluster:.6f} eV")
        print(f"Binding energy per atom: {BE_per_atom:.6f} eV")
        
        # Write log file
        log_file = os.path.join(out_dir, base + '.log')
        with open(log_file, 'w') as f:
            f.write(f"TotalEnergy {E_cluster}\\n")
            f.write(f"BE_per_atom {BE_per_atom}\\n")

        return {
            'cluster': base,
            'n_atoms': n_atoms,
            'E_cluster': E_cluster,
            'BE_per_atom': BE_per_atom
        }
        
    except Exception as e:
        print(f"LAMMPS calculation failed: {str(e)}")
        
        # Look for LAMMPS temporary files for debugging
        tmp_dirs = [d for d in os.listdir('/tmp') if d.startswith('LAMMPS-')]
        if tmp_dirs:
            print("\\nDEBUG: Looking for LAMMPS files in /tmp...")
            for tmp_dir in tmp_dirs[-2:]:  # Show last 2 directories
                full_tmp_path = f'/tmp/{tmp_dir}'
                print(f"Found LAMMPS directory: {full_tmp_path}")
                if os.path.exists(full_tmp_path):
                    files = os.listdir(full_tmp_path)
                    print(f"Files: {files}")
                    
                    # Look for log files and print their contents
                    log_files = [f for f in files if f.startswith('log_')]
                    for log_file in log_files:
                        log_path = os.path.join(full_tmp_path, log_file)
                        try:
                            with open(log_path, 'r') as f:
                                log_content = f.read()
                            print(f"\\nContents of {log_file}:")
                            print(log_content[-500:])  # Show last 500 characters
                        except Exception as log_error:
                            print(f"Could not read {log_path}: {log_error}")
        
        raise


def batch_run(cluster_dir, out_dir, results_csv, lammps_executable):
    """
    Processes all .xyz files in cluster_dir, runs binding-energy calc,
    and saves results to results_csv.
    """
    results = []
    xyz_files = sorted(glob.glob(os.path.join(cluster_dir, '*.xyz')))
    
    if not xyz_files:
        print(f"No .xyz files found in directory: {cluster_dir}", file=sys.stderr)
        return pd.DataFrame()

    failed_files = []
    for i, xyz_path in enumerate(xyz_files):
        print(f"\n[{i+1}/{len(xyz_files)}] Running on {os.path.basename(xyz_path)}...")
        try:
            res = run_binding_energy(xyz_path, out_dir, lammps_executable)
            results.append(res)
            print(f"✅ Successfully processed {os.path.basename(xyz_path)}")
        except Exception as e:
            print(f"❌ Failed to process {os.path.basename(xyz_path)}: {e}")
            failed_files.append(xyz_path)
            continue

    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(results_csv), exist_ok=True)
    df.to_csv(results_csv, index=False)
    print(f"\nResults saved to {results_csv}")
    
    # Report summary
    print(f"\n📊 Batch Processing Summary:")
    print(f"   Total files processed: {len(xyz_files)}")
    print(f"   Successful: {len(results)}")
    print(f"   Failed: {len(failed_files)}")
    
    if failed_files:
        print(f"\n❌ Failed files:")
        for f in failed_files:
            print(f"   - {os.path.basename(f)}")
    
    return df

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run GAP-20 binding energy calculations via ASE/LAMMPS in Docker')
    
    parser.add_argument('--xyz', type=str, help='Path to a single XYZ file')
    parser.add_argument('--dir', type=str, default='project/data/clusters', help='Directory of XYZ cluster files')
    parser.add_argument('--out', type=str, default='project/out', help='Output directory for data and logs')
    parser.add_argument('--results', type=str, default='project/out/results.csv', help='CSV path to save aggregated results')
    
    parser.add_argument('--lammps', type=str, default='lmp', help='LAMMPS command')
    
    args = parser.parse_args()

    if args.xyz:
        print(f"Running single file: {args.xyz}")
        if not os.path.exists(args.xyz):
             raise FileNotFoundError(f"Input xyz file not found at {args.xyz}")
        res = run_binding_energy(args.xyz, args.out, args.lammps)
        df = pd.DataFrame([res])
        os.makedirs(os.path.dirname(args.results), exist_ok=True)
        df.to_csv(args.results, index=False)
        print(f"\nResult for {args.xyz}:\n", df)
    else:
        print(f"Running batch mode on directory: {args.dir}")
        if not os.path.isdir(args.dir):
             raise NotADirectoryError(f"Input directory not found at {args.dir}")
        df = batch_run(args.dir, args.out, args.results, args.lammps)
        print("\n--- Batch run complete ---")
        print(df)
