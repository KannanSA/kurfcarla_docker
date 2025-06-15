import os
import sys
import pandas as pd
from ase import Atoms
from ase.io import read
from ase.calculators.lammpsrun import LAMMPS
import glob

# The absolute path to the main potential file inside the container.
POTENTIAL_XML_PATH = '/app/project/results/Carbon_GAP_20.xml'

# Reference energy per atom (eV) for graphite.
E_REF = -7.37

def run_binding_energy(xyz_file, out_dir, lammps_executable):
    """
    Runs a single-point energy calculation using an absolute path for the potential.
    """
    atoms = read(xyz_file)
    n_atoms = len(atoms)
    base = os.path.splitext(os.path.basename(xyz_file))[0]
    os.makedirs(out_dir, exist_ok=True)

    # Check that the main potential file exists at its absolute path.
    if not os.path.exists(POTENTIAL_XML_PATH):
        raise FileNotFoundError(f"Potential file not found at the absolute path: {POTENTIAL_XML_PATH}")

    params = {
        'pair_style': 'quip',
        # Provide the absolute path to the potential file. LAMMPS will find the
        # companion .sparseX files in the same directory.
        'pair_coeff': [f'* * {POTENTIAL_XML_PATH} "IP GAP" 1'],
        'mass': ['1 12.011']
    }

    try:
        print(f"--- Using LAMMPS command: {lammps_executable} ---")
        
        # We no longer need the 'files' argument since we are using an absolute path.
        calc = LAMMPS(
            **params,
            command=lammps_executable
        )
        
        atoms.calc = calc
        E_cluster = atoms.get_potential_energy()
        BE_per_atom = (n_atoms * E_REF - E_cluster) / n_atoms
        
    except RuntimeError as e:
        print("An unexpected LAMMPS error occurred.", file=sys.stderr)
        print("Original error message:\n", e, file=sys.stderr)
        sys.exit(1)

    # Write log.
    log_file = os.path.join(out_dir, base + '.log')
    with open(log_file, 'w') as f:
        f.write(f"TotalEnergy {E_cluster}\n")
        f.write(f"BE_per_atom {BE_per_atom}\n")

    return {
        'cluster': base,
        'n_atoms': n_atoms,
        'E_cluster': E_cluster,
        'BE_per_atom': BE_per_atom
    }

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

    for xyz_path in xyz_files:
        print(f"Running on {os.path.basename(xyz_path)}...")
        res = run_binding_energy(xyz_path, out_dir, lammps_executable)
        results.append(res)

    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(results_csv), exist_ok=True)
    df.to_csv(results_csv, index=False)
    print(f"\nResults saved to {results_csv}")
    return df

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run GAP-20 binding energy calculations via ASE/LAMMPS in Docker')
    
    parser.add_argument('--xyz', type=str, help='Path to a single XYZ file')
    parser.add_argument('--dir', type=str, default='project/data/clusters', help='Directory of XYZ cluster files')
    parser.add_argument('--out', type=str, default='project/out', help='Output directory for data and logs')
    parser.add_argument('--results', type=str, default='project/out/results.csv', help='CSV path to save aggregated results')
    
    parser.add_argument('--lammps', type=str, default='lmp_mpi', help='LAMMPS command')
    
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