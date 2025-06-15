import os
import sys
import pandas as pd
from ase import Atoms
from ase.io import read, write
import glob

# The path to the main potential file inside the container
GAP_XML = 'project/results/Carbon_GAP_20.xml'

# Reference energy per atom (eV) for graphite
E_REF = -7.37

def run_binding_energy(xyz_file, out_dir, lammps_executable):
    """
    Runs a single-point energy calculation using the GAP-20 potential
    on the cluster defined in xyz_file, computes binding energy per atom,
    writes log, and returns results dict.
    """
    atoms = read(xyz_file)
    n_atoms = len(atoms)
    base = os.path.splitext(os.path.basename(xyz_file))[0]
    os.makedirs(out_dir, exist_ok=True)
    
    # --- KEY FIX ---
    # Get the absolute path of the potential file *inside the container*.
    # The working directory of the container is '/app', so this will resolve to '/app/project/results/Carbon_GAP_20.xml'.
    # This ensures LAMMPS can always find the file, regardless of the temporary directory it runs in.
    abs_gap_xml_path = os.path.abspath(GAP_XML)

    if not os.path.exists(abs_gap_xml_path):
        raise FileNotFoundError(f"Potential file not found at absolute path: {abs_gap_xml_path}")

    params = {
        'pair_style': 'quip',
        # Provide the absolute path to the potential file in the pair_coeff command.
        'pair_coeff': [f'* * {abs_gap_xml_path} "IP GAP" 1'],
        'mass': ['1 12.011']
    }

    from ase.calculators.lammpsrun import LAMMPS
    try:
        print(f"--- Using LAMMPS command: {lammps_executable} ---")
        # We no longer need the 'files' argument, as we provide the full path.
        calc = LAMMPS(**params, command=lammps_executable)
        atoms.calc = calc
        E_cluster = atoms.get_potential_energy()
        BE_per_atom = (n_atoms * E_REF - E_cluster) / n_atoms
    except RuntimeError as e:
        msg = str(e)
        if "FileNotFoundError" in msg or "does not exist" in msg:
             print(f"ERROR: A required file was not found. This might be a potential file or another input.")
        elif (
            "pair style 'quip'" in msg
            or "ML-QUIP package" in msg
            or "Unrecognized pair style 'quip'" in msg
        ):
            print("ERROR: The LAMMPS binary in your environment does not have the ML-QUIP package enabled.")
        else:
            print(f"An unexpected LAMMPS error occurred.")
        print("Original error message:\n", msg)
        sys.exit(1)


    # Write log
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
    for fname in sorted(os.listdir(cluster_dir)):
        if not fname.endswith('.xyz'):
            continue
        xyz_path = os.path.join(cluster_dir, fname)
        print(f"Running on {fname}...")
        res = run_binding_energy(xyz_path, out_dir, lammps_executable)
        results.append(res)

    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(results_csv), exist_ok=True)
    df.to_csv(results_csv, index=False)
    print(f"Results saved to {results_csv}")
    return df

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run GAP-20 binding energy calculations via ASE/LAMMPS in Docker')
    
    parser.add_argument('--xyz', type=str, help='Path to a single XYZ file')
    parser.add_argument('--dir', type=str, default='project/data/clusters', help='Directory of XYZ cluster files')
    parser.add_argument('--out', type=str, default='project/out', help='Output directory for data and logs')
    parser.add_argument('--results', type=str, default='project/out/results.csv', help='CSV path to save aggregated results')
    
    # Call the wrapper script we created in the Dockerfile.
    parser.add_argument('--lammps', type=str, default='run_lammps.sh', help='LAMMPS command')
    
    args = parser.parse_args()

    if args.xyz and not os.path.isdir(args.dir):
        print(f"Running single file: {args.xyz}")
        if not os.path.exists(args.xyz):
             raise FileNotFoundError(f"Input xyz file not found at {args.xyz}")
        res = run_binding_energy(args.xyz, args.out, args.lammps)
        df = pd.DataFrame([res])
        df.to_csv(args.results, index=False)
        print(f"Result for {args.xyz}:\n", df)
    else:
        print(f"Running batch mode on directory: {args.dir}")
        if not os.path.isdir(args.dir):
             raise NotADirectoryError(f"Input directory not found at {args.dir}")
        df = batch_run(args.dir, args.out, args.results, args.lammps)
        print("--- Batch run complete ---")
        print(df)
