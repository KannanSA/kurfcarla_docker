#!/usr/bin/env python3
"""
Script to find and move all .xyz files into the clusters folder.
Recursively searches for .xyz files and moves them to project/data/clusters/
"""

import os
import shutil
import glob
from pathlib import Path
import hashlib

def calculate_file_hash(file_path):
    """Calculate MD5 hash of a file to check for true duplicates"""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def find_and_move_xyz_files(search_dir='.', target_dir='project/data/clusters'):
    """
    Find all .xyz files in search_dir and move them to target_dir.
    
    Args:
        search_dir (str): Directory to search for .xyz files (default: current directory)
        target_dir (str): Target directory to move files to (default: project/data/clusters)
    """
    
    # Convert to Path objects for easier handling
    search_path = Path(search_dir).resolve()
    target_path = Path(target_dir).resolve()
    
    # Create target directory if it doesn't exist
    target_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Searching for .xyz files in: {search_path}")
    print(f"Target directory: {target_path}")
    print("-" * 50)
    
    # Find all .xyz files recursively
    xyz_files = list(search_path.rglob('*.xyz'))
    
    if not xyz_files:
        print("No .xyz files found.")
        return
    
    print(f"Found {len(xyz_files)} .xyz file(s):")
    
    moved_count = 0
    skipped_count = 0
    
    for xyz_file in xyz_files:
        # Skip if file is already in target directory
        if xyz_file.parent == target_path:
            print(f"  SKIP: {xyz_file.name} (already in target directory)")
            skipped_count += 1
            continue
        
        # Create destination path
        dest_file = target_path / xyz_file.name
        
        # Handle name conflicts
        if dest_file.exists():
            # If files are identical, skip
            if xyz_file.stat().st_size == dest_file.stat().st_size:
                print(f"  SKIP: {xyz_file.name} (identical file exists in target)")
                skipped_count += 1
                continue
            else:
                # Create unique name
                counter = 1
                base_name = xyz_file.stem
                suffix = xyz_file.suffix
                while dest_file.exists():
                    dest_file = target_path / f"{base_name}_{counter}{suffix}"
                    counter += 1
                print(f"  RENAME: {xyz_file.name} -> {dest_file.name} (name conflict)")
        
        try:
            # Move the file
            shutil.move(str(xyz_file), str(dest_file))
            print(f"  MOVED: {xyz_file.name} -> {dest_file}")
            moved_count += 1
        except Exception as e:
            print(f"  ERROR: Failed to move {xyz_file.name}: {e}")
    
    print("-" * 50)
    print(f"Summary:")
    print(f"  Files moved: {moved_count}")
    print(f"  Files skipped: {skipped_count}")
    print(f"  Total processed: {len(xyz_files)}")

def remove_duplicate_files(target_dir='project/data/clusters', pattern='*_1.xyz'):
    """
    Remove duplicate files ending in _1.xyz (or other pattern) from target directory.
    
    Args:
        target_dir (str): Directory to search for duplicate files
        pattern (str): Pattern to match duplicate files (default: *_1.xyz)
    """
    target_path = Path(target_dir).resolve()
    
    if not target_path.exists():
        print(f"Target directory does not exist: {target_path}")
        return
    
    print(f"Searching for duplicate files matching '{pattern}' in: {target_path}")
    print("-" * 50)
    
    # Find all files matching the pattern
    duplicate_files = list(target_path.glob(pattern))
    
    if not duplicate_files:
        print(f"No files matching pattern '{pattern}' found.")
        return
    
    print(f"Found {len(duplicate_files)} file(s) matching pattern:")
    
    removed_count = 0
    skipped_count = 0
    
    for dup_file in duplicate_files:
        # Remove file regardless of size or content
        try:
            dup_file.unlink()
            print(f"  REMOVED: {dup_file.name}")
            removed_count += 1
        except Exception as e:
            print(f"  ERROR: Failed to remove {dup_file.name}: {e}")
            skipped_count += 1
    
    print("-" * 50)
    print(f"Summary:")
    print(f"  Files removed: {removed_count}")
    print(f"  Files skipped: {skipped_count}")
    print(f"  Total processed: {len(duplicate_files)}")

def main():
    """Main function with command line argument support"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Find and move all .xyz files to clusters folder, or remove duplicates')
    parser.add_argument('--search', '-s', type=str, default='.', 
                       help='Directory to search for .xyz files (default: current directory)')
    parser.add_argument('--target', '-t', type=str, default='project/data/clusters',
                       help='Target directory to move files to (default: project/data/clusters)')
    parser.add_argument('--dry-run', '-n', action='store_true',
                       help='Show what would be moved without actually moving files')
    parser.add_argument('--remove-duplicates', '-r', action='store_true',
                       help='Remove duplicate files ending in _1.xyz')
    parser.add_argument('--pattern', '-p', type=str, default='*_1.xyz',
                       help='Pattern to match duplicate files (default: *_1.xyz)')
    
    args = parser.parse_args()
    
    if args.remove_duplicates:
        print("DUPLICATE REMOVAL MODE")
        print("=" * 50)
        remove_duplicate_files(args.target, args.pattern)
    elif args.dry_run:
        print("DRY RUN MODE - No files will be moved")
        print("=" * 50)
        
        # Just find and list files
        search_path = Path(args.search).resolve()
        target_path = Path(args.target).resolve()
        
        print(f"Would search in: {search_path}")
        print(f"Would move to: {target_path}")
        print("-" * 50)
        
        xyz_files = list(search_path.rglob('*.xyz'))
        
        if not xyz_files:
            print("No .xyz files found.")
            return
        
        print(f"Found {len(xyz_files)} .xyz file(s) that would be moved:")
        for xyz_file in xyz_files:
            if xyz_file.parent == target_path:
                print(f"  SKIP: {xyz_file.name} (already in target directory)")
            else:
                print(f"  WOULD MOVE: {xyz_file.relative_to(search_path)} -> {target_path / xyz_file.name}")
    else:
        find_and_move_xyz_files(args.search, args.target)
        remove_duplicate_files(args.target)

if __name__ == '__main__':
    main()
