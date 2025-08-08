#!/usr/bin/env python3
"""
Pooling Strategy Calculator for Equimolar Sequencing

This script processes HSD1000 compact region table CSV files to calculate
optimal pooling strategies for equimolar sequencing.

Usage: python pooling_strategy.py <input_folder> [--max-samples MAX_SAMPLES]

The script handles both nmol/l and pmol/l molarity units, converting pmol/l to nmol/l internally.
"""

import pandas as pd
import argparse
import os
import sys
from datetime import datetime
from pathlib import Path


def load_and_process_csvs(folder_path):
    """Load all CSV files from folder and extract dimer/library data."""
    all_csv_files = list(Path(folder_path).glob("*.csv"))
    
    # Filter out previous output files (sub-pooling results)
    csv_files = [f for f in all_csv_files if not f.name.endswith('sub-pooling.csv')]
    
    if not csv_files:
        raise ValueError(f"No input CSV files found in {folder_path}")
    
    # Create file-to-plate mapping
    file_to_plate = {}
    print("\nPlate Assignment:")
    for i, csv_file in enumerate(sorted(csv_files), 1):
        file_to_plate[csv_file.name] = i
        print(f"  Plate {i:03d}: {csv_file.name}")
    print()
    
    all_samples = []
    
    for csv_file in csv_files:
        try:
            # Try different encodings to handle special characters like µ
            encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
            df = None
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(csv_file, encoding=encoding)
                    break
                except UnicodeDecodeError:
                    continue
            
            if df is None:
                print(f"Could not decode {csv_file} with any standard encoding")
                continue
            
            # Group by WellId to process each sample
            for well_id, well_data in df.groupby('WellId'):
                sample_data = {
                    'FileName': well_data['FileName'].iloc[0],
                    'Tape Well': well_id,
                    'PlateNumber': file_to_plate[csv_file.name],
                    'Dimer Conc.': None,
                    'Dimer Molarity': None,
                    'Lib Conc.': None,
                    'Lib Molarity': None
                }
                
                dimer_regions = []
                lib_regions = []
                
                # Classify regions as dimer or library
                for _, row in well_data.iterrows():
                    from_bp = row['From [bp]']
                    to_bp = row['To [bp]']
                    
                    # Dimer: shorter fragments (typically 130-160 bp)
                    if from_bp <= 160 and to_bp <= 200:
                        dimer_regions.append(row)
                    # Library: longer fragments (typically 180-350 bp)
                    elif from_bp >= 160:
                        lib_regions.append(row)
                
                # Check for multiple regions in same category
                if len(dimer_regions) > 1:
                    raise ValueError(f"Multiple dimer regions found for {well_id} in {csv_file.name}")
                if len(lib_regions) > 1:
                    raise ValueError(f"Multiple library regions found for {well_id} in {csv_file.name}")
                
                # Find the concentration column (handles different encodings of µ)
                conc_col = None
                molarity_col = None
                molarity_unit = None
                
                for col in df.columns:
                    if 'Conc.' in col:
                        conc_col = col
                    elif 'Region Molarity' in col and ('nmol/l' in col):
                        molarity_col = col
                        molarity_unit = 'nmol'
                    elif 'Region Molarity' in col and ('pmol/l' in col):
                        molarity_col = col
                        molarity_unit = 'pmol'
                
                if conc_col is None or molarity_col is None:
                    print(f"Warning: Could not find concentration/molarity columns in {csv_file}")
                    print(f"Available columns: {list(df.columns)}")
                    continue
                
                # Extract concentrations and molarities with unit conversion
                if dimer_regions:
                    sample_data['Dimer Conc.'] = dimer_regions[0][conc_col]
                    dimer_molarity = dimer_regions[0][molarity_col]
                    # Convert pmol/l to nmol/l if needed
                    if molarity_unit == 'pmol':
                        dimer_molarity = dimer_molarity / 1000
                    sample_data['Dimer Molarity'] = dimer_molarity
                
                if lib_regions:
                    sample_data['Lib Conc.'] = lib_regions[0][conc_col]
                    lib_molarity = lib_regions[0][molarity_col]
                    # Convert pmol/l to nmol/l if needed
                    if molarity_unit == 'pmol':
                        lib_molarity = lib_molarity / 1000
                    sample_data['Lib Molarity'] = lib_molarity
                
                all_samples.append(sample_data)
                
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            continue
    
    return pd.DataFrame(all_samples)


def calculate_target_ratio(df):
    """Calculate target to dimer ratio."""
    df['target ratio'] = (df['Lib Molarity'] / df['Dimer Molarity']).round(2)
    return df


def determine_pool_type(molarity):
    """Determine if pool is strong (>5nmol/l) or weak."""
    # molarity is already in nmol/l
    return "strong" if molarity > 5 else "weak"


def well_id_to_position(well_id):
    """Convert well ID (e.g., 'A1', 'H12') to position number (1-96)."""
    if not well_id or len(well_id) < 2:
        return 0
    
    # Extract row letter and column number
    row_letter = well_id[0].upper()
    try:
        col_number = int(well_id[1:])
    except ValueError:
        return 0
    
    # Convert row letter to number (A=1, B=2, ..., H=8)
    row_number = ord(row_letter) - ord('A') + 1
    
    # Validate ranges
    if row_number < 1 or row_number > 8 or col_number < 1 or col_number > 12:
        return 0
    
    # Calculate position (A1=1, A2=2, ..., A12=12, B1=13, ..., H12=96)
    position = (row_number - 1) * 12 + col_number
    return position


def calculate_pooling_strategy(df, max_samples_per_pool=48):
    """Calculate optimal pooling strategy."""
    # Sort by library molarity (strongest to weakest)
    df_sorted = df.sort_values('Lib Molarity', ascending=False).copy()
    
    # Initialize columns
    df_sorted['sub-pool number'] = 0
    df_sorted['volume added'] = 0.0
    df_sorted['sub-pool volume'] = 0.0
    df_sorted['sub-pool samples'] = 0
    df_sorted['target molarity contribution'] = 0.0
    df_sorted['notes'] = ''
    
    # Initialize new liquid handling columns
    df_sorted['SourcePlateLocation'] = ''
    df_sorted['SourceWellPosition'] = 0
    df_sorted['VolSample'] = 0.0
    df_sorted['BufferLocation'] = 'TEBuffer[001]'
    df_sorted['BufferWellPosition'] = 1
    df_sorted['VolBuffer'] = 0.0
    df_sorted['DestinationPlate'] = 'DestinationPlate[001]'
    df_sorted['DestinationWellPosition'] = 0
    
    current_pool = 1
    unassigned_samples = df_sorted.index.tolist()
    
    while unassigned_samples:
        # Start new sub-pool with strongest remaining sample
        strongest_idx = unassigned_samples[0]
        strongest_molarity = df_sorted.loc[strongest_idx, 'Lib Molarity']
        pool_type = determine_pool_type(strongest_molarity)
        
        # Volume constraints based on pool type
        min_vol_per_sample = 3 if pool_type == "strong" else 10
        max_vol_per_sample = 7 if pool_type == "strong" else 20
        
        # Start with strongest sample - use at least 3μl, preferably more
        strongest_volume = round(max(3, min_vol_per_sample), 2)
        target_moles = strongest_molarity * strongest_volume
        
        current_pool_samples = [strongest_idx]
        current_pool_volume = strongest_volume
        
        # Assign strongest sample
        df_sorted.loc[strongest_idx, 'sub-pool number'] = current_pool
        df_sorted.loc[strongest_idx, 'volume added'] = strongest_volume
        df_sorted.loc[strongest_idx, 'target molarity contribution'] = round(target_moles, 2)
        
        # Check if sample can be pooled (volume constraints)
        if strongest_volume < 1:
            df_sorted.loc[strongest_idx, 'notes'] = 'Too strong - requires <1μl'
        elif strongest_volume > 20:
            df_sorted.loc[strongest_idx, 'notes'] = 'Too weak - requires >20μl'
        
        unassigned_samples.remove(strongest_idx)
        
        # Try to add more samples to this pool
        for sample_idx in unassigned_samples[:]:
            sample_molarity = df_sorted.loc[sample_idx, 'Lib Molarity']
            
            # Calculate required volume for equimolar contribution
            required_volume = round(target_moles / sample_molarity, 2)
            
            # Check if sample fits constraints and is same pool type as the pool
            sample_pool_type = determine_pool_type(sample_molarity)
            sample_min_vol = 3 if sample_pool_type == "strong" else 10
            sample_max_vol = 7 if sample_pool_type == "strong" else 20
            
            # Only add samples of the same pool type (strong with strong, weak with weak)
            # Check if volume is within allowed range, pool won't exceed 150μl, and max samples not reached
            if (sample_pool_type == pool_type and
                sample_min_vol <= required_volume <= sample_max_vol and 
                current_pool_volume + required_volume <= 150 and
                len(current_pool_samples) < max_samples_per_pool):
                
                # Add sample to current pool
                df_sorted.loc[sample_idx, 'sub-pool number'] = current_pool
                df_sorted.loc[sample_idx, 'volume added'] = required_volume
                df_sorted.loc[sample_idx, 'target molarity contribution'] = round(target_moles, 2)
                
                current_pool_samples.append(sample_idx)
                current_pool_volume = round(current_pool_volume + required_volume, 2)
                unassigned_samples.remove(sample_idx)
            
            # Flag samples with volume issues
            elif required_volume < 3:
                df_sorted.loc[sample_idx, 'notes'] = 'Too strong - requires <3μl'
            elif required_volume > 20:
                df_sorted.loc[sample_idx, 'notes'] = 'Too weak - requires >20μl'
        
        # Update pool volume and sample count for all samples in this pool
        for sample_idx in current_pool_samples:
            df_sorted.loc[sample_idx, 'sub-pool volume'] = round(current_pool_volume, 2)
            df_sorted.loc[sample_idx, 'sub-pool samples'] = len(current_pool_samples)
        
        # Check if pool meets minimum volume requirement
        if current_pool_volume < 100:
            for sample_idx in current_pool_samples:
                current_note = df_sorted.loc[sample_idx, 'notes']
                df_sorted.loc[sample_idx, 'notes'] = f"{current_note}; Pool below 100μl minimum".strip('; ')
        
        current_pool += 1
    
    # Populate liquid handling columns
    for idx, row in df_sorted.iterrows():
        # Use stored plate number instead of extracting from filename
        plate_num = row['PlateNumber']
        df_sorted.loc[idx, 'SourcePlateLocation'] = f"SourcePlate[{plate_num:03d}]"
        
        # Convert well ID to position (1-96)
        well_position = well_id_to_position(row['Tape Well'])
        df_sorted.loc[idx, 'SourceWellPosition'] = well_position
        
        # Volume to pool is the same as volume added
        df_sorted.loc[idx, 'VolSample'] = row['volume added']
        
        # Destination well position is the sub-pool number
        df_sorted.loc[idx, 'DestinationWellPosition'] = row['sub-pool number']
    
    return df_sorted


def main():
    parser = argparse.ArgumentParser(description='Calculate pooling strategy for equimolar sequencing')
    parser.add_argument('input_folder', help='Folder containing CSV files to process')
    parser.add_argument('--max-samples', type=int, default=48, 
                       help='Maximum number of samples per sub-pool (default: 48)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_folder):
        print(f"Error: Folder {args.input_folder} does not exist")
        sys.exit(1)
    
    try:
        # Load and process CSV files
        print(f"Processing CSV files in {args.input_folder}...")
        df = load_and_process_csvs(args.input_folder)
        
        if df.empty:
            print("No data found in CSV files")
            sys.exit(1)
        
        # Calculate target ratio
        df = calculate_target_ratio(df)
        
        # Calculate pooling strategy
        df = calculate_pooling_strategy(df, args.max_samples)
        
        # Prepare output
        output_columns = [
            'FileName', 'Tape Well', 'Dimer Conc.', 'Dimer Molarity',
            'Lib Conc.', 'Lib Molarity', 'target ratio', 'sub-pool number',
            'volume added', 'target molarity contribution', 'sub-pool volume', 
            'sub-pool samples', 'notes',
            'SourcePlateLocation', 'SourceWellPosition', 'VolSample',
            'BufferLocation', 'BufferWellPosition', 'VolBuffer',
            'DestinationPlate', 'DestinationWellPosition'
        ]
        
        # Ensure all columns exist
        for col in output_columns:
            if col not in df.columns:
                df[col] = ''
        
        df_output = df[output_columns].copy()
        
        # Sort by sub-pool number and library molarity
        df_output = df_output.sort_values(['sub-pool number', 'Lib Molarity'], ascending=[True, False])
        
        # Generate output filename with timestamp and create output folder
        output_folder = os.path.join(args.input_folder, "output")
        os.makedirs(output_folder, exist_ok=True)
        
        timestamp = datetime.now().strftime('%Y-%m-%d_%H%M%S')
        output_file = os.path.join(output_folder, f"{timestamp}_sub-pooling.csv")
        
        # Save to CSV
        df_output.to_csv(output_file, index=False, encoding='utf-8')
        
        print(f"Pooling strategy saved to: {output_file}")
        print(f"Processed {len(df)} samples into {df['sub-pool number'].max()} sub-pools")
        
        # Print summary
        pool_summary = df.groupby('sub-pool number').agg({
            'sub-pool volume': 'first',
            'sub-pool samples': 'first'
        })
        
        print("\nSub-pool Summary:")
        for pool_num, row in pool_summary.iterrows():
            print(f"Pool {pool_num}: {row['sub-pool samples']} samples, {row['sub-pool volume']:.1f}ul total")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
