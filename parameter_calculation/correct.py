#!/usr/bin/env python3
"""
Script Purpose: Correct Dxy values calculated from SNP-only VCF files.

Background: 
When using popgenWindows.py on VCFs without invariant sites, the denominator 
for Dxy defaults to the number of SNPs rather than the total physical window size.
This script applies the standard correction: 
Corrected_Dxy = Raw_Dxy * (Number_of_SNPs / Physical_Window_Size).
"""

import pandas as pd
import glob
import os

# Configuration
INPUT_PATTERN = "raw_dxy_results.csv.gz"

def correct_dxy_stats():
    file_list = glob.glob(INPUT_PATTERN)
    
    if not file_list:
        print(f"[ERROR] No files matching '{INPUT_PATTERN}' found. Please check Bash script output.")
        return

    print(f"[INFO] Found {len(file_list)} file(s). Starting correction...\n")

    for f in file_list:
        print(f"Processing: {f} ...")
        
        try:
            df = pd.read_csv(f)
            
            if 'sites' not in df.columns:
                print(f"  [WARNING] Skipping: File lacks the required 'sites' column.")
                continue
            
            if 'start' in df.columns and 'end' in df.columns:
                window_sizes = df['end'] - df['start']
            else:
                print("  [WARNING] 'start'/'end' columns not found. Defaulting to a window size of 50,000.")
                window_sizes = 50000
            
            cols_to_fix = [c for c in df.columns if c.startswith('dxy_')]
            
            if not cols_to_fix:
                print(f"  [WARNING] No 'dxy_' columns found to correct.")
                continue
            
            print(f"  -> Columns to be corrected: {cols_to_fix}")

            correction_factor = df['sites'] / window_sizes
            
            for col in cols_to_fix:
                corrected_values = df[col] * correction_factor
                df[col] = corrected_values.round(5)
            
            out_name = f.replace(".gz", "") + ".corrected.csv"
            if out_name == f:
                 out_name = f + ".corrected.csv"
            
            df.to_csv(out_name, index=False, float_format='%.5f')
            print(f"  [SUCCESS] Saved corrected output to: {out_name}\n")

        except Exception as e:
            print(f"  [ERROR] Failed to process {f}. Reason: {e}\n")

if __name__ == "__main__":
    correct_dxy_stats()