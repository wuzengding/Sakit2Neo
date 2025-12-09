#!/usr/bin/env python
# workflow/scripts/prune_altanalyze.py

import sys
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Prune AltAnalyze results.")
    
    # Input/Output Arguments
    parser.add_argument("--input", required=True, help="Path to counts.original.txt")
    parser.add_argument("--annotation", required=True, help="Path to EventAnnotation.txt")
    parser.add_argument("--output_full", required=True, help="Path to output counts.original.full.txt")
    parser.add_argument("--output_pruned", required=True, help="Path to output counts.original.pruned.txt")
    parser.add_argument("--output_specific", required=True, help="Path to output counts.original.pruned.cancer-specific.txt")
    
    # Filtering Arguments
    parser.add_argument("--min_tumor_count", type=float, default=10.0, help="Minimum read count required in Tumor sample (default: 10)")
    parser.add_argument("--max_normal_count", type=float, default=0.0, help="Maximum read count allowed in Normal sample (default: 0)")
    
    args = parser.parse_args()

    print("Starting Pruning Step...")
    print("Input: " + args.input)

    try:
        # ---------------------------------------------------------
        # 1. Load Data
        # ---------------------------------------------------------
        df = pd.read_csv(args.input, sep='\t', index_col=0)
        
        # 2. Clean Index (Remove '=...')
        df.index = [item.split('=')[0] for item in df.index]
        
        # 3. Remove Duplicates and Save Full
        df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
        df.to_csv(args.output_full, sep='\t')
        print("Full matrix saved to: " + args.output_full)

        # ---------------------------------------------------------
        # 4. Filter by Event Annotation
        # ---------------------------------------------------------
        print("Filtering using annotation: " + args.annotation)
        annot_df = pd.read_csv(args.annotation, sep='\t')
        
        # Extract UIDs (format usually Gene:UID|...)
        valid_uids = []
        for item in annot_df['UID']:
            parts = item.split('|')[0].split(':')
            if len(parts) > 1:
                valid_uids.append(':'.join(parts[1:]))
            else:
                valid_uids.append(item) # Fallback
        
        # Save Pruned File (Contains Tumor AND Normal data)
        df_pruned = df.loc[df.index.isin(set(valid_uids)),:]
        df_pruned.to_csv(args.output_pruned, sep='\t')
        print("Pruned matrix saved to: " + args.output_pruned)

        # ---------------------------------------------------------
        # 5. Filter for Cancer Specificity
        # ---------------------------------------------------------
        print("Filtering for Cancer Specificity...")
        
        # Identify Columns
        tumor_cols = [c for c in df_pruned.columns if 'tumor' in c.lower() or 'Tumor' in c]
        normal_cols = [c for c in df_pruned.columns if 'normal' in c.lower() or 'Normal' in c]
        
        if not tumor_cols or not normal_cols:
            print("WARNING: Could not automatically identify Tumor/Normal columns based on names.")
            print("Columns found: " + str(df_pruned.columns.tolist()))
            print("Skipping specific filtering.")
            # Create an empty file or copy pruned to avoid pipeline crash
            df_pruned.to_csv(args.output_specific, sep='\t')
            sys.exit(0)
            
        print("Tumor Column(s): " + str(tumor_cols))
        print("Normal Column(s): " + str(normal_cols))
        
        # Logic: 
        # Keep row IF (Any Tumor > min) AND (All Normal <= max)
        
        # Get max value across all tumor cols (usually just one)
        max_tumor_vals = df_pruned[tumor_cols].max(axis=1)
        # Get max value across all normal cols
        max_normal_vals = df_pruned[normal_cols].max(axis=1)
        
        # Apply Filter
        condition_tumor = max_tumor_vals >= args.min_tumor_count
        condition_normal = max_normal_vals <= args.max_normal_count
        
        df_specific = df_pruned.loc[condition_tumor & condition_normal, :]
        
        # Save Result
        df_specific.to_csv(args.output_specific, sep='\t')
        print("Cancer Specific matrix saved to: " + args.output_specific)
        print("Events remaining: " + str(len(df_specific)))

    except Exception as e:
        print("Error during pruning: " + str(e))
        sys.exit(1)

if __name__ == "__main__":
    main()