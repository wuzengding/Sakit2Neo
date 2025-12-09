#!/usr/bin/env python
# workflow/scripts/prune_altanalyze.py

import sys
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Prune AltAnalyze results.")
    parser.add_argument("--input", required=True, help="Path to counts.original.txt")
    parser.add_argument("--annotation", required=True, help="Path to EventAnnotation.txt")
    parser.add_argument("--output_full", required=True, help="Path to output counts.original.full.txt")
    parser.add_argument("--output_pruned", required=True, help="Path to output counts.original.pruned.txt")
    
    args = parser.parse_args()

    print("Starting Pruning Step...")
    print("Input: " + args.input)

    try:
        # 1. Load Data
        df = pd.read_csv(args.input, sep='\t', index_col=0)
        
        # 2. Clean Index (Remove '=...')
        df.index = [item.split('=')[0] for item in df.index]
        
        # 3. Remove Duplicates and Save Full
        df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
        df.to_csv(args.output_full, sep='\t')
        print("Full matrix saved to: " + args.output_full)

        # 4. Filter by Event Annotation
        print("Filtering using annotation: " + args.annotation)
        annot_df = pd.read_csv(args.annotation, sep='\t')
        
        # Extract UIDs (format usually Gene:UID|...)
        # We handle potential errors if format differs slightly
        valid_uids = []
        for item in annot_df['UID']:
            parts = item.split('|')[0].split(':')
            if len(parts) > 1:
                valid_uids.append(':'.join(parts[1:]))
            else:
                valid_uids.append(item) # Fallback
        
        # Filter
        df_pruned = df.loc[df.index.isin(set(valid_uids)),:]
        df_pruned.to_csv(args.output_pruned, sep='\t')
        
        print("Success. Pruned file saved to: " + args.output_pruned)

    except Exception as e:
        print("Error during pruning: " + str(e))
        sys.exit(1)

if __name__ == "__main__":
    main()