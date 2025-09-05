#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filters a VCF file based on a list of selected variants from a review Excel report.

This script reads the "Somatic & Noise Event Report" sheet, identifies all
variants marked with 'yes' in the 'Manual_Select' column, and creates a new,
smaller VCF file containing only the full records for these selected variants.
The output VCF is indexed and ready for phasing analysis.
"""

import pandas as pd
import pysam
import argparse
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report_xlsx', required=True, help='Input Excel report file.')
    parser.add_argument('--input_vcf', required=True, help='The original VEP-annotated VCF file.')
    parser.add_argument('--output_vcf', required=True, help='Path for the output filtered VCF file.')
    args = parser.parse_args()

    # --- Step 1: Read the Excel report and get selected variants ---
    logging.info(f"Reading selected variants from {args.report_xlsx}...")
    try:
        df = pd.read_excel(args.report_xlsx, sheet_name="Somatic & Noise Event Report")
        df['Manual_Select'] = df['Manual_Select'].astype(str).str.lower()
        selected_df = df[df['Manual_Select'] == 'yes']
    except Exception as e:
        logging.error(f"Failed to read Excel file: {e}")
        sys.exit(1)

    if selected_df.empty:
        logging.warning("No variants marked 'yes' found in the report. Creating an empty VCF.")
        # Create an empty VCF with just the header
        with pysam.VariantFile(args.input_vcf) as vcf_in:
            with pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header) as vcf_out:
                pass # Just write the header
        pysam.tabix_index(args.output_vcf, preset="vcf")
        sys.exit(0)
    
    # Create a set of unique variant coordinates for fast lookup.
    # Key is a tuple: (chrom, pos, ref, alt)
    target_variants = set()
    for _, row in selected_df.iterrows():
        # Pysam positions are 1-based, so no conversion needed from Excel
        target_variants.add((str(row['CHROM']), int(row['POS']), str(row['REF']), str(row['ALT'])))

    logging.info(f"Identified {len(target_variants)} unique variants to select.")

    # --- Step 2: Iterate through the input VCF and write selected records ---
    logging.info(f"Filtering {args.input_vcf} to create {args.output_vcf}...")
    records_written = 0
    with pysam.VariantFile(args.input_vcf) as vcf_in:
        # Create the output VCF with the same header as the input
        with pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                # Check if any allele in this record is in our target set
                is_target_record = False
                for alt_allele in record.alts:
                    variant_tuple = (record.chrom, record.pos, record.ref, alt_allele)
                    if variant_tuple in target_variants:
                        is_target_record = True
                        break # Found a match, no need to check other alleles
                
                if is_target_record:
                    vcf_out.write(record)
                    records_written += 1

    logging.info(f"Wrote {records_written} records to the filtered VCF.")

    # --- Step 3: Index the output VCF for phasing tools ---
    logging.info(f"Indexing {args.output_vcf} with tabix...")
    pysam.tabix_index(args.output_vcf, preset="vcf")

    logging.info("Filtering and indexing complete.")


if __name__ == "__main__":
    main()