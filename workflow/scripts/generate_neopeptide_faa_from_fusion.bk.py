#!/usr/bin/env python3
import csv
import argparse
import sys

def clean_sequence(seq):
    """Removes stop codon asterisks from the sequence."""
    return seq.replace('*', '')

def get_ref_part(left_12):
    """
    Processes the Gene1 part for the Reference sequence.
    Arriba uses lowercase letters for junction mismatches. In the Reference 
    sequence, we replace these non-wildtype residues with 'A' to maintain 
    the alignment with the Variant sequence.
    """
    res = ""
    for char in left_12:
        if char.islower():
            res += "A"
        else:
            res += char
    return res

def process_fusion_peptides(input_tsv, output_faa, confidence_levels):
    """
    Parses Arriba output and generates neo-peptides with Ref/Var pairs.
    """
    target_conf = [c.lower() for c in confidence_levels]
    
    try:
        with open(input_tsv, 'r', encoding='utf-8') as tsv_file, \
             open(output_faa, 'w', encoding='utf-8') as faa_file:
            
            # Arriba TSV files use '#gene1' as the header for the first column
            reader = csv.DictReader(tsv_file, delimiter='\t')
            entry_count = 0

            for row in reader:
                # 1. Filter by Confidence level
                confidence = row.get('confidence', '').lower()
                if confidence not in target_conf:
                    continue
                
                gene1 = row.get('#gene1', 'Unknown')
                gene2 = row.get('gene2', 'Unknown')
                peptide_raw = row.get('peptide_sequence', '')
                
                # Process only if sequence is present and contains the fusion marker '|'
                if not peptide_raw or '|' not in peptide_raw:
                    continue
                
                # 2. Split and Clean sequence parts
                parts = peptide_raw.split('|')
                left_full = parts[0]
                right_full = clean_sequence(parts[1])
                
                # Logic: Extract last 12aa of Gene1 (left part)
                left_12 = left_full[-12:] if len(left_full) >= 12 else left_full
                
                # Check for out-of-frame (presence of lowercase in original right part)
                is_out_of_frame = any(c.islower() for c in parts[1])
                
                var_seq_right = ""
                if not is_out_of_frame:
                    # In-frame: Take first 12aa of the right gene
                    var_seq_right = right_full[:12]
                else:
                    # Out-of-frame: Take all characters of the novel/processed peptide sequence
                    var_seq_right = right_full
                
                # 3. Construct Final Sequences
                # Var: Merge parts and convert all residues (including junction mismatches) to uppercase
                var_seq = (left_12 + var_seq_right).upper()
                
                # Ref: Clean Gene1 lowercase residues and pad right side with 'A's to match Var length
                ref_left = get_ref_part(left_12)
                ref_seq = ref_left + ('A' * len(var_seq_right))
                
                # 4. Generate fusion ID (e.g., fusion001) and write 4-line FASTA format
                entry_count += 1
                fusion_id = f"fusion{entry_count:03d}"
                
                faa_file.write(f">Ref|-\n{ref_seq}\n")
                faa_file.write(f">Var|-|{gene1}-{gene2}|{fusion_id}|-\n{var_seq}\n")
                
        print(f"Success: {entry_count} entries processed into '{output_faa}'.")

    except FileNotFoundError:
        print(f"Error: Input file '{input_tsv}' not found.")
    except Exception as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert Arriba fusion results to a neo-peptide FASTA file."
    )
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True, 
                        help="Input Arriba fusion TSV file")
    parser.add_argument("-o", "--output", required=True, 
                        help="Output FAA file path")
    
    # Optional arguments
    parser.add_argument("-c", "--confidence", nargs='+', default=["high", "medium"],
                        help="Confidence levels to extract (e.g., high medium). Default: high medium")

    args = parser.parse_args()
    print(args)
    process_fusion_peptides(args.input, args.output, args.confidence)

if __name__ == "__main__":
    main()
