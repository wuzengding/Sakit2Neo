#!/usr/bin/env python3
import argparse
import sys
import re
import gzip
import subprocess
from pathlib import Path


class VCFReformatter:
    def __init__(self):
        self.header_replacements = {
            '##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">': 
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
            '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">': 
                None,  # Will be skipped
            '##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">': 
                '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">': 
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
            '##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">': 
                '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher\'s Exact Test to detect strand bias.">'
        }

    def reformat_vcf(self, input_vcf: Path, output_vcf: Path):
        """Reformat the VCF file according to the specified requirements."""
        with input_vcf.open('r') as infile, output_vcf.open('w') as outfile:
            for line in infile:
                if line.startswith('##'):
                    # Process header lines
                    replaced = False
                    for old, new in self.header_replacements.items():
                        if old in line:
                            if new is not None:
                                outfile.write(new + '\n')
                            replaced = True
                            break
                    if not replaced:
                        outfile.write(line)
                
                elif line.startswith('#CHROM'):
                    # Keep the header line as is
                    outfile.write(line)
                
                else:
                    # Process data lines
                    fields = line.strip().split('\t')
                    
                    # Process INFO field
                    info = fields[7]
                    if info.startswith('DP=') and 'SOMATIC' not in info:
                        info_parts = [p for p in info.split(';') if not p.startswith('DP=')]
                        info = ';'.join(info_parts)
                    fields[7] = info
                    
                    # Process FORMAT and sample fields
                    if len(fields) > 8:
                        fields = self._process_sample_fields(fields)
                    
                    outfile.write('\t'.join(fields) + '\n')

    def _process_sample_fields(self, fields):
        """Process FORMAT and sample columns."""
        format_fields = fields[8].split(':')
        samples = fields[9:]
        
        # Get indices of old fields
        try:
            rd_idx = format_fields.index('RD')
            ad_idx = format_fields.index('AD')
            freq_idx = format_fields.index('FREQ')
            dp4_idx = format_fields.index('DP4')
        except ValueError as e:
            raise ValueError(f"Missing required FORMAT field in VCF: {e}")

        # Define new format fields and process each sample
        new_format = []
        new_samples = []
        
        # Keep GT, GQ, DP if present
        for fmt in ['GT', 'GQ', 'DP']:
            if fmt in format_fields:
                fmt_idx = format_fields.index(fmt)
                new_format.append(fmt)
                for i, sample in enumerate(samples):
                    if i >= len(new_samples):
                        new_samples.append([])
                    new_samples[i].append(sample.split(':')[fmt_idx])
        
        # Add new AD field (combine RD and AD)
        new_format.append('AD')
        for i, sample in enumerate(samples):
            sample_fields = sample.split(':')
            rd = sample_fields[rd_idx] if rd_idx < len(sample_fields) else '0'
            ad = sample_fields[ad_idx] if ad_idx < len(sample_fields) else '0'
            new_samples[i].append(f"{rd},{ad}")
        
        # Add AF field (convert FREQ to float)
        new_format.append('AF')
        for i, sample in enumerate(samples):
            sample_fields = sample.split(':')
            freq = sample_fields[freq_idx].replace('%', '') if freq_idx < len(sample_fields) else '0'
            try:
                af = float(freq) / 100
                new_samples[i].append(f"{af:.4f}")
            except ValueError:
                new_samples[i].append('0.0')
        
        # Add SB field (same as DP4 but properly typed)
        new_format.append('SB')
        for i, sample in enumerate(samples):
            sample_fields = sample.split(':')
            dp4 = sample_fields[dp4_idx].split(',') if dp4_idx < len(sample_fields) else ['0'] * 4
            new_samples[i].append(','.join(dp4[:4]))  # Ensure exactly 4 values
        
        # Update the fields
        fields[8] = ':'.join(new_format)
        fields[9:] = [':'.join(sample) for sample in new_samples]
        
        return fields

    def compress_and_index_vcf(self, input_vcf: Path, output_vcf_gz: Path):
        """Compress with bgzip and index with tabix."""
        try:
            # 使用 bgzip 压缩（而非 Python 的 gzip）
            subprocess.run(
                ['bgzip', '-f', '-c', str(input_vcf)],  # -c: 输出到标准输出
                check=True,
                stdout=open(output_vcf_gz, 'wb'),  # 重定向到文件
                stderr=subprocess.PIPE
            )
        
            # 检查文件是否为空
            if output_vcf_gz.stat().st_size == 0:
                raise RuntimeError("Compressed VCF file is empty after bgzip")

            # 索引
            subprocess.run(
                ['tabix', '-p', 'vcf', str(output_vcf_gz)],
                check=True,
                stderr=subprocess.PIPE
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to process VCF: {e.stderr.decode()}") from e



def main():
    parser = argparse.ArgumentParser(
        description="Reformat VarScan2 VCF to standardized format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'input_vcf',
        type=Path,
        help='Input VCF file to reformat'
    )
    parser.add_argument(
        'output_vcf',
        type=Path,
        help='Output VCF file (uncompressed)'
    )
    parser.add_argument(
        'output_vcf_gz',
        type=Path,
        help='Output compressed VCF file (.gz)'
    )
    parser.add_argument(
        '--keep-intermediate',
        action='store_true',
        help='Keep the intermediate uncompressed VCF file'
    )
    
    args = parser.parse_args()
    
    reformatter = VCFReformatter()
    
    try:
        # Step 1: Reformat the VCF
        print(f"Reformatting {args.input_vcf}...")
        reformatter.reformat_vcf(args.input_vcf, args.output_vcf)
        
        # Step 2: Compress and index
        print(f"Compressing and indexing to {args.output_vcf_gz}...")
        reformatter.compress_and_index_vcf(args.output_vcf, args.output_vcf_gz)
        
        # Step 3: Clean up if requested
        if not args.keep_intermediate:
            print(f"Removing intermediate file {args.output_vcf}...")
            args.output_vcf.unlink()
        
        print("Successfully completed!")
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()