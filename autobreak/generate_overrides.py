"""
generate_sequence_overrides.py

Generate sequence overrides to optimize problematic DNA strands based on thermodynamic analysis.
This script reads the strand_data.xlsx file and generates an -overrides.txt file with suggested
mutations to improve strand stability.

Output format (compatible with parse_override_input in apply_overrides script):

    # Comments starting with '#'
    # ...
    vh_num[idx][direction]: SEQUENCE

where:
    - vh_num  = HelixID (integer)
    - idx     = StartIdx (integer, cadnano index along that helix)
    - direction = 'fwd' or 'rev'
    - SEQUENCE  = optimized DNA sequence (A/T/C/G), spaces allowed
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import sys


class SequenceOptimizer:
    """
    Analyzes strand thermodynamic data and generates sequence overrides
    to improve problematic regions.
    """
    
    def __init__(self, strand_data_path, output_path=None):
        """
        Initialize the optimizer.
        
        Args:
            strand_data_path: Path to the strand_data.xlsx file
            output_path: Path for output overrides file (default: auto-generate)
        """
        self.strand_data_path = Path(strand_data_path)
        
        if output_path:
            self.output_path = Path(output_path)
        else:
            # Auto-generate output path in same directory
            base_name = self.strand_data_path.stem.replace('_strand_data', '')
            self.output_path = self.strand_data_path.parent / f'{base_name}_overrides.txt'
        
        # Load data
        self.df = pd.read_excel(strand_data_path, sheet_name='StrandData')
        
        print(f"Loaded {len(self.df)} strands from {strand_data_path}")
        
        # Ensure numeric columns are properly typed
        numeric_cols = [
            'Tf', 'ProbFold', 'dGtotal', 'dGhyb', 'dGloop', 'dGconc',
            'Length', 'HelixID', 'StartIdx', 'EndIdx', 'maxTm'
        ]
        for col in numeric_cols:
            if col in self.df.columns:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')
        
        # Show Tf statistics
        if 'Tf' in self.df.columns:
            print(f"  Tf range: {self.df['Tf'].min():.1f} to {self.df['Tf'].max():.1f}")
            print(f"  Tf mean: {self.df['Tf'].mean():.1f}")
            print(f"  Strands with Tf < 45: {len(self.df[self.df['Tf'] < 45.0])}")
    
    def identify_problematic_strands(self, tf_threshold=45.0, use_dg=False,
                                     dg_threshold=-5.0):
        """
        Identify strands that need optimization based on thermodynamic criteria.
        
        Args:
            tf_threshold: Minimum Tf (folding temperature)
            use_dg: If True, also filter by dG threshold (default: False)
            dg_threshold: Maximum dG in kcal/mol (only used if use_dg=True)
            
        Returns:
            DataFrame of problematic strands sorted by Tf (lowest first)
        """
        # Primary filter: Tf only
        problematic = self.df[self.df['Tf'] < tf_threshold].copy()
        
        # Optional secondary filter: binding energy
        if use_dg and 'dGtotal' in self.df.columns:
            problematic = problematic[problematic['dGtotal'] > dg_threshold].copy()
        
        # Sort by lowest Tf first (most problematic)
        problematic = problematic.sort_values('Tf', ascending=True)
        
        print(f"\nIdentified {len(problematic)} problematic strands:")
        print(f"  - Tf < {tf_threshold}: {len(problematic)}")
        if use_dg:
            print(f"  - Also filtered by dGtotal > {dg_threshold}: enabled")
        if len(problematic) > 0:
            print(f"  - Tf range: {problematic['Tf'].min():.1f} to {problematic['Tf'].max():.1f}")
            print(f"  - Length range: {problematic['Length'].min()} to {problematic['Length'].max()} bp")
        
        return problematic
    
    def generate_gc_optimized_sequence(self, original_seq, target_gc=0.5):
        """
        Generate an optimized sequence with better GC content.
        
        Strategy:
        - Replace weak AT pairs with GC pairs where possible
        - Maintain sequence diversity
        - Avoid long runs of same base
        
        Args:
            original_seq: Original DNA sequence
            target_gc: Target GC content (0.0 to 1.0)
            
        Returns:
            Optimized sequence string
        """
        seq = list(original_seq.upper())
        if len(seq) == 0:
            return original_seq
        
        current_gc = sum(1 for b in seq if b in 'GC') / len(seq)
        
        # If already at target, return original
        if abs(current_gc - target_gc) < 0.05:
            return original_seq
        
        # Need to increase GC content
        if current_gc < target_gc:
            # Replace A->G and T->C preferentially
            replacements = {'A': 'G', 'T': 'C'}
            needed = int((target_gc - current_gc) * len(seq))
            
            # Find positions to replace, avoiding runs
            positions = []
            for i, base in enumerate(seq):
                if base in replacements:
                    # Check neighbors to avoid creating long GC runs
                    prev_gc = (i > 0 and seq[i-1] in 'GC')
                    next_gc = (i < len(seq)-1 and seq[i+1] in 'GC')
                    
                    if not (prev_gc and next_gc):
                        positions.append(i)
            
            # Replace up to needed amount
            for i in positions[:needed]:
                seq[i] = replacements[seq[i]]
        
        # Need to decrease GC content
        else:
            replacements = {'G': 'A', 'C': 'T'}
            needed = int((current_gc - target_gc) * len(seq))
            
            positions = []
            for i, base in enumerate(seq):
                if base in replacements:
                    prev_at = (i > 0 and seq[i-1] in 'AT')
                    next_at = (i < len(seq)-1 and seq[i+1] in 'AT')
                    
                    if not (prev_at and next_at):
                        positions.append(i)
            
            for i in positions[:needed]:
                seq[i] = replacements[seq[i]]
        
        return ''.join(seq)
    
    def calculate_complement(self, seq):
        """Calculate the Watson-Crick complement of a sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(b, b) for b in seq.upper())
    
    def generate_overrides(self, problematic_strands, target_gc=0.55, max_overrides=100):
        """
        Generate sequence overrides for problematic strands.
        
        Args:
            problematic_strands: DataFrame of strands to optimize
            target_gc: Target GC content for optimized sequences
            max_overrides: Maximum number of overrides to generate
            
        Returns:
            List of override dicts with keys:
                HelixID, StartIdx, EndIdx, Direction,
                OriginalSeq, OptimizedSeq, OriginalTf,
                OriginalProbFold, OriginalGC, OptimizedGC, Length
        """
        overrides = []
        
        for idx, row in problematic_strands.head(max_overrides).iterrows():
            helix_id = row['HelixID']
            start_idx = min(row['StartIdx'], row['EndIdx'])
            end_idx = max(row['StartIdx'], row['EndIdx'])
            original_seq = row['Sequence']
            direction = row['Direction']
            
            # Generate optimized sequence
            optimized_seq = self.generate_gc_optimized_sequence(original_seq, target_gc)
            
            # If sequence changed, add to overrides
            if optimized_seq != original_seq:
                overrides.append({
                    'HelixID': int(helix_id),
                    'StartIdx': int(start_idx),
                    'EndIdx': int(end_idx),
                    'Direction': direction,
                    'OriginalSeq': original_seq,
                    'OptimizedSeq': optimized_seq,
                    'OriginalTf': float(row['Tf']),
                    'OriginalProbFold': float(row['ProbFold']),
                    'OriginalGC': sum(1 for b in original_seq if b in 'GC') / len(original_seq),
                    'OptimizedGC': sum(1 for b in optimized_seq if b in 'GC') / len(optimized_seq),
                    'Length': int(row['Length'])
                })
        
        print(f"\nGenerated {len(overrides)} sequence overrides")
        return overrides
    
    def _normalize_direction(self, direction_str):
        """
        Map whatever is in the 'Direction' column into 'fwd' or 'rev'
        (required by parse_override_input).
        """
        if not isinstance(direction_str, str):
            return 'fwd'
        d = direction_str.strip().lower()
        if d.startswith('f'):
            return 'fwd'
        if d.startswith('r'):
            return 'rev'
        # default to forward if unknown
        return 'fwd'
    
    def write_overrides_file(self, overrides):
        """
        Write overrides to:
        1) A machine-readable overrides file (self.output_path) in the format:

            HelixID[Idx][Direction]: OptimizedSequence

        e.g.  41[109][fwd]: GTT

        2) A separate annotation file '<basename>_annotated.txt' that keeps all
        the detailed comments, like:

            # Helix 41: 109-111 (fwd, 3 bp)
            # Original: Tf=24.3, ProbFold=0.202, GC=0.0%
            # Optimized: GC=33.3%
            # Original:  ATT
            # Optimized: GTT
            41[109][fwd]: GTT
        """
        # 1) Machine-readable overrides file for the mutation/apply script
        overrides_path = self.output_path
        # 2) Separate annotation file with all the verbose comments
        annotated_path = (
            overrides_path.parent
            / f"{overrides_path.stem}_annotated{overrides_path.suffix}"
        )

        # ---- Write minimal, parseable overrides file ----
        with open(overrides_path, 'w') as f:
            f.write("# Sequence overrides for DNA Origami optimization\n")
            f.write("# Generated by generate_sequence_overrides.py\n")
            f.write("#\n")
            f.write("# Format: HelixID[Idx][Direction]: OptimizedSequence\n")
            f.write("# Direction: fwd = forward strand, rev = reverse strand\n")
            f.write("#\n")
            f.write(f"# Total overrides: {len(overrides)}\n")
            f.write("#\n\n")

            for override in overrides:
                helix = override['HelixID']
                idx = override['StartIdx']
                direction = self._normalize_direction(override['Direction'])
                seq = override['OptimizedSeq']

                # This is the only line the apply script cares about
                f.write(f"{helix}[{idx}][{direction}]: {seq}\n")

        # ---- Write full annotated file ----
        with open(annotated_path, 'w') as f:
            f.write("# Detailed sequence overrides for DNA Origami optimization\n")
            f.write("# Generated by generate_sequence_overrides.py\n")
            f.write("# This file contains human-readable annotations.\n")
            f.write("# The last line of each block is in the same format as the\n")
            f.write("# machine-readable overrides file and can also be parsed.\n")
            f.write("#\n")
            f.write("# Example block:\n")
            f.write("# Helix 41: 109-111 (fwd, 3 bp)\n")
            f.write("# Original: Tf=24.3, ProbFold=0.202, GC=0.0%\n")
            f.write("# Optimized: GC=33.3%\n")
            f.write("# Original:  ATT\n")
            f.write("# Optimized: GTT\n")
            f.write("# 41[109][fwd]: GTT\n")
            f.write("#\n")
            f.write(f"# Total overrides: {len(overrides)}\n")
            f.write("#\n\n")

            for override in overrides:
                helix = override['HelixID']
                idx = override['StartIdx']
                direction = self._normalize_direction(override['Direction'])
                seq = override['OptimizedSeq']

                f.write(
                    f"# Helix {helix}: {override['StartIdx']}-{override['EndIdx']} "
                    f"({direction}, {override['Length']} bp)\n"
                )
                f.write(
                    f"# Original: Tf={override['OriginalTf']:.1f}, "
                    f"ProbFold={override['OriginalProbFold']:.3f}, "
                    f"GC={override['OriginalGC']*100:.1f}%\n"
                )
                f.write(
                    f"# Optimized: GC={override['OptimizedGC']*100:.1f}%\n"
                )
                f.write(f"# Original:  {override['OriginalSeq']}\n")
                f.write(f"# Optimized: {seq}\n")
                # Also include the same override line here for convenience
                f.write(f"{helix}[{idx}][{direction}]: {seq}\n\n")

        print(f"\nWrote overrides to: {overrides_path}")
        print(f"Wrote annotated overrides to: {annotated_path}")

    def generate_summary_report(self, overrides):
        """Generate a summary report of the optimizations (printed to stdout)."""
        if not overrides:
            print("\nNo overrides generated.")
            return
        
        print("\n" + "="*80)
        print("OPTIMIZATION SUMMARY")
        print("="*80)
        
        avg_original_tf = np.mean([o['OriginalTf'] for o in overrides])
        avg_original_prob = np.mean([o['OriginalProbFold'] for o in overrides])
        avg_original_gc = np.mean([o['OriginalGC'] for o in overrides])
        avg_optimized_gc = np.mean([o['OptimizedGC'] for o in overrides])
        
        print(f"\nOriginal strand statistics (average):")
        print(f"  Tf:       {avg_original_tf:.1f}")
        print(f"  ProbFold: {avg_original_prob:.3f}")
        print(f"  GC%%:      {avg_original_gc*100:.1f}%%")
        
        print(f"\nOptimized strand statistics (average):")
        print(f"  GC%%:      {avg_optimized_gc*100:.1f}%%")
        
        print(f"\nGC content change: {(avg_optimized_gc - avg_original_gc)*100:+.1f}%%")
        
        # Show top 5 most improved (by lowest original Tf)
        print("\nTop 5 strands by lowest original Tf:")
        for i, override in enumerate(sorted(overrides, key=lambda x: x['OriginalTf'])[:5], 1):
            print(
                f"  {i}. Helix {override['HelixID']} "
                f"[{override['StartIdx']}-{override['EndIdx']}]: "
                f"Tf={override['OriginalTf']:.1f}, "
                f"ProbFold={override['OriginalProbFold']:.3f}"
            )


def main():
    parser = argparse.ArgumentParser(
        description='Generate sequence overrides to optimize problematic DNA strands',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - optimize strands with Tf < 45°C (default)
  python generate_sequence_overrides.py data/design_strand_data.xlsx
  
  # Stricter Tf threshold
  python generate_sequence_overrides.py data/design_strand_data.xlsx --tf-threshold 50
  
  # Also consider binding energy
  python generate_sequence_overrides.py data/design_strand_data.xlsx --use-dg --dg-threshold -6.0
  
  # Specify output file and limit overrides
  python generate_sequence_overrides.py data/design_strand_data.xlsx -o custom_overrides.txt --max-overrides 50
        """
    )
    
    parser.add_argument('strand_data', type=str,
                        help='Path to strand_data.xlsx file from autobreak analysis')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output path for overrides file (default: auto-generate)')
    parser.add_argument('--tf-threshold', type=float, default=45.0,
                        help='Minimum acceptable Tf in Celsius (default: 45.0)')
    parser.add_argument('--use-dg', action='store_true',
                        help='Also filter by dG threshold (default: off, use only Tf)')
    parser.add_argument('--dg-threshold', type=float, default=-5.0,
                        help='Maximum acceptable dG in kcal/mol (only used with --use-dg)')
    parser.add_argument('--target-gc', type=float, default=0.55,
                        help='Target GC content for optimized sequences (default: 0.55)')
    parser.add_argument('--max-overrides', type=int, default=100,
                        help='Maximum number of overrides to generate (default: 100)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.strand_data).exists():
        print(f"Error: Input file not found: {args.strand_data}")
        sys.exit(1)
    
    if not 0.0 <= args.target_gc <= 1.0:
        print("Error: target-gc must be between 0.0 and 1.0")
        sys.exit(1)
    
    # Run optimization
    print("="*80)
    print("DNA ORIGAMI SEQUENCE OPTIMIZATION")
    print("="*80)
    
    optimizer = SequenceOptimizer(args.strand_data, args.output)
    
    # Identify problematic strands
    problematic = optimizer.identify_problematic_strands(
        tf_threshold=args.tf_threshold,
        use_dg=args.use_dg,
        dg_threshold=args.dg_threshold
    )
    
    if len(problematic) == 0:
        print("\nNo problematic strands found with current thresholds.")
        print("All strands meet the specified criteria.")
        return
    
    # Generate overrides
    overrides = optimizer.generate_overrides(
        problematic,
        target_gc=args.target_gc,
        max_overrides=args.max_overrides
    )
    
    if len(overrides) == 0:
        print("\nNo sequence changes needed - all sequences already optimal.")
        return
    
    # Write output file (in vh_num[idx][direction]: SEQ format)
    optimizer.write_overrides_file(overrides)
    
    # Generate summary
    optimizer.generate_summary_report(overrides)
    
    print("\n" + "="*80)
    print("NEXT STEPS:")
    print("="*80)
    print(f"1. Review the overrides in: {optimizer.output_path}")
    print("2. Apply overrides using your mutation/apply-overrides script")
    print("3. Re-run autobreak analysis to verify improvements")
    print("="*80)


if __name__ == "__main__":
    main()