"""
generate_sequence_overrides.py

Generate sequence overrides to optimize problematic DNA strands based on
thermodynamic analysis, write them in a format compatible with the
override-application script, and optionally apply the overrides to
produce a full corrected scaffold sequence.

Outputs:

1) overrides file (machine-readable):
   HelixID[Idx][Direction]: OptimizedSequence
   e.g.  41[109][fwd]: GTT

2) annotated overrides file (human-readable + same override lines)

3) OPTIONAL: full corrected scaffold sequence .txt if cadnano design,
   original sequence, and output path are provided.
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import sys
import re

import cadnano
from cadnano.document import Document


# ----------------------------------------------------------------------
# Helper functions for applying overrides (same logic as first script)
# ----------------------------------------------------------------------

def parse_override_input(file_name):
    """Parses a text file containing strand data into a dictionary.

    Expected line format (one per override):

        vh_num[idx][direction]: SEQUENCE

    where:
        - vh_num: integer helix / virtual helix ID
        - idx: integer base index within that helix
        - direction: 'fwd' or 'rev'
        - SEQUENCE: DNA sequence (A/T/C/G, spaces allowed)

    Comment lines starting with '#' are ignored.

    Returns:
        dict: keys = (vh_num, idx, is_fwd), values = sequence (no spaces).
    """
    override_data = {}
    with open(file_name, 'r') as file:
        for line in file:
            # Example: 12[189][fwd]: AGATC GACGCT  GAA AGCGTC GATCT
            match = re.match(
                r'(\d+)\[(\d+)\]\[(fwd|rev)\]:[\s\u2014]+([ATCG\s]+)',
                line
            )
            if match:
                vh_num, idx, direction, sequence = match.groups()
                is_fwd = True if direction == 'fwd' else False
                override_data[(int(vh_num), int(idx), is_fwd)] = (
                    sequence.replace(" ", "").rstrip()
                )
    return override_data


def get_base_length_from_start(part, vh_num, is_forward, idx):
    """
    Given a cadnano Part, find the specified base (vh_num, strand direction, idx)
    and calculate its distance from the start of its oligo (1-based).

    Returns:
        int or None: 1-based position of this base along the scaffold/oligo,
                     or None if no strand is found.
    """
    strand = part.getStrand(is_forward, vh_num, idx)
    if strand is None:
        print(f"No strand found at ({vh_num}, {is_forward}, {idx})")
        return None

    length = 0
    for s in strand.oligo().strand5p().generator3pStrand():
        if s == strand:
            length += idx - s.idx5Prime() + 1
            break
        length += s.totalLength()
    return length


def apply_sequence_overrides(input_cadnano, input_sequence, input_overrides, output_sequence):
    """
    Applies sequence overrides to a scaffold sequence using a cadnano design.

    Args:
        input_cadnano (str): Path to the input cadnano JSON file.
        input_sequence (str): Path to the original full scaffold sequence .txt.
        input_overrides (str): Path to overrides file in vh_num[idx][direction]: SEQ format.
        output_sequence (str): Path to write the modified full scaffold sequence (.txt).
    """
    # Load cadnano design
    app = cadnano.app()
    doc = app.document = Document()
    doc.readFile(input_cadnano)
    part = doc.activePart()

    # Load original sequence (full scaffold)
    with open(input_sequence, "r") as f:
        original_sequence = f.read().strip()

    # Parse overrides
    overrides = parse_override_input(input_overrides)

    # Apply overrides to the original sequence
    modified_sequence = list(original_sequence)
    for (vh_num, idx, is_forward), override_seq in overrides.items():
        start_idx = get_base_length_from_start(part, vh_num, is_forward, idx)
        if start_idx is None:
            continue

        start_idx0 = start_idx - 1  # convert to 0-based index
        end_idx0 = start_idx0 + len(override_seq)
        modified_sequence[start_idx0:end_idx0] = override_seq

    # Save modified sequence
    with open(output_sequence, "w") as f:
        f.write("".join(modified_sequence))


# ----------------------------------------------------------------------
# Optimizer class (unchanged logic, new write_overrides_file behavior)
# ----------------------------------------------------------------------

class SequenceOptimizer:
    """
    Analyzes strand thermodynamic data and generates sequence overrides
    to improve problematic regions.
    """
    def __init__(self, strand_data_path, output_path=None):
        self.strand_data_path = Path(strand_data_path)

        if output_path:
            self.output_path = Path(output_path)
        else:
            base_name = self.strand_data_path.stem.replace('_strand_data', '')
            self.output_path = self.strand_data_path.parent / f'{base_name}_overrides.txt'

        # Load data
        self.df = pd.read_excel(strand_data_path, sheet_name='StrandData')

        print(f"Loaded {len(self.df)} strands from {strand_data_path}")

        numeric_cols = [
            'Tf', 'ProbFold', 'dGtotal', 'dGhyb', 'dGloop', 'dGconc',
            'Length', 'HelixID', 'StartIdx', 'EndIdx', 'maxTm'
        ]
        for col in numeric_cols:
            if col in self.df.columns:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')

        if 'Tf' in self.df.columns:
            print(f"  Tf range: {self.df['Tf'].min():.1f} to {self.df['Tf'].max():.1f}")
            print(f"  Tf mean: {self.df['Tf'].mean():.1f}")
            print(f"  Strands with Tf < 45: {len(self.df[self.df['Tf'] < 45.0])}")

    def identify_problematic_strands(self, tf_threshold=45.0, use_dg=False,
                                     dg_threshold=-5.0):
        problematic = self.df[self.df['Tf'] < tf_threshold].copy()

        if use_dg and 'dGtotal' in self.df.columns:
            problematic = problematic[problematic['dGtotal'] > dg_threshold].copy()

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
        seq = list(original_seq.upper())
        if len(seq) == 0:
            return original_seq

        current_gc = sum(1 for b in seq if b in 'GC') / len(seq)

        if abs(current_gc - target_gc) < 0.05:
            return original_seq

        if current_gc < target_gc:
            replacements = {'A': 'G', 'T': 'C'}
            needed = int((target_gc - current_gc) * len(seq))

            positions = []
            for i, base in enumerate(seq):
                if base in replacements:
                    prev_gc = (i > 0 and seq[i-1] in 'GC')
                    next_gc = (i < len(seq)-1 and seq[i+1] in 'GC')
                    if not (prev_gc and next_gc):
                        positions.append(i)

            for i in positions[:needed]:
                seq[i] = replacements[seq[i]]
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
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(b, b) for b in seq.upper())

    def generate_overrides(self, problematic_strands, target_gc=0.55, max_overrides=100):
        overrides = []

        for _, row in problematic_strands.head(max_overrides).iterrows():
            helix_id = row['HelixID']
            start_idx = min(row['StartIdx'], row['EndIdx'])
            end_idx = max(row['StartIdx'], row['EndIdx'])
            original_seq = row['Sequence']
            direction = row['Direction']

            optimized_seq = self.generate_gc_optimized_sequence(original_seq, target_gc)

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
        if not isinstance(direction_str, str):
            return 'fwd'
        d = direction_str.strip().lower()
        if d.startswith('f'):
            return 'fwd'
        if d.startswith('r'):
            return 'rev'
        return 'fwd'

    def write_overrides_file(self, overrides):
        """
        Write overrides to:
        1) machine-readable overrides file (self.output_path)
        2) annotated overrides file '<basename>_annotated.txt'
        """
        overrides_path = self.output_path
        annotated_path = (
            overrides_path.parent
            / f"{overrides_path.stem}_annotated{overrides_path.suffix}"
        )

        # --- minimal overrides file ---
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
                f.write(f"{helix}[{idx}][{direction}]: {seq}\n")

        # --- annotated overrides file ---
        with open(annotated_path, 'w') as f:
            f.write("# Detailed sequence overrides for DNA Origami optimization\n")
            f.write("# Generated by generate_sequence_overrides.py\n")
            f.write("# This file contains human-readable annotations.\n")
            f.write("# The last line of each block is in the same format as the\n")
            f.write("# machine-readable overrides file and can also be parsed.\n")
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
                f.write(f"{helix}[{idx}][{direction}]: {seq}\n\n")

        print(f"\nWrote overrides to: {overrides_path}")
        print(f"Wrote annotated overrides to: {annotated_path}")

    def generate_summary_report(self, overrides):
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

        print("\nTop 5 strands by lowest original Tf:")
        for i, override in enumerate(
            sorted(overrides, key=lambda x: x['OriginalTf'])[:5], 1
        ):
            print(
                f"  {i}. Helix {override['HelixID']} "
                f"[{override['StartIdx']}-{override['EndIdx']}]: "
                f"Tf={override['OriginalTf']:.1f}, "
                f"ProbFold={override['OriginalProbFold']:.3f}"
            )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate and optionally apply sequence overrides to a DNA origami scaffold',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - generate overrides (Tf < 45°C), no application
  python generate_sequence_overrides.py data/design_strand_data.xlsx

  # Generate overrides and apply them to produce a corrected scaffold
  python generate_sequence_overrides.py data/design_strand_data.xlsx \
      --cadnano dbp.q1.json \
      --sequence 2_scaffold_3364.txt \
      --output-sequence corrected_scaffold.txt
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

    # NEW: optional arguments to directly apply overrides
    parser.add_argument('--cadnano', type=str,
                        help='Path to cadnano JSON design file (optional; needed to apply overrides)')
    parser.add_argument('--sequence', type=str,
                        help='Path to original full scaffold sequence .txt (optional; needed to apply overrides)')
    parser.add_argument('--output-sequence', type=str,
                        help='Path to write full corrected scaffold sequence .txt (optional)')

    args = parser.parse_args()

    if not Path(args.strand_data).exists():
        print(f"Error: Input file not found: {args.strand_data}")
        sys.exit(1)

    if not 0.0 <= args.target_gc <= 1.0:
        print("Error: target-gc must be between 0.0 and 1.0")
        sys.exit(1)

    print("="*80)
    print("DNA ORIGAMI SEQUENCE OPTIMIZATION")
    print("="*80)

    optimizer = SequenceOptimizer(args.strand_data, args.output)

    problematic = optimizer.identify_problematic_strands(
        tf_threshold=args.tf_threshold,
        use_dg=args.use_dg,
        dg_threshold=args.dg_threshold
    )

    if len(problematic) == 0:
        print("\nNo problematic strands found with current thresholds.")
        print("All strands meet the specified criteria.")
        return

    overrides = optimizer.generate_overrides(
        problematic,
        target_gc=args.target_gc,
        max_overrides=args.max_overrides
    )

    if len(overrides) == 0:
        print("\nNo sequence changes needed - all sequences already optimal.")
        return

    # 1) write overrides + annotated version
    optimizer.write_overrides_file(overrides)

    # 2) summary report
    optimizer.generate_summary_report(overrides)

    # 3) optionally apply overrides to produce full corrected sequence
    if args.cadnano or args.sequence or args.output_sequence:
        if not (args.cadnano and args.sequence and args.output_sequence):
            print("\nWARNING: To apply overrides and produce a corrected scaffold,")
            print("you must provide --cadnano, --sequence, and --output-sequence together.")
        else:
            overrides_path = optimizer.output_path  # machine-readable overrides file
            for path_str, label in [
                (args.cadnano, "cadnano design"),
                (args.sequence, "original sequence"),
                (overrides_path, "overrides file"),
            ]:
                p = Path(path_str) if label != "overrides file" else overrides_path
                if not p.exists():
                    raise FileNotFoundError(f"{label} not found: {p}")

            apply_sequence_overrides(
                args.cadnano,
                args.sequence,
                str(overrides_path),
                args.output_sequence
            )
            print(f"\nFull corrected scaffold sequence written to: {args.output_sequence}")

    print("\n" + "="*80)
    print("DONE")
    print("="*80)


if __name__ == "__main__":
    main()