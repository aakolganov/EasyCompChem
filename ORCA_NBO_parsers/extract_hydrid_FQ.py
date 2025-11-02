"""
This module processes vibrational data from computational chemistry output (.out) files,
reads hydride data from a CSV file, and computes vibrational frequency information
for specific atomic modes. The processed data is then saved into a new CSV file
containing best-mode frequencies for hydrides.

The module includes functions for loading vibrational data, parsing hydride fields,
and analyzing atom-specific vibrational modes for their frequency within pre-defined ranges.

Constants:
    MIN_FREQ (float): Minimum vibrational frequency considered IR-active (in cm^-1).
    MAX_FREQ (float): Maximum vibrational frequency considered IR-active (in cm^-1).

Functions:
    load_vibrational_data(filepath):
        Extracts vibrational frequencies and mode displacement matrix from the .out file.

    parse_hydride_field(field):
        Parses a formatted string of hydride indices and identifies separators.

    extract_best_mode_for_atom(matrix, freqs, atom_number):
        Determines the best vibrational mode for a specific atom based on amplitude
        and falls within the defined IR-active frequency range.

    process_row(row, hydride_cols, out_folder):
        Processes each row of the CSV, matches it to an .out file, computes vibrational
        frequencies for hydrides groups, and updates the row with new values.

Script Execution:
    Reads inputs from input_csv, processes the vibrational data using the specified folder
    out_folder, and saves the processed results to output_csv.
"""

## TODO - make the script read hydrides without input file also

import csv
import numpy as np
import sys
import os
import re

def load_vibrational_data(filepath):
    """
    Open the .out file and extract:
      - Vibrational frequencies (numpy array)
      - Full displacement (normal mode) matrix (numpy array)
    """
    freqs = []
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        sys.exit(f"Error opening {filepath}: {e}")

    # --- Extract frequencies ---
    in_freq_block = False
    for line in lines:
        if "VIBRATIONAL FREQUENCIES" in line:
            in_freq_block = True
            continue
        if in_freq_block:
            if not line.strip() or set(line.strip()) == set("-"):
                continue
            m = re.match(r"\s*(\d+):\s*([\d\.\-Ee]+)", line)
            if m:
                freqs.append(float(m.group(2)))
            else:
                if freqs:
                    in_freq_block = False

    # --- Parse displacement matrix blocks ---
    i, n = 0, len(lines)
    while i < n and "NORMAL MODES" not in lines[i]:
        i += 1
    if i == n:
        sys.exit(f"Could not find 'NORMAL MODES' section in {filepath}")
    i += 1

    # Skip to the first data header
    while i < n and not lines[i].lstrip().startswith(tuple('0123456789')):
        i += 1

    blocks = []
    block_nrows = None
    while i < n and lines[i].lstrip()[0].isdigit():
        # Skip header line
        i += 1
        data_lines = []

        # Collect data lines
        while i < n and lines[i].strip():
            tokens = lines[i].strip().split()
            # New header if all tokens are ints
            if data_lines and all(tok.isdigit() for tok in tokens):
                break
            data_lines.append(lines[i].rstrip())
            i += 1

        # Parse block
        block_data = []
        for dl in data_lines:
            toks = dl.split()
            if len(toks) < 2:
                continue
            block_data.append([float(x) for x in toks[1:]])

        if block_data:
            arr = np.array(block_data)
            if block_nrows is None:
                block_nrows = arr.shape[0]
            elif arr.shape[0] != block_nrows:
                sys.exit(f"Inconsistent row counts in blocks in {filepath}")
            blocks.append(arr)

        # Skip blank lines
        while i < n and not lines[i].strip():
            i += 1

    if not blocks:
        sys.exit(f"No displacement matrix data found in {filepath}")

    matrix = np.hstack(blocks)
    return matrix, np.array(freqs)

def parse_hydride_field(field):
    """
    Parse a hydride field preserving original grouping separator.
    Returns a list of index groups and the original separator.
    """
    f = field.strip()
    if not f or f.lower() == 'n/a':
        return None, None

    # Determine group separator (; or ,)
    if ';' in f:
        group_sep = ';'
    else:
        group_sep = ','

    group_strs = [grp.strip() for grp in f.split(group_sep) if grp.strip()]
    groups = []
    for grp in group_strs:
        parts = [p.strip() for p in grp.split(',') if p.strip()]
        groups.append([int(x) for x in parts])

    return groups, group_sep

MIN_FREQ = 800.0  # cm^-1
MAX_FREQ = 2700.0  # cm^-1

def extract_best_mode_for_atom(matrix, freqs, atom_number):
    """
    For a given atom (1-based), find the mode with maximal displacement within the
    IR-active range. Iterate down-ranked modes until one falls within [MIN_FREQ, MAX_FREQ].
    If none qualify, return the top mode regardless of frequency.
    """
    start_row = (atom_number - 1) * 3
    block = matrix[start_row:start_row+3, :]
    amps = np.sum(block**2, axis=0)

    # Sort mode indices by descending amplitude
    sorted_idx = np.argsort(amps)[::-1]

    # Look for first mode within desired range
    for idx in sorted_idx:
        freq_val = freqs[idx]
        if MIN_FREQ <= freq_val <= MAX_FREQ:
            return freq_val

    # Fallback: return highest-amplitude mode
    best_idx = sorted_idx[0]
    return freqs[best_idx]

def process_row(row, hydride_cols, out_folder):
    """
    Process one CSV row: locate its .out file, extract best-mode frequencies.
    Returns a dict mapping each hydride column to its new value.
    """
    # Determine base name
    fname = (row.get('Filename') or row.get('Basename') or '').strip()
    if not fname.lower().endswith('.xyz'):
        fname += '.xyz'
    print (fname)
    base = os.path.splitext(fname)[0]
    candidates = [f for f in os.listdir(out_folder)
                  if f.startswith(base) and f.endswith('.out')]
    if candidates:
        matrix, freqs = load_vibrational_data(os.path.join(out_folder, candidates[0]))
        new_vals = {}

        for col in hydride_cols:
            parsed, sep = parse_hydride_field(row[col])
            if parsed is None:
                new_vals[col] = row[col]
            else:
                group_outputs = []
                for group in parsed:
                    freqs_list = [f"{extract_best_mode_for_atom(matrix, freqs, atom):.2f}"
                                  for atom in group]
                    group_outputs.append(','.join(freqs_list))
                new_vals[col] = sep.join(group_outputs)

        return new_vals
    else:
        #sys.exit(f"Could not find .out file for {fname}")
        return('n/a')

if __name__ == '__main__':

    input_csv   = 'EFG_data_indeces_only.csv'
    out_folder  = 'data_FQ'
    output_csv  = 'Vib_FQ_data.csv'
    hydride_cols = [
        'Hydrides - mu1H',
        'Hydrides - mu2H',
        'Hydrides - mu3H'
    ]

    # Read CSV
    with open(input_csv, newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))

    # Process each row
    for row in rows:
        new_values = process_row(row, hydride_cols, out_folder)
        if new_values == 'n/a':
                row = 'n/a'
        else:
            for col, val in new_values.items():
                row[col] = val

    # Write output
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Extraction complete. Output written to {output_csv}")