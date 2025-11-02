"""
Extraction of second-order perturbation energies from NBO output files,
optimized for sigma P-C* orbitals.

This module contains functions to parse second-order perturbation entries
from NBO output files and extract LP orbital indices for specific oxygen atoms.
The results are returned as pandas DataFrames, which can be used for further analysis or saved as CSV files.
"""

import re
import pandas as pd
from pathlib import Path


# Original regex for second-order perturbation entries
pattern = re.compile(
    r"""
    ^\s*\d+\.\s+LP\s*\(\s*(?P<lp_num>\d+)\)\s*O\s*(?P<o_idx>\d+)
    \s+\d+\.\s+BD\*\(\s*1\)\s*
    (?: # two alternative orderings:
    C\s*(?P<c_idx1>\d+)-\s*P\s*(?P<p_idx1>\d+)
    | P\s*(?P<p_idx2>\d+)-\s*C\s*(?P<c_idx2>\d+)
    )
    \s+(?P<E2>[0-9]+\.?[0-9]*)\s+(?P<ENL_EL>[0-9]+\.?[0-9]*)\s+(?P<F>[0-9]+\.?[0-9]*)
    """,
    re.VERBOSE
)

# Another regex to extract LP orbital indices for any oxygen
lp_orbital_pattern = re.compile(
    r"""
    ^\s*(?P<orbital_idx>\d+)\.\s+LP\s*\(\s*(?P<lp_num>\d+)\)\s*O\s*(?P<o_idx>\d+)
    """,
    re.VERBOSE
)

def parse_nbo_second_order(file_path):
    """
    Parse an NBO output file for LP(O) -> BD*(C-P) second-order perturbation entries.
    Returns a DataFrame with columns: ['file','lp_num','o_idx','c_idx','p_idx','E2']
    """
    records = []
    with open(file_path, 'r') as f:
        in_block = False
        for line in f:
            if 'SECOND ORDER PERTURBATION THEORY ANALYSIS' in line:
                in_block = True
                continue
            if in_block:
                if line.strip() == '' or line.startswith(' within unit'):
                    continue
                if line.strip().startswith('---') or line.strip().startswith('Total'):
                    break
                m = pattern.match(line)
                if m:
                    d = m.groupdict()
                    lp_num = int(d['lp_num'])
                    # prefer the first branch, otherwise fallback to the second
                    if d['c_idx1']:
                        c_idx = int(d['c_idx1']); p_idx = int(d['p_idx1'])
                    else:
                        c_idx = int(d['c_idx2']); p_idx = int(d['p_idx2'])
                    E2 = float(d['E2'])
                    records.append({
                        'Filename': Path(file_path).stem,
                        'lp_num': lp_num,
                        'c_idx': c_idx,
                        'E2': E2
                    })
    return pd.DataFrame(records)

def extract_lp_orbital_indices(file_path):
    """
    Extract LP orbital indices for the oxygen with the highest index.
    Returns a DataFrame with one row per file, with orbital indices listed comma-separated.
    """
    lp_records = []

    with open(file_path, 'r') as f:
        in_block = False
        for line in f:
            if 'SECOND ORDER PERTURBATION THEORY ANALYSIS' in line:
                in_block = True
                continue
            if in_block:
                if line.strip() == '' or line.startswith(' within unit'):
                    continue
                if line.strip().startswith('---') or line.strip().startswith('Total'):
                    break

                # Match any LP orbital (not just LP->BD* interactions)
                m = lp_orbital_pattern.match(line)
                if m:
                    d = m.groupdict()
                    lp_records.append({
                        'Filename': Path(file_path).stem,
                        'o_idx': int(d['o_idx']),
                        'lp_num': int(d['lp_num']),
                        'orbital_idx': int(d['orbital_idx'])
                    })

    if not lp_records:
        return pd.DataFrame()

    # Convert to DataFrame and find the highest oxygen index
    df = pd.DataFrame(lp_records)
    max_o_idx = df['o_idx'].max()

    # Filter for only the highest oxygen index
    last_o_lps = df[df['o_idx'] == max_o_idx].copy()

    # Group by filename only and concatenate ALL orbital indices with commas
    grouped = last_o_lps.groupby(['Filename', 'o_idx']).agg({
        'orbital_idx': lambda x: ','.join(map(str, sorted(x.unique())))
    }).reset_index()

    # Rename the column to be more descriptive
    grouped = grouped.rename(columns={'orbital_idx': 'orbital_indices'})

    return grouped

#example usage

if __name__ == '__main__':
    # 1. Get all .out files in the specified folder
    files = list(Path('data_HLG/CMO_PBE0_NBO').glob('*.out'))

    # 2. Parse each file for second-order perturbations
    dfs = [parse_nbo_second_order(f) for f in files]
    all_df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    # 3. Extract LP orbital indices for the last oxygen
    lp_orbital_dfs = [extract_lp_orbital_indices(f) for f in files]
    lp_orbital_df = pd.concat(lp_orbital_dfs, ignore_index=True) if lp_orbital_dfs else pd.DataFrame()

    # 4. Prepare summary per file (original functionality)
    summary_records = []
    lps = [1, 2, 3]

    if not all_df.empty:
        for fname, group in all_df.groupby('Filename'):
            c_list = sorted(group['c_idx'].unique())
            label_map = {c: f'C{i+1}' for i, c in enumerate(c_list)}
            rec = {'Filename': fname}

            for lp in lps:
                lp_group = group[group['lp_num'] == lp]
                if lp_group.empty:
                    for c in c_list:
                        rec[f'LP{lp}-{label_map[c]}'] = None
                    rec[f'LP{lp}-sum'] = None
                else:
                    for c in c_list:
                        val = lp_group.loc[lp_group['c_idx'] == c, 'E2']
                        rec[f'LP{lp}-{label_map[c]}'] = val.iloc[0] if not val.empty else None
                    rec[f'LP{lp}-sum'] = lp_group['E2'].sum()

            rec['Total-sum'] = group['E2'].sum() if not group.empty else None
            summary_records.append(rec)

    summary_df = pd.DataFrame(summary_records)

    #5. Write original summary to CSV
    summary_df.to_csv('nbo_LP_CP_summary.csv', index=False, na_rep='')

    # Write LP orbital indices to separate CSV
    if not lp_orbital_df.empty:
        lp_orbital_df.to_csv('lp_orbital_indices_last_oxygen.csv', index=False)
        print(f"LP orbital indices for last oxygen saved to 'lp_orbital_indices_last_oxygen.csv'")
        print(f"Found {len(lp_orbital_df)} LP entries for the highest oxygen indices")
    else:
        print("No LP orbital indices found")

    # Display results
    print("Summary DataFrame:")
    print(summary_df)

    if not lp_orbital_df.empty:
        print("\nLP Orbital Indices for Last Oxygen:")
        print(lp_orbital_df)


