"""
Module for extracting hybridization and polarization data from NBO analysis in
output files, focusing on specific atom pairs like P–O and Donor–O.

The module reads NBO output files, identifies relevant orbital data blocks, and
extracts key metrics for P–O and Donor–O type bonds based on provided donor atom
information. The results are saved into a CSV file summarizing the bond properties.

Additional functions and regular expressions are used to process the file content
and aggregate required data attributes.

Functions:
    - safe_float: Converts a string to a float safely, returning None on failure.
    - extract_single_nlmo_row: Extracts hybridization/polarization data for
      specific bond types (P–O or Donor–O) from an output file.
"""

import pandas as pd
import re

# Helper functions
def safe_float(s):
    try:
        return float(s)
    except:
        return None

# Relaxed regex for orbital headers
bd_line_pattern = re.compile(
    r"\(\s*(?P<occupancy>\d+\.\d+)\)\s+(?P<bd_percent>\d+\.\d+)% BD\s+\( 1\)\s*(?P<atom1>\w+)\s+(?P<idx1>\d+)-\s+(?P<atom2>\w+)\s+(?P<idx2>\d+)"
)

# Function to extract P–O and Donor–O data per file
def extract_single_nlmo_row(file_path, filename, donor_tag):
    donor_element, donor_index = donor_tag.split()
    result_row = {"Filename": filename}

    with open(file_path, "r") as f:
        content = f.read()

    section_match = re.search(
        r"Hybridization/Polarization Analysis of NLMOs in NAO Basis:\s+NLMO / Occupancy / Percent from Parent NBO / Atomic Hybrid Contributions(.*?)\n\n",
        content, re.DOTALL
    )
    if not section_match:
        return result_row

    nlmo_block = section_match.group(1)

    blocks = re.split(
        r"\n\s*(\d+\.\s*\(2\.00000\).*?BD\s+\( 1\)\s*\w+\s+\d+-\s+\w+\s+\d+)",
        nlmo_block
    )
    entries = []
    for i in range(1, len(blocks), 2):
        header = blocks[i]
        body = blocks[i + 1] if i + 1 < len(blocks) else ''
        entries.append(header + '\n' + body)

    for entry in entries:
        bd_header_match = bd_line_pattern.search(entry)
        if not bd_header_match:
            continue

        atom1, atom2 = bd_header_match.group("atom1"), bd_header_match.group("atom2")
        idx1, idx2 = bd_header_match.group("idx1"), bd_header_match.group("idx2")

        atom_pairs = {(atom1, idx1), (atom2, idx2)}
        bond_type = None
        if {("P", idx1), ("O", idx2)} <= atom_pairs or {("P", idx2), ("O", idx1)} <= atom_pairs:
            bond_type = "P–O"
        elif {(donor_element, donor_index), ("O", idx1)} <= atom_pairs or {(donor_element, donor_index), ("O", idx2)} <= atom_pairs:
            bond_type = "Donor–O"
        else:
            continue

        occupancy = safe_float(bd_header_match.group("occupancy"))
        bd_percent = safe_float(bd_header_match.group("bd_percent"))

        atom_matches = re.findall(
            r"(\d+\.\d+)%\s+(\w+)\s+(\d+)\s+s\(\s*(\d+\.\d+)%\)p\s*([\d\.]+)\(\s*(\d+\.\d+)%\)d\s*[\d\.]+\(\s*(\d+\.\d+)%\)",
            entry
        )

        donor_contrib = o_contrib = p_contrib = None
        donor_ps = donor_d_percent = None
        o_ps = o_d_percent = None
        p_ps = p_d_percent = None

        for match in atom_matches:
            percent, elem, idx, s_pct, p_idx, p_pct, d_pct = match
            percent = safe_float(percent)
            p_idx = safe_float(p_idx)
            d_pct = safe_float(d_pct)

            if elem == donor_element and idx == donor_index:
                donor_contrib = percent
                donor_ps = p_idx
                donor_d_percent = d_pct
            elif elem == "O":
                o_contrib = percent
                o_ps = p_idx
                o_d_percent = d_pct
            elif elem == "P":
                p_contrib = percent
                p_ps = p_idx
                p_d_percent = d_pct

        if bond_type == "P–O" and None not in (p_contrib, o_contrib):
            result_row.update({
                "PO_Occupancy": occupancy,
                "PO_BD%": bd_percent,
                "PO_O_contribution_%": o_contrib,
                "PO_P_contribution_%": p_contrib,
                "PO_O:P_ratio": o_contrib / p_contrib if p_contrib else None,
                "PO_O_p/s": o_ps,
                "PO_O_d%": o_d_percent,
                "PO_P_p/s": p_ps,
                "PO_P_d%": p_d_percent
            })

        elif bond_type == "Donor–O" and None not in (donor_contrib, o_contrib):
            result_row.update({
                "DO_Occupancy": occupancy,
                "DO_BD%": bd_percent,
                "DO_O_contribution_%": o_contrib,
                "DO_Donor_contribution_%": donor_contrib,
                "DO_Donor:O_ratio": donor_contrib / o_contrib if o_contrib else None,
                "DO_O_p/s": o_ps,
                "DO_O_d%": o_d_percent,
                "DO_Donor_p/s": donor_ps,
                "DO_Donor_d%": donor_d_percent
            })

    return result_row

if __name__ == "__main__":
    donor_df = pd.read_csv("Donor_Data_Rel_NLMO.csv")
    donor_map = dict(zip(donor_df["Filename"], donor_df["Acceptor_atom"]))

    # Run extraction on all .out files in directory
    all_rows = []
    from pathlib import Path

    folder = Path("data_HLG/CMO_PBE0_NBO")
    for path in folder.glob("*.out"):
        filename = path.name
        if filename not in donor_map:
            continue
        donor_tag = donor_map[filename]
        row = extract_single_nlmo_row(path, filename, donor_tag)
        all_rows.append(row)

    # Save to CSV
    df = pd.DataFrame(all_rows)
    df.to_csv("nlmo_relativistic_bond_summary.csv", index=False)
    df