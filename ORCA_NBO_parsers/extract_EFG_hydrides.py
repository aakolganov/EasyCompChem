"""
Contains utility functions for processing ORCA output files and calculating
properties related to electric field gradients (EFG) of hydride nuclei.

Provides functionality to extract (e**2qQ, eta) values, raw EFG matrices, and compute
canonical orientation vectors.

Functions:
- extract_eqQ_and_eta_values: Reads an ORCA output file and extracts quadrupole coupling
  constants (e**2qQ) and asymmetry parameters (eta) for hydride nuclei.
- process_hydride_csv: Processes a CSV file containing hydride indices and computes the
  corresponding EFG values, then writes the results to a new CSV file.
- get_canonical_orientation: Diagonalizes an EFG matrix and computes the canonical
  orientation in a right-handed coordinate system.
- extract_raw_efg_matrices: Extracts the raw EFG matrices from an ORCA output file for
  all hydride nuclei.
"""

## TODO - make the script read hydrides without input file also

import os
import pandas as pd
import re
import numpy as np

def extract_eqQ_and_eta_values(file_path):
    """
    Extracts e**2qQ and eta values from an ORCA output file for each hydride nucleus.

    Args:
        file_path (str): Path to the ORCA output file.

    Returns:
        dict: A dictionary with hydride indices and their corresponding (e**2qQ, eta) values.
    """
    with open(file_path, 'r') as file:
        content = file.readlines()

    efg_values = {}
    current_nucleus = None
    current_eqQ = None

    for line in content:
        # Match "Nucleus [Index]H"
        nucleus_match = re.search(r'Nucleus\s+(\d+)H', line)
        if nucleus_match:
            current_nucleus = int(nucleus_match.group(1))  # Extract nucleus index
            current_eqQ = None  # Reset for new nucleus

        # Match "e**2qQ = [value] MHz"
        if current_nucleus is not None:
            eqQ_match = re.search(r'e\*\*2qQ\s+=\s+([-\d.]+)\s+MHz', line)
            if eqQ_match:
                current_eqQ = np.abs(float(eqQ_match.group(1)))

            # Match "eta = [value]"
            eta_match = re.search(r'eta\s+=\s+([-\d.]+)', line)
            if eta_match and current_eqQ is not None:
                eta_value = float(eta_match.group(1))
                efg_values[current_nucleus] = (current_eqQ, eta_value)
                current_nucleus = None  # Reset for the next block
                current_eqQ = None

    return efg_values

def process_hydride_csv(input_csv_path,
                        output_csv_path,
                        orca_folder_path,
                        include_eta=True):
    """
    Processes a CSV of hydride indices, extracts e**2qQ (and optionally eta),
    and writes to a new CSV.

    Args:
      include_eta (bool): if False, only e**2qQ is output for all hydride types.
    """
    df_in  = pd.read_csv(input_csv_path)
    df_out = df_in.copy()

    for idx, row in df_in.iterrows():
        fname = row['Filename']
        if pd.isna(fname):
            continue
        orca_file = os.path.join(orca_folder_path,
                                 fname.replace('.xyz','_input.inp.out'))
        if not os.path.isfile(orca_file):
            continue

        efg_map = extract_eqQ_and_eta_values(orca_file)

        for col in ['Hydrides - mu1H','Hydrides - mu2H','Hydrides - mu3H','H2_coord']:
            raw = row[col]
            if pd.isna(raw):
                continue

            sep = ';' if ';' in raw else ','
            entries = [e.strip() for e in raw.split(sep)]
            out_entries = []

            for ent in entries:
                idxs = [int(i) for i in ent.split(',')]
                if include_eta and col in ['Hydrides - mu2H','Hydrides - mu3H']:
                    # (e2qQ, eta) per index
                    pieces = []
                    for i in idxs:
                        val = efg_map.get(i-1)
                        if val:
                            pieces.append(f"({val[0]:.6f},{val[1]:.6f})")
                        else:
                            pieces.append("(N/A,N/A)")
                    out_entries.append(','.join(pieces))
                else:
                    # only e2qQ
                    pieces = []
                    for i in idxs:
                        val = efg_map.get(i-1)
                        if val:
                            pieces.append(f"{val[0]:.6f}")
                        else:
                            pieces.append("N/A")
                    out_entries.append(','.join(pieces))

            df_out.at[idx, col] = sep.join(out_entries)

    df_out.to_csv(output_csv_path, index=False)


# Example usage:
process_hydride_csv(
    input_csv_path  = "EFG_data_indeces_only.csv",
    output_csv_path = "EFG_data.csv",
    orca_folder_path= "data_EFG",
    include_eta     = False   # set True to include (e2qQ, eta) for mu2H/mu3H
)

import os
import pandas as pd
import re
import numpy as np

def get_canonical_orientation(raw_efg_matrix):
    """
    Diagonalizes a raw EFG matrix and returns the canonical orientation vectors.

    The canonical orientation is defined by sorting the eigenvectors such that the
    eigenvalue with the largest absolute value corresponds to the Z-axis. The
    resulting orientation forms a right-handed coordinate system.

    Args:
        raw_efg_matrix (np.ndarray): A 3x3 numpy array of the raw EFG tensor.

    Returns:
        dict: A dictionary with the new canonical orientation vectors.
              Format: {'X': [x1,x2,x3], 'Y': [y1,y2,y3], 'Z': [z1,z2,z3]}
    """
    if raw_efg_matrix is None:
        return {}

    # Step 1: Perform eigendecomposition of the symmetric EFG matrix
    eigvals, eigvecs = np.linalg.eigh(raw_efg_matrix)

    # Step 2: Find the index of the eigenvalue with the maximum absolute value
    max_abs_idx = np.argmax(np.abs(eigvals))

    # Step 3: Define the new Z-axis from the corresponding eigenvector
    z_axis = eigvecs[:, max_abs_idx]

    # Step 4: Identify the other two eigenvectors
    other_indices = [i for i in range(3) if i != max_abs_idx]
    v1 = eigvecs[:, other_indices[0]]
    v2 = eigvecs[:, other_indices[1]]

    # Step 5: Create a right-handed coordinate system
    # We assign one of the remaining vectors to the X-axis and define the Y-axis
    # using the cross product to ensure a right-handed system (X x Y = Z).
    x_axis = v1
    y_axis = np.cross(z_axis, x_axis)

    # The calculated y_axis should be parallel or anti-parallel to v2.
    # We check the dot product to ensure consistency and flip if necessary.
    if np.dot(y_axis, v2) < 0:
        y_axis = -y_axis

    # The new orientation matrix R has [x_axis, y_axis, z_axis] as columns.
    # We need to ensure it's a proper rotation (determinant = +1).
    # Since we constructed it with a cross product, it should be right-handed.

    return {
        'X': x_axis.tolist(),
        'Y': y_axis.tolist(),
        'Z': z_axis.tolist()
    }


def extract_raw_efg_matrices(file_path):
    """
    Extracts Raw EFG matrices from an ORCA output file for each hydride nucleus.

    Args:
        file_path (str): Path to the ORCA output file.

    Returns:
        dict: A dictionary with hydride indices and their corresponding raw EFG matrix.
              Format: {nucleus_index: np.ndarray}
    """
    with open(file_path, 'r') as file:
        content = file.readlines()

    raw_efg_matrices = {}
    current_nucleus = None
    efg_section = False
    efg_matrix_lines = []

    for line in content:
        # Match "Nucleus [Index]H" to identify the start of a new atom's data
        nucleus_match = re.search(r'Nucleus\s+(\d+)H', line)
        if nucleus_match:
            # If we were in an EFG section, something went wrong, so we reset
            if efg_section:
                efg_section = False
                efg_matrix_lines = []
            current_nucleus = int(nucleus_match.group(1))

        # Check if we're entering the Raw EFG matrix section
        if current_nucleus is not None and 'Raw EFG matrix' in line:
            efg_section = True
            efg_matrix_lines = []
            continue

        # If in the EFG section, collect the matrix lines
        if efg_section and len(efg_matrix_lines) < 3:
            # Regex to capture the three float values in a line
            values_match = re.findall(r'([-]?\d+\.\d+)', line)
            if len(values_match) == 3:
                efg_matrix_lines.append([float(v) for v in values_match])

            # Once we have 3 lines, the matrix is complete
            if len(efg_matrix_lines) == 3:
                raw_efg_matrices[current_nucleus] = np.array(efg_matrix_lines)
                current_nucleus = None
                efg_section = False
                efg_matrix_lines = []

    return raw_efg_matrices


def format_orientation_output(orientation_dict):
    """
    Formats orientation vectors into a compact string.

    Args:
        orientation_dict (dict): Dictionary containing X, Y, Z orientation vectors.

    Returns:
        str: Formatted orientation string "X(...);Y(...);Z(...)". Returns 'N/A' if dict is invalid.
    """
    if not orientation_dict or any(key not in orientation_dict for key in ['X', 'Y', 'Z']):
        return 'N/A'

    x_str = f"X({','.join(f'{val:.6f}' for val in orientation_dict['X'])})"
    y_str = f"Y({','.join(f'{val:.6f}' for val in orientation_dict['Y'])})"
    z_str = f"Z({','.join(f'{val:.6f}' for val in orientation_dict['Z'])})"
    return f"{x_str};{y_str};{z_str}"


def process_hydride_orientations_csv(input_csv_path, output_csv_path, orca_folder_path):
    """
    Processes a CSV, extracts raw EFG matrices, calculates canonical orientations,
    and saves them to a new CSV file.

    Args:
        input_csv_path (str): Path to the input CSV file.
        output_csv_path (str): Path to save the output CSV file.
        orca_folder_path (str): Path to the folder containing ORCA output files.
    """
    input_data = pd.read_csv(input_csv_path)
    output_data = input_data.copy()

    # Create new columns for the orientation data to avoid issues
    hydride_cols = [col for col in input_data.columns if 'Hydrides' in col]
    for col in hydride_cols:
         output_data[f"{col}_orientations"] = ''

    for idx, row in input_data.iterrows():
        filename = row['Filename']
        if pd.isna(filename):
            continue

        orca_file_path = os.path.join(orca_folder_path, filename.replace('.xyz', '_input.inp.out'))
        if not os.path.isfile(orca_file_path):
            print(f"Warning: ORCA output file not found for {filename}, skipping.")
            continue

        # Extract all raw EFG matrices from the ORCA file at once
        raw_efg_matrices = extract_raw_efg_matrices(orca_file_path)

        for col in hydride_cols:
            if pd.isna(row[col]):
                continue

            entries = str(row[col]).split(';') if ';' in str(row[col]) else str(row[col]).split(',')
            separator = ';' if ';' in str(row[col]) else ','

            orientation_results = []
            for entry in entries:
                if not entry: continue
                indices = list(map(int, entry.split(',')))

                entry_orientations = []
                for hydride_idx in indices:
                    # Get the raw EFG matrix for the current hydride
                    raw_matrix = raw_efg_matrices.get(hydride_idx)

                    # Calculate the canonical orientation from the raw matrix
                    orientation_dict = get_canonical_orientation(raw_matrix)

                    # Format the result for output
                    entry_orientations.append(format_orientation_output(orientation_dict))

                orientation_results.append(','.join(entry_orientations))

            # Store the formatted string in the correct new column
            new_col_name = f"{col}_orientations"
            output_data.loc[idx, new_col_name] = separator.join(orientation_results)

    # Save the output to a new CSV file
    output_data.to_csv(output_csv_path, index=False)
    print(f"Canonical EFG orientations extracted and saved to {output_csv_path}")

# Example usage:
if __name__ == "__main__":
    # --- PLEASE CONFIGURE YOUR PATHS HERE ---
    input_csv_path = "EFG_data_indeces_only.csv"
    output_csv_path = 'EFG_orientations_canonical.csv'
    orca_folder_path = 'data_EFG'

    process_hydride_orientations_csv(input_csv_path, output_csv_path, orca_folder_path)

