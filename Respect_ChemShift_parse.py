"""
Module for extracting chemical shielding data from ReSpect.cs files and saving the results
to a CSV file.

This module defines functions to process .cs files in a specified folder, extract
specific shielding components (SIGMA_11, SIGMA_22, SIGMA_33), as well as isotropic
and directional components for three rows (DIA, PARA, SUM) from the NMR shielding
tensor table. The results are saved to a CSV file for further analysis.
"""

import re, glob, os
import pandas as pd


def extract_shielding_info(folder=".", output="shielding.csv"):
    """
    Extracts chemical shielding data from .cs files in the specified folder.

    For each file the following values are extracted:
      - SIGMA_11, SIGMA_22, SIGMA_33 (from the "Principal values of the NMR shielding tensor" table)
      - DIA row: isotropic value and x, y, z values
      - PARA row: isotropic value and x, y, z values
      - SUM row: isotropic value and x, y, z values

    The results are saved as a CSV file with the columns:
      Filename, SIGMA_11, SIGMA_22, SIGMA_33,
      DIA_iso, DIA_x, DIA_y, DIA_z,
      PARA_iso, PARA_x, PARA_y, PARA_z,
      SUM_iso, SUM_x, SUM_y, SUM_z

    Parameters:
      folder (str): Folder containing the .cs files (default: current directory).
      output (str): Name of the output CSV file (default: "shielding.csv").

    Returns:
      results (list of tuples): Each tuple contains the extracted data for one file.
    """
    # Define a regex pattern for a floating‐point number (supports scientific notation)
    float_re = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"

    # Patterns for the principal shielding values
    sigma11_pat = re.compile(r"SIGMA_11\s+" + float_re)
    sigma22_pat = re.compile(r"SIGMA_22\s+" + float_re)
    sigma33_pat = re.compile(r"SIGMA_33\s+" + float_re)

    # Patterns for the shielding table rows.
    # Each row has an isotropic value plus x, y, z components.
    dia_pat = re.compile(r"^\s*DIA\s+" + float_re + r"\s+" + float_re + r"\s+" + float_re + r"\s+" + float_re, re.MULTILINE)
    para_pat = re.compile(r"^\s*PARA\s+" + float_re + r"\s+" + float_re + r"\s+" + float_re + r"\s+" + float_re, re.MULTILINE)
    sum_pat = re.compile(r"^\s*SUM\s+" + float_re + r"\s+" + float_re + r"\s+" + float_re + r"\s+" + float_re, re.MULTILINE)

    results = []
    # Loop over all .cs files in the specified folder
    for filepath in glob.glob(os.path.join(folder, "*cs")):
        electronic_data = {}
        with open(filepath, "r") as file:
            content = file.read()

        # Extract sigma values (take last occurrence if multiple are found)
        e_matches = sigma11_pat.findall(content)
        sigma11_val = e_matches[-1] if e_matches else None

        e_matches = sigma22_pat.findall(content)
        sigma22_val = e_matches[-1] if e_matches else None

        e_matches = sigma33_pat.findall(content)
        sigma33_val = e_matches[-1] if e_matches else None

        # Extract DIA row: isotropic value, x, y, z
        dia_match = dia_pat.search(content)
        if dia_match:
            dia_iso = dia_match.group(1)
            dia_x   = dia_match.group(2)
            dia_y   = dia_match.group(3)
            dia_z   = dia_match.group(4)
        else:
            dia_iso = dia_x = dia_y = dia_z = None

        # Extract PARA row: isotropic value, x, y, z
        para_match = para_pat.search(content)
        if para_match:
            para_iso = para_match.group(1)
            para_x   = para_match.group(2)
            para_y   = para_match.group(3)
            para_z   = para_match.group(4)
        else:
            para_iso = para_x = para_y = para_z = None

        # Extract SUM row: isotropic value, x, y, z
        sum_match = sum_pat.search(content)
        if sum_match:
            sum_iso = sum_match.group(1)
            sum_x   = sum_match.group(2)
            sum_y   = sum_match.group(3)
            sum_z   = sum_match.group(4)
        else:
            sum_iso = sum_x = sum_y = sum_z = None

        results.append((
            os.path.basename(filepath),
            sigma11_val, sigma22_val, sigma33_val,
            dia_iso, dia_x, dia_y, dia_z,
            para_iso, para_x, para_y, para_z,
            sum_iso, sum_x, sum_y, sum_z
        ))

    # Write the results to the specified CSV file
    header = ("Filename,SIGMA_11,SIGMA_22,SIGMA_33,"
              "DIA_iso,DIA_x,DIA_y,DIA_z,"
              "PARA_iso,PARA_x,PARA_y,PARA_z,"
              "SUM_iso,SUM_x,SUM_y,SUM_z\n")

    with open(output, "w") as out_file:
        out_file.write(header)
        for row in results:
            row_str = ",".join([str(item) if item is not None else "" for item in row])
            out_file.write(row_str + "\n")

    print(f"Shielding data extracted and saved to '{output}'.")
    return results

if __name__ == 'main':
    results = extract_shielding_info(folder="data_RESPECT/all_outs", output="rel_shielding_data.csv")

    df = pd.DataFrame(results, columns=["Filename", "SIGMA_11", "SIGMA_22", "SIGMA_33",
                                        "DIA_iso", "DIA_x", "DIA_y", "DIA_z",
                                        "PARA_iso", "PARA_x", "PARA_y", "PARA_z",
                                        "SUM_iso", "SUM_x", "SUM_y", "SUM_z"])
    df