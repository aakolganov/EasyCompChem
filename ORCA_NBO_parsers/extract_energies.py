import os, re, glob

def extract_energies(folder=".", output="energies.csv"):
    """
    Extracts electronic energy and Gibbs free energy values from ORCA .out files in the specified folder.

    Parameters:
        folder (str): Path to the folder containing .out files.
        output (str): Filename for the output CSV file.

    Returns:
        results (list of tuples): A list of tuples, each containing
                                  (filename, electronic energy, free energy).
    """
    # Compile regular expressions for energy extraction
    electronic_energy_pattern = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
    free_energy_pattern = re.compile(r"Final Gibbs free energy\s+.*?(-?\d+\.\d+)\s+Eh")

    results = []
    search_pattern = os.path.join(folder, "*.out")

    # Loop over all .out files in the specified folder
    for filepath in glob.glob(search_pattern):
        electronic_energy = None
        free_energy = None

        with open(filepath, "r") as file:
            content = file.read()

            # Extract the last occurrence of each energy type
            e_matches = electronic_energy_pattern.findall(content)
            if e_matches:
                electronic_energy = e_matches[-1]

            f_matches = free_energy_pattern.findall(content)
            if f_matches:
                free_energy = f_matches[-1]

        # Store the results using the basename of the file for clarity
        results.append((os.path.basename(filepath), electronic_energy, free_energy))

    # Write the results to the specified CSV file
    with open(output, "w") as out_file:
        out_file.write("Filename,Electronic Energy,Free Energy\n")
        for filename, e_energy, g_energy in results:
            out_file.write(f"{filename},{e_energy},{g_energy}\n")

    print(f"Results saved to {output}")
    return results

if __name__ == 'main':
    # Extract energies from .out files in the specified folder and save results'
    energies_acids = extract_energies(folder='data_from_ORCA/all_outs_acids', output="acids_ens.csv")

    # Optionally, display the results
    # for item in results:
    #   print(item)