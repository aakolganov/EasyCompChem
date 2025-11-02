import glob, os
#extract homo-lumo gap from the ORCA stuff


def extract_homo_lumo_gaps(folder=".", output="homo_lumo_gaps.csv"):
    """
    Extracts HOMO-LUMO gaps from the *last* 'ORBITAL ENERGIES' section in ORCA .out files.

    Parameters:
        folder (str): Path to the folder containing .out files.
        output (str): Filename for the output CSV file.

    Returns:
        results (list of tuples): Each tuple contains (filename, HOMO energy, LUMO energy, gap in eV)
    """
    # Match all ORBITAL ENERGIES sections (up to SPIN DOWN)
    orbital_section_pattern = re.compile(
        r"-{10,}\s*ORBITAL ENERGIES\s*-{10,}(.*?)(?:SPIN DOWN ORBITALS|(?=-{10,}))",
        re.DOTALL
    )
    orbital_energy_line_pattern = re.compile(
        r"\s*\d+\s+([01]\.0000)\s+-?\d+\.\d+\s+(-?\d+\.\d+)"
    )

    results = []
    search_pattern = os.path.join(folder, "*.out")

    for filepath in glob.glob(search_pattern):
        homo_energy = None
        lumo_energy = None
        gap = None

        with open(filepath, "r") as file:
            content = file.read()

            # Find all orbital energy blocks and select the last one
            orbital_matches = orbital_section_pattern.findall(content)
            if orbital_matches:
                last_orbital_block = orbital_matches[-1]
                orbital_lines = orbital_energy_line_pattern.findall(last_orbital_block)

                for occ, energy in orbital_lines:
                    energy = float(energy)
                    if occ == "1.0000":
                        homo_energy = energy
                    elif occ == "0.0000" and homo_energy is not None:
                        lumo_energy = energy
                        gap = lumo_energy - homo_energy
                        break

        results.append((os.path.basename(filepath), homo_energy, lumo_energy, gap))

    # Write to CSV
    with open(output, "w") as out_file:
        out_file.write("Filename,HOMO (eV),LUMO (eV),HOMO-LUMO Gap (eV)\n")
        for filename, homo, lumo, gap in results:
            out_file.write(f"{filename},{homo},{lumo},{gap}\n")

    print(f"HOMO-LUMO gaps saved to {output}")
    return results

if __name__ == 'main':
    extract_homo_lumo_gaps(folder='data_from_ORCA/all_outs_adducts', output="homo_lumo_gaps.csv")
