import os, shutil

def create_crest_files(folder, charge):
    # Define the CREST script template
    crest_script_template = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=24GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 48
#SBATCH --partition=genoa
#SBATCH --time=48:00:00
#SBATCH -J CREST

~/bin/xtb {filename} --gfn2 --chrg {charge} --opt -T 48 > xtb_output.log
~/bin/crest xtbopt.xyz --gfn2 --chrg {charge} -T 48 > crest_output.log
"""

    # Get all .xyz files in the specified folder
    xyz_files = [f for f in os.listdir(folder) if f.endswith('.xyz')]

    for xyz_file in xyz_files:
        # Create a folder for each .xyz file
        folder_name = os.path.splitext(xyz_file)[0]
        new_folder_path = os.path.join(folder, folder_name)
        os.makedirs(new_folder_path, exist_ok=True)

        # Create the CREST script
        crest_script_content = crest_script_template.format(filename=xyz_file, charge=charge)
        script_path = os.path.join(new_folder_path, "run_crest.sh")

        with open(script_path, 'w') as script_file:
            script_file.write(crest_script_content)

        # Move the .xyz file into the folder
        os.rename(os.path.join(folder, xyz_file), os.path.join(new_folder_path, xyz_file))

    print("Folders and CREST scripts generated successfully.")


def extract_xyz_from_CREST(crest_folder, output_folder):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Process each folder in CREST_FOLDER
    for folder in os.listdir(crest_folder):
        folder_path = os.path.join(crest_folder, folder)

        # Skip if it's not a folder
        if not os.path.isdir(folder_path):
            continue

        # Check if the folder name ends with "_OK"

        # Copy and rename crest_best.xyz
        source_file = os.path.join(folder_path, "crest_best.xyz")
        if os.path.exists(source_file):
            new_name = folder + ".xyz"
            destination_file = os.path.join(output_folder, new_name)
            shutil.copy(source_file, destination_file)

    print(f"Processed all folders. Results are in '{output_folder}'.")

#example_usage

if __name__ == "__main__":
    create_crest_files('path/to/your/xyzs/', 0)
    crest_folder = 'path/to/your/crest_output/'
    output_folder = 'path/to/your/output_folder'
    extract_xyz_from_CREST(crest_folder, output_folder)
