"""
Process and extract phosphorus chemical shift data matrices from ORCA output files.

This module processes ORCA computational chemistry output files (`*.out`)
to extract diagonalized chemical shift matrices for phosphorus nuclei. It then
computes differences from standard or BEA reference values and saves the results
to a CSV file. The module supports custom reference values as input.

The function handles the identification of matrices and optional reference
values, performs calculations, and outputs results in a structured format.
"""

import os, re, glob, pandas as pd

def extract_chem_shift_matrix_ORCA(folder=".",
                                   output="diag_matrix_diff_P.csv",
                                   ref_std=None,
                                   ref_bea=None):
    """
    Scan every *.out in *folder*, locate the diagonalised sT*s matrix for
    phosphorus nuclei, subtract reference values, and write a CSV.

    * If the filename contains 'BEA' (case-insensitive) the BEA reference
      set is used, otherwise the standard reference set.
    * Provide custom references with *ref_std* / *ref_bea*
      (each a dict with keys 'sDSO','sPSO','Total', value=list[4 floats]).
    """

    # ---------- reference data ----------------------------------------------
    ref_std = ref_std or {               # “normal” files
        "sDSO":  [962.856,  965.926, 971.039, 966.607],
        "sPSO":  [-758.375, -760.021, -551.868, -690.088],
        "Total": [204.480,  205.905, 419.171, 276.519]
    }

    ref_bea = ref_std # or {               # files whose name contains BEA
    #     "sDSO":  [964.279, 967.333  ,972.002, 967.871],   # <-- dummy placeholders
    #     "sPSO":  [ -767.717, -770.103, -559.517,  -699.112],
    #     "Total": [196.562, 197.230   , 412.485 , 268.759]
    # }

    # ---------- regex helpers -----------------------------------------------
    float_re = r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
    nucleus_header_pat = re.compile(
        r"^\s*-+\s*\n\s*Nucleus\s+(\d{1,3}P)\s*:\s*\n\s*-+",
        re.MULTILINE
    )
    row_pat = {name: re.compile(
        fr"{name}\s+{float_re}\s+{float_re}\s+{float_re}\s+iso=\s+{float_re}")
        for name in ("sDSO", "sPSO", "Total")
    }

    records = []

    # ---------- loop over all ORCA *.out files -------------------------------
    for path in glob.glob(os.path.join(folder, "*.out")):
        with open(path, "r", errors="ignore") as f:
            text = f.read()

        blocks = list(nucleus_header_pat.finditer(text))
        if not blocks:
            continue  # no phosphorus in this file

        # choose reference set
        ref = ref_bea if re.search("bea", os.path.basename(path), flags=re.I) else ref_std

        # walk through every P block
        for i, hdr in enumerate(blocks):
            start = hdr.start()
            end   = blocks[i+1].start() if i+1 < len(blocks) else len(text)
            chunk = text[start:end]

            diffs = {}
            for name, pat in row_pat.items():
                m = pat.search(chunk)
                if m:
                    vals  = [float(m.group(g)) for g in range(1,5)]
                    diffs[name] = [ref[name][k] - vals[k] for k in range(4)]
                else:
                    diffs[name] = [None]*4

            records.append({
                "Filename": os.path.basename(path),
                #"Nucleus":  hdr.group(1),
                "d_sDSO_v1":   diffs["sDSO"][0], "d_sDSO_v2":   diffs["sDSO"][1],
                "d_sDSO_v3":   diffs["sDSO"][2], "d_sDSO_iso":  diffs["sDSO"][3],
                "d_sPSO_v1":   diffs["sPSO"][0], "d_sPSO_v2":   diffs["sPSO"][1],
                "d_sPSO_v3":   diffs["sPSO"][2], "d_sPSO_iso":  diffs["sPSO"][3],
                "d_Total_v1":  diffs["Total"][0],"d_Total_v2":  diffs["Total"][1],
                "d_Total_v3":  diffs["Total"][2],"d_Total_iso": diffs["Total"][3],
            })

    # ---------- write CSV ----------------------------------------------------
    df = pd.DataFrame(records)
    df.to_csv(output, index=False)
    print(f"{len(df)} phosphorus block(s) processed → '{output}'")
    return records

if __name__ == "__main__":
    extract_chem_shift_matrix_ORCA(folder="data_ORCA_NMR/all_outs_ZORA",output="NMR_ORCA_ZORA.csv")