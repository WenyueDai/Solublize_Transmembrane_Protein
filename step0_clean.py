import os
import subprocess
import json
from pymol import cmd  # type: ignore
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# conda activate mlfold
# To download pdb directly from pdb database
# curl -O https://files.rcsb.org/download/6LN2.pdb

folder_with_pdbs = "/home/eva/0_solubilize_transmembrane/6LN2" #Change it to the folder where pdb is downloaded
clean_pdbs_dir = os.path.join(folder_with_pdbs, "clean") #The folder with cleaned pdb
chains_to_keep = "A,B,C,D,E,F,G,H" #Chain to keep during pdb clean up. Change to A will delete all the other chains.
chain_to_cal_hydro_sasa = "A" #Chain to calculate hydrophobic and sasa, can analysis multiple chain like A,B

os.makedirs(clean_pdbs_dir, exist_ok=True)

# ========== STRUCTURE CLEANING FUNCTIONS ==========

def convert_cif_to_pdb(input_cif, output_pdb):
    """Use pymol load and save to convert cit to pdb

    Args:
        input_cif (_type_): xxx.cif file
        output_pdb (_type_): xxx.pdb file

    Returns:
        _type_: pdb file
    """
    try:
        cmd.load(input_cif, "structure")
        cmd.save(output_pdb, "structure")
        cmd.delete("all")
        return True
    except Exception as e:
        print(f"CIF → PDB failed: {input_cif}, {e}")
        cmd.delete("all")
        return False

def reformat_pdb_with_pdbfixer(input_pdb, output_pdb, ph=7.0):
    try:
        print(f"Fixing: {input_pdb}")
        fixer = PDBFixer(filename=input_pdb)

        # Find and remove nonstandard residues
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        # Remove waters
        fixer.removeHeterogens(keepWater=False)

        # Identify and add missing atoms
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Add hydrogens at specified pH
        fixer.addMissingHydrogens(pH=ph)

        # Sanity check
        if fixer.topology.getNumAtoms() == 0:
            print(f"No atoms found in {input_pdb}, skipping.")
            return False

        # Write the fixed structure
        with open(output_pdb, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        print(f"Saved fixed PDB to {output_pdb}")
        return True

    except Exception as e:
        print(f"PDBFixer error for {input_pdb}: {e}")
        return False

def clean_and_renumber_pdb(input_file, output_file, chains_to_keep=None):
    try:
        cmd.load(input_file, "protein")
        
        # Remove salt
        # cmd.remove("resn HOH or resn NA+CL+MG+CA+ZN or hetatm")
        cmd.remove("resn NA+CL+MG+CA+ZN or hetatm")
        
        # Delete the non-main alternative
        # Some atoms might save in more positions (side chain that wiggles between two shapes etc)
        cmd.remove("not alt ''+A")

        # If chain filtering is requested
        if chains_to_keep:
            chains = [c.strip() for c in chains_to_keep.replace(",", " ").split()]
            chain_filter = " or ".join([f"chain {c}" for c in chains])
            cmd.remove(f"not ({chain_filter})")

        atoms = cmd.get_model("all").atom
        if not atoms:
            print(f"No atoms found: {input_file}")
            cmd.delete("all")
            return False

        # Renumber residue index by starting from 1
        starting_residue = int(atoms[0].resi)
        cmd.alter("all", f"resi=str(int(resi) - {starting_residue - 1})")
        
        # Renumber atom serial number in ATOM/HETATM line
        atom_counter = 1
        cmd.alter("all", "serial=atom_counter; atom_counter+=1", space={'atom_counter': atom_counter})
        cmd.sort()
        cmd.save(output_file, "protein")
        cmd.delete("all")
        return os.path.exists(output_file)

    except Exception as e:
        print(f"PyMOL cleanup failed: {input_file}, {e}")
        cmd.delete("all")
        return False

def batch_process_structures(input_folder, clean_folder, chains_to_keep=None):
    for file in os.listdir(input_folder):
        if not (file.endswith(".pdb") or file.endswith(".cif")):
            continue
        input_path = os.path.join(input_folder, file)
        base = os.path.splitext(file)[0]
        intermediate_pdb = os.path.join(clean_folder, base + "_intermediate.pdb")
        cleaned_pdb = os.path.join(clean_folder, base + ".pdb")
        is_converted = False

        if file.endswith(".cif"):
            print(f"Converting CIF: {file}")
            if not convert_cif_to_pdb(input_path, intermediate_pdb):
                continue
            is_converted = True
        else:
            intermediate_pdb = input_path

        print(f"Reformatting: {intermediate_pdb}")
        if not reformat_pdb_with_pdbfixer(intermediate_pdb, intermediate_pdb, ph=7.0):
            continue

        print(f"PyMOL cleanup to: {cleaned_pdb}")
        success = clean_and_renumber_pdb(intermediate_pdb, cleaned_pdb, chains_to_keep)
        if success:
            print(f"Saved: {cleaned_pdb}")
        else:
            print(f"Failed: {cleaned_pdb}")

        if is_converted and os.path.exists(intermediate_pdb):
            os.remove(intermediate_pdb)
            print(f"Removed: {intermediate_pdb}")

    print("All structures processed.\n")
    
batch_process_structures(folder_with_pdbs, clean_pdbs_dir, chains_to_keep)

# === ANALYZE RELATIVE SASA + HYDROPHOBICITY ===
# Eisenberg hydrophobicity scale
def analyze_sasa_and_hydrophobicity(obj_name, output_txt_path, chains_to_analyze, sasa_threshold=0.5, hydro_threshold=0, hydro_sasa_threthod=0.2):
    hydrophobicity_eisenberg = {
        'ALA': 0.25, 'ARG': -1.76, 'ASN': -0.64, 'ASP': -0.72,
        'CYS': 0.04, 'GLN': -0.69, 'GLU': -0.62, 'GLY': 0.16,
        'HIS': -0.40, 'ILE': 0.73, 'LEU': 0.53, 'LYS': -1.10,
        'MET': 0.26, 'PHE': 0.61, 'PRO': -0.07, 'SER': -0.26,
        'THR': -0.18, 'TRP': 0.37, 'TYR': 0.02, 'VAL': 0.54
    }

    # Convert string like "B,C,D" to list ['B', 'C', 'D']
    chains = [c.strip() for c in chains_to_analyze.split(',') if c.strip()]
    chain_selection = " or ".join([f"chain {c}" for c in chains])
    selection_str = f"{obj_name} and ({chain_selection})"

    # Step 1: Calculate per-residue relative SASA
    sasa_dic = cmd.get_sasa_relative(selection=selection_str, state=1)

    # Step 2: Get unique residues (resi, resn, chain)
    unique_residues = set()
    cmd.iterate(selection_str, "unique_residues.add((resi, resn, chain))", space={"unique_residues": unique_residues})

    # Step 3: Process each residue
    residue_data = []
    chain_resi_list = {}
    skipped_missing_sasa = 0
    skipped_below_threshold = 0
    skipped_nonhydrophobic = 0

    for resi, resn, chain in sorted(unique_residues, key=lambda x: (int(x[0]), x[1], x[2])):
        sasa = None
        for (model, segi, c, r), s in sasa_dic.items():
            if c == chain and r == resi:
                sasa = s
                break

        if sasa is None:
            print(f"Skipping {chain}:{resi} ({resn}) — SASA not found")
            skipped_missing_sasa += 1
            continue

        if sasa < sasa_threshold:
            skipped_below_threshold += 1
            continue

        resn_upper = resn.upper()
        hydro_score = hydrophobicity_eisenberg.get(resn_upper)
        if hydro_score is None:
            continue

        if hydro_score <= hydro_threshold:
            skipped_nonhydrophobic += 1
            continue

        hydro_sasa = hydro_score * sasa

        if hydro_sasa > hydro_sasa_threthod:
            print(f"Keeping {chain}:{resi} ({resn}) — SASA={sasa:.2f}, Hydrophobicity={hydro_score:.2f}")
            residue_data.append((chain, resn_upper, resi, sasa, hydro_score, hydro_sasa))

    print(f"\nSummary:")
    print(f" Kept residues: {len(residue_data)}")
    print(f" Skipped (no SASA): {skipped_missing_sasa}")
    print(f" Skipped (below SASA threshold): {skipped_below_threshold}")
    print(f" Skipped (non-hydrophobic): {skipped_nonhydrophobic}")

    residue_data.sort(key=lambda x: (-x[5], -x[4], -x[3], x[0], int(x[2]), x[1]))
    for chain, _, resi, _, _, _ in residue_data:
        chain_resi_list.setdefault(chain, []).append(int(resi))
    for chain in sorted(chain_resi_list):
        resi_list = " ".join(str(r) for r in chain_resi_list[chain])
        print(f'  "{chain}": [{resi_list}]')
        resi_list_pymol = "+".join(str(r) for r in chain_resi_list[chain])
        print(f'  "{chain}": [{resi_list_pymol}]')

    with open(output_txt_path, 'w') as f:
        f.write("sasa_threshold=0.5, hydro_threshold=0, hydro_sasa_threthod=0.2")
        f.write("Chain:Resi:Resn\tRela_SASA*Hydro\tRelative_SASA\tHydrophobicity\n")
        for chain, resn_upper, resi, sasa, hydro_score, hydro_sasa in residue_data:
            f.write(f"{chain}:{resi}:{resn_upper}\t{hydro_sasa:.2f}\t{sasa:.2f}\t{hydro_score:.2f}\n")
        for chain in sorted(chain_resi_list):
            resi_list = " ".join(str(r) for r in chain_resi_list[chain])
            f.write(f'  "{chain}": [{resi_list}]')

    print(f"[DEBUG] Written results to: {output_txt_path}")
    cmd.delete("all")


# === RUN ANALYSIS FOR EACH CLEANED PDB ===

# Loop through cleaned files
for cleaned_file in os.listdir(clean_pdbs_dir):
    if cleaned_file.endswith(".pdb"):
        pdb_path = os.path.join(clean_pdbs_dir, cleaned_file)
        obj_name = os.path.splitext(cleaned_file)[0]
        output_txt = pdb_path.replace(".pdb", "_sasa_hydro.txt")
        cmd.load(pdb_path, obj_name)
        analyze_sasa_and_hydrophobicity(obj_name, output_txt, chains_to_analyze=chain_to_cal_hydro_sasa)


