import os
import subprocess
import json
from pymol import cmd  # type: ignore
from pdbfixer import PDBFixer
from openmm.app import PDBFile

"""This is solublempnn to find the sequence for both capping, and the 14AA extended from the capping sequence
"""

#conda activate mlfold
# curl -O https://files.rcsb.org/download/7SSZ.pdb

import os
import subprocess
import json
from pymol import cmd  # type: ignore
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# ========== USER SETTINGS ==========

folder_with_pdbs = "/home/eva/0_solubilize_transmembrane"
clean_pdbs_dir = os.path.join(folder_with_pdbs, "clean")
output_dir = os.path.join(folder_with_pdbs, "output")
chains_to_keep = "A,B,C,D,E,F,G,H" # Write down all chain "A,B,C,D,E,F,G,H" to keep all chain, or individual chain

# Manually define fixed residue indices per chain
fixed_positions_per_chain = {
    "A": [15, 16, 17, 18, 19, 20]
}

mpnn_base = "/home/eva/ProteinMPNN"
parse_script = os.path.join(mpnn_base, "helper_scripts/parse_multiple_chains.py")
assign_script = os.path.join(mpnn_base, "helper_scripts/assign_fixed_chains.py")
run_script = os.path.join(mpnn_base, "protein_mpnn_run.py")

num_seq = "10"
sampling_temp = "0.2"
omit_aas = "CWY"
seed = "37"
batch_size = "1"

os.makedirs(clean_pdbs_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# ========== STRUCTURE CLEANING FUNCTIONS ==========

def convert_cif_to_pdb(input_cif, output_pdb):
    try:
        cmd.load(input_cif, "structure")
        cmd.save(output_pdb, "structure")
        cmd.delete("all")
        return True
    except Exception as e:
        print(f"CIF → PDB failed: {input_cif}, {e}")
        cmd.delete("all")
        return False

def reformat_pdb_with_pdbfixer(input_pdb, output_pdb):
    try:
        fixer = PDBFixer(filename=input_pdb)
        if fixer.topology.getNumAtoms() == 0:
            print(f"No atoms in {input_pdb}, skipping.")
            return False
        with open(output_pdb, "w") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        return True
    except Exception as e:
        print(f"PDBFixer error: {input_pdb}, {e}")
        return False

def clean_and_renumber_pdb(input_file, output_file, chains_to_keep=None):
    try:
        cmd.load(input_file, "protein")
        cmd.remove("resn HOH or resn NA+CL+MG+CA+ZN or hetatm")
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

        starting_residue = int(atoms[0].resi)
        cmd.alter("all", f"resi=str(int(resi) - {starting_residue - 1})")
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

def batch_process_structures(input_folder, clean_folder, chains_to_keep="A"):
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
        if not reformat_pdb_with_pdbfixer(intermediate_pdb, intermediate_pdb):
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

    print("🎉 All structures processed.\n")

# ========== STEP 1: CLEAN PDBs ==========
batch_process_structures(folder_with_pdbs, clean_pdbs_dir, chains_to_keep)

# ========== STEP 2: PARSE PDB CHAINS ==========
parsed_jsonl = os.path.join(output_dir, "parsed_pdbs.jsonl")
assigned_jsonl = os.path.join(output_dir, "assigned_pdbs.jsonl")
fixed_jsonl = os.path.join(output_dir, "fixed_positions_caps.jsonl")

print("Parsing chains...")
subprocess.run([
    "python", parse_script,
    "--input_path", clean_pdbs_dir,
    "--output_path", parsed_jsonl
], check=True)

# ========== STEP 3: MANUAL FIXED POSITIONS ==========
fixed_dict = {}
assigned_dict = {}

chain_list = [c.strip() for c in chains_to_keep.split(",") if c.strip()]

with open(parsed_jsonl) as f:
    parsed_entries = [json.loads(line.strip()) for line in f]

for entry in parsed_entries:
    name = entry["name"]
    fixed_dict[name] = {}
    assigned_dict[name] = [chain_list, []]
    for c in chain_list:
        seq = entry.get(f"seq_chain_{c}", "")
        if not seq:
            print(f"Chain {c} not found in {name}, skipping.")
            continue
        sequence_length = len(seq)
        fixed_residues = [r for r in fixed_positions_per_chain.get(c, []) if 1 <= r <= sequence_length]
        fixed_dict[name][c] = sorted(fixed_residues)
        print(f"Fixed positions for {name} chain {c}: {fixed_dict[name][c]}")

with open(fixed_jsonl, "w") as f:
    f.write(json.dumps(fixed_dict) + "\n")

# ========== STEP 4: ASSIGN CHAIN IDs ==========
print("Assigning chains...")
subprocess.run([
    "python", assign_script,
    "--input_path", parsed_jsonl,
    "--output_path", assigned_jsonl,
    "--chain_list", chains_to_keep.replace(",", " ")
], check=True)

# ========== STEP 5: RUN PROTEINMPNN ==========
print("Running ProteinMPNN...")
subprocess.run([
    "python", run_script,
    "--jsonl_path", parsed_jsonl,
    "--chain_id_jsonl", assigned_jsonl,
    "--fixed_positions_jsonl", fixed_jsonl,
    "--out_folder", output_dir,
    "--use_soluble_model",
    "--num_seq_per_target", num_seq,
    "--sampling_temp", sampling_temp,
    "--omit_AAs", omit_aas,
    "--seed", seed,
    "--batch_size", batch_size
], check=True)

print("ProteinMPNN run complete.")
