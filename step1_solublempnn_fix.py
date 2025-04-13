import os
import subprocess
import json

# conda activate mlfold

# To select the residue to mutate, in pymol, open clean.pdb
# iterate sele and name CA, resi
# We use "--use_soluble_model" as solubleMPNN

# ========== USER SETTINGS ==========

folder_with_pdbs = "/home/eva/0_solubilize_transmembrane/7KP4_2" # The same path from step0, no change
clean_pdbs_dir = os.path.join(folder_with_pdbs, "clean") # The same path from step0, no change
output_dir = os.path.join(folder_with_pdbs, "output")
os.makedirs(output_dir, exist_ok=True)

# Always include all the chain in clean.pdb, otherwise aftermpnn the non-defined chain will be deleted
chain_to_design = "A B" 

# Manually define fixed residue indices per chain TO DESIGN
# If you want to define fixed resideu to EXCLUDE FROM DESIGN, TURN OFF "--specify_non_fixed" in Make fixed position step
#only designing residues listed here, chain A
fixed_positions_per_chain = "36 2 32 85 143 67 92 82 4 172 126 62 66 15 130 89 81 12 129 123 116 119 140 95 19 176 58 134 23 5 47 165 22 78 75 33 56 158, ," 

mpnn_base = "/home/eva/ProteinMPNN"
parse_script = os.path.join(mpnn_base, "helper_scripts/parse_multiple_chains.py")
assign_script = os.path.join(mpnn_base, "helper_scripts/assign_fixed_chains.py")
fix_position_script = os.path.join(mpnn_base, "helper_scripts/make_fixed_positions_dict.py")
run_script = os.path.join(mpnn_base, "protein_mpnn_run.py")

num_seq = "20"
sampling_temp = "0.2"
omit_aas = "CWY" # Residue you dont want to use
seed = "42" # To ensure the reproducability of solublempnn
batch_size = "1"

parsed_jsonl = os.path.join(output_dir, "parsed_pdbs.jsonl")
assigned_jsonl = os.path.join(output_dir, "assigned_pdbs.jsonl")
fixed_jsonl = os.path.join(output_dir, "fixed_positions.jsonl")

print("Parsing chains...")
subprocess.run([
    "python", parse_script,
    "--input_path", clean_pdbs_dir,
    "--output_path", parsed_jsonl
], check=True)

print("Assign fixed chains...")
subprocess.run([
    "python", assign_script,
    "--input_path", parsed_jsonl,
    "--output_path", assigned_jsonl,
    "--chain_list", chain_to_design,
], check=True)

print("Make fixed position...")
subprocess.run([
    "python", fix_position_script,
    "--input_path", parsed_jsonl,
    "--output_path", fixed_jsonl,
    "--chain_list", chain_to_design,
    "--position_list", fixed_positions_per_chain,
    "--specify_non_fixed"
], check=True)

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
