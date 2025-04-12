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