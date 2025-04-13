import os
import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, RestrictToRepacking

pyrosetta.init("-mute all")

# conda activate getcontact

# ========== USER CONFIG ==========
fasta_path = "/home/eva/0_solubilize_transmembrane/7KP4_2/output/seqs/7KP4.fa"
pdb_path = "/home/eva/0_solubilize_transmembrane/7KP4_2/7KP4.pdb"
output_folder = "/home/eva/0_solubilize_transmembrane/6LN2/output/seqs/mutated_pdb"

def one_letter_to_three(letter):
    table = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    return table.get(letter)

def parse_fasta_sequences_with_chain_split(fasta_path):
    """
    Parse sequences and headers, splitting sequence by chain using '/'
    """
    sequences = []
    current_header = ""
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_header = line.strip()
            else:
                full_sequence = line.strip()
                split_by_chain = full_sequence.split('/')
                sequences.append((current_header, split_by_chain))
    return sequences

def extract_sample_number(header):
    """
    Extract sample number from header like:
    >T=0.2, sample=4, score=1.9392 ...
    """
    for part in header.split(','):
        part = part.strip()
        if part.startswith("sample="):
            return part.split('=')[1]
    return "wt"


def get_chain_order_from_pdb(pdb_path):
    pose = pose_from_pdb(pdb_path)
    pdb_info = pose.pdb_info()
    print(f"pdb_info: {pdb_info}")
    chain_order = []
    seen = set()
    for i in range(1, pose.total_residue() + 1):
        chain = pdb_info.chain(i)
        if chain not in seen:
            chain_order.append(chain)
            seen.add(chain)
    return chain_order

def mutate_pose_chain(pose, target_sequence, chain_id):
    pdb_info = pose.pdb_info()
    residue_ids = [i for i in range(1, pose.total_residue() + 1) if pdb_info.chain(i) == chain_id]
    
    if len(residue_ids) != len(target_sequence):
        raise ValueError(f"Mismatch for chain {chain_id}: {len(residue_ids)} in PDB vs {len(target_sequence)} in sequence")

    for i, resi in enumerate(residue_ids):
        target_aa = target_sequence[i]
        current_aa = pose.residue(resi).name1()
        if target_aa != current_aa:
            aa3 = one_letter_to_three(target_aa)
            if aa3:
                new_res = ResidueFactory.create_residue(pose.residue_type_set_for_pose().name_map(aa3))
                pose.replace_residue(resi, new_res, True)

def repack_pose(pose):
    scorefxn = get_score_function()
    mover = PackRotamersMover(scorefxn)
    tf = TaskFactory()
    tf.push_back(InitializeFromCommandline())
    tf.push_back(RestrictToRepacking())
    mover.task_factory(tf)
    mover.apply(pose)

def process_pdb_with_chain_sequences(pdb_path, chain_sequences, chain_ids, output_path):
    pose = pose_from_pdb(pdb_path)
    
    for chain_id, seq in zip(chain_ids, chain_sequences):
        mutate_pose_chain(pose, seq, chain_id)

    repack_pose(pose)
    pose.dump_pdb(output_path)
    print(f"Saved mutated structure: {output_path}")

def batch_mutate_all_fasta_sequences(fasta_path, pdb_path, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    all_chain_sequences = parse_fasta_sequences_with_chain_split(fasta_path)
    chain_ids = get_chain_order_from_pdb(pdb_path)

    for i, (header, chain_seqs) in enumerate(all_chain_sequences):
        if len(chain_seqs) != len(chain_ids):
            print(f"Sequence {i+1} has {len(chain_seqs)} chains but PDB has {len(chain_ids)} — skipping.")
            continue

        sample_id = extract_sample_number(header)
        output_filename = f"{os.path.splitext(os.path.basename(pdb_path))[0]}_sample{sample_id}.pdb"
        out_path = os.path.join(output_folder, output_filename)

        try:
            process_pdb_with_chain_sequences(pdb_path, chain_seqs, chain_ids, out_path)
        except Exception as e:
            print(f"Failed to process sequence sample={sample_id}: {e}")


batch_mutate_all_fasta_sequences(fasta_path, pdb_path, output_folder)
