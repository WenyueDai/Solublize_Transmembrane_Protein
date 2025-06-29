Great — let’s break this README into a **friendly step-by-step tutorial style**, with clear sections for each stage and beginner-friendly guidance. Here’s a polished version you can use:

---

# **Transmembrane Protein Solubilization & Repacking Pipeline**

This project provides a pipeline for redesigning and repacking transmembrane proteins to make them soluble and structurally optimized.

---

# **Pipeline Steps Overview**

| Step | Script                     | Purpose                                                                             |
| ---- | -------------------------- | ----------------------------------------------------------------------------------- |
| 0    | `step0_clean.py`           | Clean and renumber PDB/CIF files; identify residues with high hydrophobicity × SASA |
| 1    | `step1_solublempnn_fix.py` | Use ProteinMPNN (soluble model) to redesign the sequence, fixing specified residues |
| 2    | `step2_pyrosetta_map.py`   | Repack the redesigned sequence into the cleaned PDB structure using PyRosetta       |
| 3    | *(planned)*                | CHAI-1 structural prediction                                                        |
| 4    | *(planned)*                | OpenMM or membrane positioning                                                      |
| 5    | *(planned)*                | Pipeline integration with dynamic residue cutoff selection                          |
| 6    | *(planned)*                | Full reproducibility tracking                                                       |

---

# **Beginner-Friendly Step-by-Step Guide**

---

## **Step 0 — Set Up Your Linux Development Environment (WSL + Conda + VS Code)**

**Windows 11** users can use WSL (Windows Subsystem for Linux) to get a Linux-like environment:

1. **Open PowerShell** (as administrator):

   ```powershell
   wsl --install
   ```

   (this installs Ubuntu by default)

2. **Restart your computer**

3. In the Ubuntu terminal, update packages:

   ```bash
   sudo apt update && sudo apt upgrade -y
   ```

4. **Download and install VS Code** from [Visual Studio Code](https://code.visualstudio.com)

   * During installation, check:

     * *Add to PATH*
     * *Add “Open with Code” to context menu*
     * *Register VS Code as default editor*

5. In VS Code:

   * Press `Ctrl+Shift+X`
   * Install these extensions:

     * **WSL**
     * **Python**
     * **Jupyter**
     * **Pylance**

1. You are now ready to code in Linux via VS Code!

---

## **Step 1 — Prepare Your Conda Environment for PDB Cleaning & ProteinMPNN**

1. In your Ubuntu terminal:

   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   source ~/.bashrc
   ```

2. Create your project folder:

   ```bash
   mkdir -p ~/protein_project && cd ~/protein_project
   ```

3. Set up the environment:

   ```bash
   conda create -n mlfold python=3.10 -y
   conda activate mlfold
   ```

4. Install required packages:

   ```bash
   conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
   pip install numpy biopython pandas tqdm matplotlib
   conda install -c conda-forge openmm pdbfixer pymol-open-source
   ```

5. **Download ProteinMPNN**:

   ```bash
   git clone https://github.com/dauparas/ProteinMPNN.git
   ```

---

## **Step 2 — Clean Your PDB/CIF Files (`step0_clean.py`)**

1. Place your downloaded `.pdb` or `.cif` file in your working folder (e.g., `~/protein_project/inputs`).

   * You can download directly, for example:

     ```bash
     curl -O https://files.rcsb.org/download/7kp4.pdb
     ```

2. Edit the `step0_clean.py` script to point to your folder (e.g., `/home/eva/protein_project/inputs`) and set:

   * `chains_to_keep`
   * `sasa_threshold`
   * `hydro_threshold`

3. Run the cleaner:

   ```bash
   conda activate mlfold
   python step0_clean.py
   ```

1. *This will create cleaned structures in a `clean` folder and report residues with high hydrophobicity × SASA.*

---

## **Step 3 — Redesign with ProteinMPNN (`step1_solublempnn_fix.py`)**

1. Edit `step1_solublempnn_fix.py`:

   * Set `folder_with_pdbs` to the cleaned folder path
   * Adjust `chain_to_design` and `fixed_positions_per_chain` as needed

2. Run:

   ```bash
   conda activate mlfold
   python step1_solublempnn_fix.py
   ```

1. *This will redesign sequences with the ProteinMPNN soluble model and output new sequences to an `output` folder.*

---

## **Step 4 — Repack with PyRosetta (`step2_pyrosetta_map.py`)**

1. Install PyRosetta following instructions for [PyRosetta Windows 10+](https://www.pyrosetta.org/downloads/windows-10)

2. Create and activate the environment:

   ```bash
   conda create -n pyrosetta python=3.10 -y
   conda activate pyrosetta
   pip install pyrosetta-installer
   python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
   ```

3. Edit `step2_pyrosetta_map.py`:

   * Update `fasta_path` to the `.fa` from step 1
   * Update `pdb_path` to the cleaned structure

4. Run:

   ```bash
   conda activate pyrosetta
   python step2_pyrosetta_map.py
   ```

1. *This repacks the redesigned sequence into the PDB coordinates.*

---

## **Step 5 — Future Plans**

1. Integrate CHAI-1 structural prediction
2. Add OpenMM or membrane-positioning functionality
3. Allow dynamic cutoff residue selection across the entire workflow
4. Improve reproducibility tracking

---

# **Basic Command-Line Reminders**

* **activate environment**:

  ```bash
  conda activate mlfold
  ```
* **list files including hidden**:

  ```bash
  ls -a
  ```
* **check GPU**:

  ```bash
  nvidia-smi
  ```
* **find all .pdb recursively**:

  ```bash
  find . -name "*.pdb"
  ```
* **print working directory**:

  ```bash
  pwd
  ```
* **make a new folder**:

  ```bash
  mkdir myfolder
  ```

---
