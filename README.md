
# Description
1. step0_clean.py clean and renumber pdb/cif
2. step1_solublempnn_fix.py Use solubleMPNN to redesign transmembrane protein (non-fixed position)
3. step2_pyrosetta_map.py Use pyrosetta to repack sequence to the input pdb
4. Optional for future: chai-1 structural prediction
5. Optional OPENMM OR Positioning of proteins in membranes

# Instruction for new starters
## step 0. Common bash script used in terminal
1. conda create -n $ENV_NAME python=$PYTHON_VER -y
2. conda activate $ENV_NAME
3. conda deactivate $ENV_NAME
4. cd test # go to test folder
5. cd .. # go one directory level up
6. cd ~ # go to your home directory
5. ls -a # list all files (including hidden ones)
6. pwd # Print the current working directory (where you are).
6. find . -name "*.pdb" # Find All .pdb Files Recursively
7. mkdir test_folder # Make a test_folder folder
8. nvidia-smi # Check GPU availablility


## step 1. VScode with wsl
1. Search for powershell, run as administrator.
2. wsl --install #Run this to install wsl and ubuntus (default)
3. Set ubuntus as default terminal backend (only work for window 11 above)
4. reboot the computer
5. In ubuntus terminal run: sudo apt update && sudo apt upgrade -y #update system to latest package
6. Go to: https://code.visualstudio.com, download VScode
7. During install: "Add to PATH""Add ‘Open with Code’ to context menu""Register VS Code as default editor"
8. In VScode, go to extension tab: Ctrl+Shift+X
9. Search for WSL, Python, Jupyter, Pylance
10. In ubuntu terminal, create test folder: 
mkdir -p ~/test_project cd ~/test_project
11. Download conda wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
12. Run installer: bash Miniconda3-latest-Linux-x86_64.sh
13. Restart terminal: source ~/.bashrc
14. Set up git: sudo apt install git -y
git config --global user.name "Your Name"
git config --global user.email "your@email.com"

## step 2. conda env mlfold for pdb preparation and solublempnn
1. conda create -n mlfold python=3.10 -y
2. conda activate mlfold 
3. conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
4. pip install numpy biopython pandas tqdm matplotlib
5. conda install -c conda-forge openmm pdbfixer -y
6. conda install -c conda-forge pymol-open-source -y
7. Download ProteinMPNN to local: git clone https://github.com/dauparas/ProteinMPNN.git
8. To run the step0_clean.py and step1_solublempnn.py:
- conda activate mlfold
- python /home/eva/0_solubilize_transmembrane/step0_clean.py
- python /home/eva/0_solubilize_transmembrane/step1_solublempnn_fix.py

## step 3. conda env getcontact for pyrosetta repack
Go to https://www.pyrosetta.org/downloads/windows-10 for pyrosetta
1. conda create -n pyrosetta python=3.10 -y
2. conda activate pyrosetta
3. pip install pyrosetta-installer 
4. python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
5. To run the step2_pyrosetta_map.py:
- conda activate pyrosetta
- python /home/eva/0_solubilize_transmembrane/step2_pyrosetta_map.py
