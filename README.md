
# Description
1. step0_clean.py clean and renumber pdb/cif
2. step1_solublempnn_fix.py Use solubleMPNN to redesign transmembrane protein (non-fixed position)
3. step2_pyrosetta_map.py Use pyrosetta to repack sequence to the input pdb
4. Optional for future: chai-1 structural prediction

# Instruction for new starters
## step 1. VScode with wsl

1. Search for powershell, run as administrator.
2. wsl --install #Run this to install wsl and ubuntus (default)
3. Set ubuntus as default terminal backend
4. reboot
5. In ubuntus terminal run: sudo apt update && sudo apt upgrade -y #update system to latest package
6. Go to: https://code.visualstudio.com, download VScode
7. During install: "Add to PATH""Add ‘Open with Code’ to context menu""Register VS Code as default editor"
8. In VScode, go to extension tab: Ctrl+Shift+X
9. Search for WSL, Python, Jupyter, Pylance
10. In ubuntu terminal, create test folder: 
mkdir -p ~/test_project
cd ~/test_project
11. Download conda wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
12. Run installer: bash Miniconda3-latest-Linux-x86_64.sh
13. Restart terminal: source ~/.bashrc
14. Set up git: sudo apt install git -y
git config --global user.name "Your Name"
git config --global user.email "your@email.com"

## step 2. conda env mlfold for pdb preparation and solublempnn
conda create -n mlfold python=3.10 -y
conda activate mlfold
conda install pytorch torchvision torchaudio cudatoolkit=11.3 -c pytorch
pip install numpy biopython pandas tqdm matplotlib

## step 3. conda env getcontact for pyrosetta repack
Go to https://www.pyrosetta.org/downloads/windows-10 and fill out the form for a free academic license for PyRosetta
conda install pyrosetta --channel https://<username>:<password>@conda.gralab.jhu.edu
