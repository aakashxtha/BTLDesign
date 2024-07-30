# BTL Protein Design Project

## Overview
This project focuses on computational protein design using PyRosetta, specifically targeting the 1BTL protein structure. It includes scripts for designing mutations around specified residues, running parallel design jobs, and analyzing the results.

## Table of Contents
1. [Files](#files)
2. [Dependencies](#dependencies)
3. [Installation](#installation)
4. [Usage](#usage)
   - [Basic Design](#basic-design)
   - [Parallel Design](#parallel-design)
   - [Design from Relaxed Structure](#design-from-relaxed-structure)
5. [Output](#output)
6. [Design Process](#design-process)
7. [Notes](#notes)
8. [License](#license)

## Files
- **Main Scripts**:
  - `design.py`: Main script for protein design
  - `designparallel.py`: Parallel version of the design script
  - `designfromrelax.py`: Design script starting from a relaxed structure
  - `fastdesign.py`: Script for fast design (needs refinement)
- **Utility Scripts**:
  - `pymol-run.py` and `pymolscp.py`: PyMOL script writers
- **Job Management**:
  - `parallel.sh`: Shell script for running parallel designs
  - `submit.sh` and `submitp.sh`: Job submission scripts
- **Input Files**:
  - `1BTL.pdb`: Original PDB file
  - `1BTL.clean.pdb`: Cleaned PDB file
  - `relaxed_1BTL.pdb`: Relaxed structure PDB file

## Dependencies
- PyRosetta (2023 version)
- Python 3.x
- SLURM workload manager (for parallel execution)

## Installation
1. Clone the repository:
   ```
   git clone https://github.com/aakashxtha/BTLDesign.git
   cd BTLDesign
   ```
2. Ensure PyRosetta is installed in your environment. The scripts use the PyRosetta environment located at `/packages/envs/pyrosetta-2023/` on the SOL computer.

   (Optional) For command line runs only:
   ```
   module load mamba/latest
   source activate pyrosetta-2023
   ```

## Usage

### Basic Design
Command line:
```
python design.py -t <target_residue> -f <pdb_filename>
```
- `<target_residue>`: Can be in PDB numbering (e.g., `44A`) or pose numbering (e.g., `19`)
- `<pdb_filename>`: Input PDB file

Example:
```
python design.py -t 44A -f 1BTL.pdb
```

Using SLURM:
```
./parallel.sh 1BTL <target_residue> 1
```
Estimated runtime: 1-2 hours on a standard desktop computer. Generates 10 design structures and a score file (.fasc).

### Parallel Design
Using SLURM:
```
./parallel.sh <pdb_file> <target_residue> <number_of_runs>
```
Example:
```
./parallel.sh 1BTL 44A 10
```
This submits 10 parallel jobs, each running `designparallel.py`. Each run generates 10 design structures, totaling 100 designs.
Estimated runtime: 2-3 hours on HPC (configured for SOL supercomputer).

### Design from Relaxed Structure
```
python designfromrelax.py -t <target_residue>
```
or
```
sbatch submit.sh <target_residue>
```

## Output
Results are generated in the `Outputs/mutating_res<target>` directory:
- `1BTL_design_<design_num>.pdb` or `1BTL_design_<run_num>_<design_num>.pdb`: Designed PDB structures
- `1BTL_design.fasc` or `1BTL_design_<run_num>.fasc`: Design scores log
- `mutations.txt`: Log of mutations made in each design
- `radius.txt`: Log of design and repack shell radii used

## Design Process
1. Clean the input PDB
2. Relax the structure
3. Set up design parameters (design shell, repack shell, etc.)
4. Perform design around the target residue
5. Minimize and relax the designed structure
6. Output the results

## Notes
- Default approach uses a fixed backbone
- Certain important residues of 1BTL are restricted to repacking only
- Design shell radius is randomly chosen between 8-12Å for each iteration
- Repack shell extends 4Å beyond the design shell

## License
See the `LICENSE` file for details.

---

For more detailed information about each script and its specific functionality, please refer to the individual script files.
