import os
import re

def parse_mutations(mutation_file):
    mutations = {}
    with open(mutation_file, 'r') as f:
        for line in f:
            parts = line.strip().split(' - Mutations: ')
            if len(parts) == 2:
                pdb_file, mutation_list = parts
                mutations[pdb_file] = mutation_list.split(', ')
    return mutations

def generate_pymol_script():
    # Get all PDB files in the current directory
    pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]
    
    # Find the clean PDB file
    clean_pdb = next((f for f in pdb_files if 'clean' in f), None)
    if not clean_pdb:
        raise ValueError("No clean PDB file found.")
    
    # Parse mutations from mutation.txt
    mutations = parse_mutations('mutations.txt')
    
    # Generate PyMOL commands
    commands = [
        "# Load and align structures",
        f"load {clean_pdb}, clean",
        "hide everything, all",
        "show cartoon, all",
        "set cartoon_transparency, 0.5",
    ]
    
    for pdb_file in pdb_files:
        if pdb_file != clean_pdb:
            design_name = os.path.splitext(pdb_file)[0]
            commands.extend([
                f"load {pdb_file}, {design_name}",
                f"align {design_name}, clean"
            ])
    
    commands.append("\n# Show mutations")
    
    # Create selections for each design's mutations
    for pdb_file, mutation_list in mutations.items():
        if pdb_file in pdb_files and pdb_file != clean_pdb:
            design_name = os.path.splitext(pdb_file)[0]
            res_nums = [re.search(r'\d+', mutation).group() for mutation in mutation_list]
            res_selection = '+'.join(res_nums)
            commands.extend([
                f"select mut_{design_name}, {design_name} and resi {res_selection}",
                f"show sticks, mut_{design_name}",
                f"set cartoon_transparency, 0, mut_{design_name}"
            ])
    
    # Set up view
    commands.extend([
        "\n# Set up view",
        "center clean",
        "zoom all",
        "set ray_shadows, 0",
        "set depth_cue, 0"
    ])
    
    # Write commands to a PyMOL script file
    output_file = 'align_and_show_all_mutations.pml'
    with open(output_file, 'w') as f:
        f.write('\n'.join(commands))
    
    print(f"PyMOL script has been generated: {output_file}")

if __name__ == "__main__":
    generate_pymol_script()