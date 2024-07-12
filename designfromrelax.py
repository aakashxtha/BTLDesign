import os
import random
import argparse
import pyrosetta

def main(target=None):
    # Initialize PyRosetta with custom flags
    pyrosetta.init('''
        -relax:default_repeats 5
        -relax:constrain_relax_to_start_coords
        -relax:coord_constrain_sidechains
        -relax:ramp_constraints false
        -score::hbond_params correct_params
        -no_his_his_pairE
        -extrachi_cutoff 1
        -multi_cool_annealer 10
        -ex1 -ex2
        -use_input_sc
        -flip_HNQ
        -ignore_unrecognized_res
        -relax:coord_cst_stdev 0.5
    ''')

    # Load the cleaned PDB and relaxed pdb
    btl = pyrosetta.pose_from_pdb("1BTL.clean.pdb")
    original_pose = btl.clone()
    relaxed_pose = pyrosetta.pose_from_pdb("relaxed_1BTL.pdb")
    scorefxn = pyrosetta.get_fa_scorefxn()
    print("Original score: ",scorefxn(original_pose))
    print("Relaxed score: ",scorefxn(relaxed_pose))

    # Setup directory structure
    main_dir = os.getcwd()
    output_dir = os.path.join("Outputs", f"mutating_res{target}")
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)

    # FastRelax 
    # fr = pyrosetta.rosetta.protocols.relax.FastRelax()
    # fr.set_scorefxn(scorefxn)
    # fr.apply(btl)
    # print(scorefxn.score(btl))
    # relaxed_pose = btl.clone()
    # btl.dump_pdb(f"relaxed_{filename}")
    
    # Setup MoveMap for flexible backbone
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)
    mm.set_jump(True)
    
    # Movers Setup
    min_bb = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    min_bb.score_function(scorefxn)
    min_bb.movemap(mm)
    min_bb.min_type('dfpmin_armijo_nonmonotone')
    min_bb.tolerance(0.001)
    min_bb.max_iter(1000)

    multi_min = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()
    multi_min.set_mover(min_bb)
    multi_min.set_scorefxn(scorefxn)
    multi_min.set_max_accepted_trials(10)
    multi_min.set_sampletype("low")
    multi_min.set_temperature(0.6)
    multi_min.set_recover_low(True)
    multi_min.set_preapply(False)
    
    # Job distributor for designs
    job = pyrosetta.toolbox.py_jobdistributor.PyJobDistributor('1BTL_design', 3, scorefxn)
    job.native_pose = original_pose
    pose = pyrosetta.Pose()

    while not job.job_complete:
        pose.assign(relaxed_pose)
        design_around(pose, target, scorefxn)# design
        i = print_mutations(original_pose, pose) #prints the mutations
        if len(i) == 0:
            continue
        multi_min.apply(pose)# minimization
        print(scorefxn.score(pose))
        # fr.apply(pose) # relax
        # print(scorefxn.score(pose))
        job.output_decoy(pose)
    os.chdir(main_dir)
    
    
def design_around(pose, target, scorefxn):
    """Setups the design around protocol and design the pose accordingly."""
    #Selectors
    target_residue = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(str(target))
    design_radius = random.randint(8, 12)
    design_shell = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(target_residue, design_radius, False)
    repack_shell = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(target_residue, design_radius + 4, True)
    # Important residues that should not be designed
    imp_residues = "2,19,20,37,42,45,48,51,82,97,105,107,120,141,142,149,158,194,201,209,218,226,230,235"
    imp_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(imp_residues)
    
    # Setup TaskFactory
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    # Task operations
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), imp_residue_selector, False))
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), design_shell, True))
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), repack_shell, True))
    
    # Mover Setup
    prm = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
    prm.task_factory(tf)
    prm.score_function(scorefxn)
    prm.apply(pose)    
    print(scorefxn.score(pose))
    print_radius(design_radius)
    
def print_mutations(original_pose, designed_pose):
    """Prints and logs the mutations between the original and designed poses."""
    original_seq = original_pose.sequence()
    designed_seq = designed_pose.sequence()
    mutations = []
    for i in range(0, original_pose.total_residue()):
        if original_seq[i] != designed_seq[i]:
            resid = original_pose.pdb_info().pose2pdb(i)
            mutations.append(f"{original_seq[i]}{resid.split()[0]}{designed_seq[i]}")
            # mutations.append(f"{original_seq[i]}{i+1}{designed_seq[i]}")
    mutations_str = "Mutations: " + ", ".join(mutations)
    print(mutations_str)
    with open("mutations.txt", "a") as mutation_file:
        if not len(mutations) == 0:
            mutation_file.write(mutations_str + "\n")
        else:
            mutation_file.write("\n")
        mutation_file.close()
    return mutations

def print_radius(radius):
    """Prints and logs the radius of design shell and repack shell."""
    radius_str = f"Design shell: {radius}A      Repack shell: {radius+4}A"
    print(radius_str)
    with open("radius.txt", "a") as radius_file:
        radius_file.write(radius_str + "\n")
        radius_file.close()            
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", type=str, help="Target residue to mutate as integer.")
    args = parser.parse_args()

    # Run protocol
    main(target=args.target)

# python designfromrelax.py -t $target
        
