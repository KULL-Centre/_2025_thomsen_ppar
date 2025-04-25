import mdtraj as md
import numpy as np
import sys

traj_type = str(sys.argv[1])
Zn_domain_traj = np.arange(2116, 3123) #Atoms for Zn domain in backmapped traj
Zn_domain_xtal = np.arange(63, 1070) #Atoms for Zn domain in xtal structure
output_traj = f'{traj_type}/prodrun_AAbackmapped_withZn.xtc'
output_top = f'{traj_type}/prodrun_AAbackmapped_withZn.gro'
xtal_pdb = 'PPARg_3E00_withZn_withH.pdb'
timestep=1000

traj = md.load(f'{traj_type}/prodrun_AAbackmapped.xtc', top=f'{traj_type}/prodrun_AAbackmapped.gro')

#Load xtal structure with Zn
xtal = md.load(xtal_pdb)

#Load xtal structure with Zn
Zn_xtal = xtal.top.select('name ZN')

for frame in range(len(traj)):
    
    print(frame)
    
    #Superpose xtal structure to traj frame
    traj_frame = traj[frame]
    aligned_xtal = md.Trajectory.superpose(xtal, traj_frame, frame=0, atom_indices=Zn_domain_xtal, ref_atom_indices=Zn_domain_traj)

    #Slice out the Zn atoms from xtal structure and add to traj frame
    aligned_Zn = aligned_xtal.atom_slice(Zn_xtal)
    traj_frame_withZn = md.Trajectory.stack(traj_frame, aligned_Zn)
    
    #Update time step
    traj_frame_withZn.time = frame*timestep
    
    #Add to output traj
    if frame==0:
        new_traj = traj_frame_withZn
    else:
        new_traj = md.Trajectory.join(new_traj, traj_frame_withZn)
        
md.Trajectory.save_xtc(new_traj, output_traj)
md.Trajectory.save_gro(new_traj[0], output_top)
