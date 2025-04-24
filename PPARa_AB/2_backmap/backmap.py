import os
import sys
import mdtraj as md

pulchra = ' /projects/prism/people/hzr104/envs/pulchra/bin/pulchra'
cg_traj = str(sys.argv[1]) #'PPARg_AB/traj.xtc'
cg_top = str(sys.argv[2]) #'PPARg_AB/top.pdb'
aa_xtc_out = str(sys.argv[3])
aa_pdb_out = str(sys.argv[4])
aa_gro_out = str(sys.argv[5])

#Load CG trajectory
traj = md.load(cg_traj, top=cg_top)


def fix_topology(t,seq):
    cgtop = md.Topology()
    cgchain = cgtop.add_chain()
    for res in seq:
        cgres = cgtop.add_residue(res, cgchain)
        cgtop.add_atom('CA', element=md.element.carbon, residue=cgres)
    traj = md.Trajectory(t.xyz, cgtop, t.time, t.unitcell_lengths, t.unitcell_angles)
    traj = traj.superpose(traj, frame=0)
    return traj

seq = []
for residue in traj.topology.residues:
    seq.append(str(residue)[0:3])

traj_fixed = fix_topology(traj,seq)

#Loop over CG trajectory frames
for i in range(len(traj)):

    #Save frame to pdb file
    md.Trajectory.save_pdb(traj_fixed[i], 'CG_frame.pdb')
    #Backmap
    os.system(f'{pulchra} CG_frame.pdb')
    #Load backmapped frame
    aa_frame = md.load('CG_frame.rebuilt.pdb', top='CG_frame.rebuilt.pdb')
    #Remove saved files for CG frame and backmapped frame
    os.system('rm CG_frame.pdb CG_frame.rebuilt.pdb')
    #Change time of backmapped frame
    aa_frame.time == i*1000
    #Add new backmapped frame to trajectory
    if i==0:
        aa_traj = aa_frame
    else:
        aa_traj = md.Trajectory.join(aa_traj, aa_frame)

md.Trajectory.save_xtc(aa_traj, aa_xtc_out)
md.Trajectory.save_pdb(aa_traj[0], aa_pdb_out)
md.Trajectory.save_gro(aa_traj[0], aa_gro_out)
