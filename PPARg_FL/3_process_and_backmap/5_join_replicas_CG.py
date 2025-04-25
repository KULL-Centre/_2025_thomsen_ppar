import mdtraj as md

replicas = 2

#Make full traj
timestep = 1000

frame_out = 0
for rep in range(1, replicas+1):
    traj = md.load(f'PPARg_rep{rep}/prodrun_nopbc.xtc', top=f'PPARg_rep{rep}/PRO_CG.gro')

    for i in range(len(traj)):
        frame = traj[i]
        frame.time = timestep*frame_out
        frame_out += 1

        print(frame.time)

        if rep==1 and i==0:
            traj_joined = frame
        else:
            traj_joined = traj_joined.join(frame)

traj_joined.save('prodrun_joinedreplicas.xtc')
