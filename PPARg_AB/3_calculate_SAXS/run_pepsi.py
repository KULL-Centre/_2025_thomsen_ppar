import os, sys
import numpy as np
import mdtraj as md

pepsi_path = '/storage1/thomasen/software/Pepsi-SAXS'
drho = 1.0
r0 = 1.025
exp_SAXS = str(sys.argv[1])
AA_traj = str(sys.argv[2])
AA_top = str(sys.argv[3])
outfile_calc = str(sys.argv[4])

traj = md.load(AA_traj,top=AA_top)

#Run Pepsi with fixed average params for each frame and get chi2
for i in range(len(traj)):
    
    #Save traj frame
    traj[i].save_pdb('AA_frame.pdb')
    
    #Calculate SAXS profile with pepsi
    outfile = f'SAXS_frame{i}_fixparams.fit'
    command = pepsi_path + ' AA_frame.pdb ' + exp_SAXS + ' -o ' + outfile + ' --I0 1.0 --scaleFactor 1.0 --dro ' + str(drho) + ' --r0_min_factor ' + str(r0) + ' --r0_max_factor ' + str(r0) + ' --r0_N 1'
    os.system(command)
    
    #Read SAXS profile and add to running sum
    q,dI,Ifit = np.genfromtxt(outfile,skip_header=6,skip_footer=0,usecols=[0,2,3],unpack=True)
    
    #Write header and q-values to BME calc file on first iteration
    if i==0:
        header = "# label"
        for q_value in q:
            header += " \t %e" % q_value
        header += " \n"
        
        with open(outfile_calc,'w') as f:
            f.write(header)
    
    #Write SAXS profile to BME calc file
    frame_line = f"frame_{i+1}"
    for Ifit_value in Ifit:
        frame_line += " \t %e" % Ifit_value
    frame_line += '\n'
    with open(outfile_calc,'a') as f:
        f.write(frame_line)

    #Clean up files
    os.system(f'rm AA_frame.pdb SAXS_frame{i}_fixparams.*')

