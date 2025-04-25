import os, sys

start_frame = int(sys.argv[1])
frame_skip = int(sys.argv[2])
nr_frames = 40001
pepsi_path = '/home/projects/ku_10001/apps/Pepsi-SAXS'
exp_SAXS = f'../../PPARg_lowsalt_SAXSexpt.dat'
drho = 1.0
r0 = 1.025

backmapped_dir='../AA_frames_withZn'

#Run Pepsi with fixed average params for each frame and get chi2
for i in range(start_frame,nr_frames,frame_skip):
    
    while os.path.isfile(f'SAXS_frame{i}.fit') == False:    

        #Calculate SAXS profile with pepsi
        outfile = f'SAXS_frame{i}.fit'
        command = pepsi_path + f' {backmapped_dir}/AA_frame{i}.pdb ' + exp_SAXS + ' -o ' + outfile + ' --I0 1.0 --scaleFactor 1.0 --dro ' + str(drho) + ' --r0_min_factor ' + str(r0) + ' --r0_max_factor ' + str(r0) + ' --r0_N 1'
        os.system(command)
    

