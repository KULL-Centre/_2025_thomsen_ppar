import numpy as np
import sys

bme_dir='/lindorffgrp-isilon/thomasen/software/BME'
sys.path.append(bme_dir)
import BME as BME
import pickle as pkl

def save_pickle(filename, pickle_obj):
    with open(filename, 'wb') as f:
        pkl.dump(pickle_obj, f)

def load_pickle(filename):
    with open(filename, 'rb') as f:
        loaded_obj = pkl.load(f)
        
    return loaded_obj

#BME reweighting
thetas = np.geomspace(0.01,10000,20)
    
exp_file = 'PPARg_ABonly_SAXSexpt.dat'
calc_file = 'calc_SAXS.dat'

chi2_vs_theta = []
phi_vs_theta = []
weights_vs_theta = []

for theta in thetas: 
    
    print(f'theta {theta}')
    
    # initialize. A name must be specified 
    rew = BME.Reweight(f'BME_theta{theta}')

    # load the experimental and calculated datasets
    rew.load(exp_file,calc_file)

    # fit the data 
    chi2_before, chi2_after, phi, calc0, calc_rew = rew.ibme(theta=theta,iterations=30,ftol=0.001,offset=True)
    weights = rew.get_ibme_weights()
    
    chi2_vs_theta.append(chi2_after)
    phi_vs_theta.append(phi)
    weights_vs_theta.append(weights[-1])
    
reweighting = {}
reweighting['chi2'] = chi2_vs_theta
reweighting['phi'] = phi_vs_theta
reweighting['weights'] = weights_vs_theta
reweighting['theta'] = thetas

save_pickle('reweighting.pkl', reweighting)
