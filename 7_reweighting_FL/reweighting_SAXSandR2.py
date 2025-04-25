import numpy as np
import sys

bme_dir='/lindorffgrp-isilon/thomasen/software/BME'
sys.path.append(bme_dir)
import BME as BME
import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

R2_fmod_err_rel = float(sys.argv[1])
thetas = np.geomspace(10,10000,10)
scale_offset_iterations = 15

SAXS_exp_file = '../PPARg_lowsalt_SAXSexpt.dat'
SAXS_calc_file = '../calc_SAXS_free_joined.dat'
R2_exp_file = '../AB_seq_DBD_hinge_LBD_contacts_CA_allresis_R2diff_BMEexpt.dat'
R2_calc_file = '../AB_seq_DBD_hinge_LBD_contacts_CA_allresis_R2diff_BMEcalc.dat'
initial_weights_file = '../../6_reweighting_dAB/weights/weights.pkl'

def save_pickle(filename, pickle_obj):
    with open(filename, 'wb') as f:
        pkl.dump(pickle_obj, f)

def load_pickle(filename):
    with open(filename, 'rb') as f:
        loaded_obj = pkl.load(f)
        
    return loaded_obj

def fit_scale_offset(sim, exp, exp_err):
    #Get weight for each point based on exp error
    sample_weight=1.0/(exp_err**2)
    
    #Linear regression
    reg = LinearRegression(fit_intercept=True).fit(sim.reshape(-1,1),exp.reshape(-1,1),sample_weight=sample_weight)
    r_value = reg.score(sim.reshape(-1,1),exp.reshape(-1,1),sample_weight=sample_weight)
    slope,intercept = reg.coef_[0],reg.intercept_
    
    sim_fit = sim*slope+intercept
    
    return sim_fit, slope, intercept, r_value

def load_BME_files(exp_file, calc_file):
    df_exp = pd.read_csv(exp_file,sep="\s+",header=None,comment="#")
    df_calc = pd.read_csv(calc_file,sep="\s+",header=None,comment="#")
    
    return df_exp, df_calc

SAXS_load = BME.Reweight('SAXS_load')
SAXS_load.load(SAXS_exp_file,SAXS_calc_file)
SAXS_exp = SAXS_load.get_experiment()
SAXS_calc = SAXS_load.get_calculated()

R2_load = BME.Reweight('R2_load')
R2_load.load(R2_exp_file,R2_calc_file)
R2_exp = R2_load.get_experiment()
R2_calc = R2_load.get_calculated()

initial_weights = load_pickle(initial_weights_file)

chi2_vs_theta = []
phi_vs_theta = []
weights_vs_theta = []
chi2_vs_iterations_vs_theta = []

print(R2_exp.shape)
print(R2_exp[:,0])
print(R2_exp[:,1])

#Replace deltaR2 error with propagated exp err and forward model err
R2_fmod_err = R2_exp[:,0]*R2_fmod_err_rel
R2_exp[:,1] = np.sqrt(np.square(R2_exp[:,1]) + np.square(R2_fmod_err))

for theta in thetas:
    chi2_vs_iteration = []
    weights=initial_weights    
    for i in range(scale_offset_iterations):
        SAXS_calc_avg = np.average(SAXS_calc, axis=0, weights=weights)
        SAXS_calc_avg_fit, SAXS_scale, SAXS_offset, r_value = fit_scale_offset(SAXS_calc_avg, SAXS_exp[...,0], SAXS_exp[...,1])
        SAXS_calc_scale_offset = SAXS_calc*SAXS_scale + SAXS_offset

        R2_calc_avg = np.average(R2_calc, axis=0, weights=weights)
        R2_calc_avg_fit, R2_scale, R2_offset, r_value = fit_scale_offset(R2_calc_avg, R2_exp[...,0], R2_exp[...,1])
        R2_calc_scale_offset = R2_calc*R2_scale + R2_offset
        
        # initialize. A name must be specified
        rew = BME.Reweight(f'BME_theta{theta}', w0=initial_weights)

        # load the experimental and calculated datasets
        rew.load_array('SAXS', SAXS_exp, SAXS_calc_scale_offset)
        rew.load_array('SAXS', R2_exp, R2_calc_scale_offset)

        chi2_before, chi2_after, phi = rew.fit(theta=theta)
        weights = rew.get_weights()

        chi2_vs_iteration.append(chi2_after)
        
    chi2_vs_theta.append(chi2_after)
    phi_vs_theta.append(phi)
    weights_vs_theta.append(weights)
    chi2_vs_iterations_vs_theta.append(chi2_vs_iteration)

reweighting = {}
reweighting['chi2'] = chi2_vs_theta
reweighting['phi'] = phi_vs_theta
reweighting['weights'] = weights_vs_theta
reweighting['theta'] = thetas
reweighting['chi2_vs_scaleoffsetiterations'] = chi2_vs_iterations_vs_theta
reweighting['R2_fmod_err'] = R2_fmod_err

save_pickle('reweighting.pkl', reweighting)
