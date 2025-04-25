import sys
import numpy as np
import pickle as pkl
import mdtraj as md

AB_resis = np.arange(0, 137)
AB_shorter_resis = np.arange(0, 127)
DBD_hinge_LBD_resis = np.arange(137, 505)
DBD_resis = np.arange(137, 205)
LBD_resis = np.arange(234, 505)
D33 = np.array([32])

contacts_cutoff=1.1

def save_pickle(filename, pickle_obj):
    with open(filename, 'wb') as f:
        pkl.dump(pickle_obj, f)

def load_pickle(filename):
    with open(filename, 'rb') as f:
        loaded_obj = pkl.load(f)
        
    return loaded_obj

#Make array of all resi combinations for input to mdtraj distance/contacts calcuations
def make_resi_pairs(resi_range_1, resi_range_2):
    resis_list_1 = np.concatenate([resi_range_1]*len(resi_range_2),axis=None)
    resis_list_2 = np.repeat(resi_range_2, len(resi_range_1), axis=None)
    resi_pairs = np.append(resis_list_1.reshape(-1,1), resis_list_2.reshape(-1,1), axis=1)
    
    return resi_pairs

#Function to calculate contacts for each resi in group1 to all residues in group2
#For each protein returns a list along the sequence of group1 with average nr contacts to group2
def get_seq_contacts(group1_resis, group2_resis, contacts_cutoff, traj):

    #Calculate contacts between each residue in sequence and group
    seq_contacts = []
    for resi in group1_resis:

        #Get resi pairs
        resi_pairs = make_resi_pairs([resi], group2_resis)  
        #Calculate distances
        resi_distances = md.compute_contacts(traj, contacts=resi_pairs, scheme='closest')
        #Count contacts
        resi_contacts = np.count_nonzero(resi_distances[0]<contacts_cutoff, axis=1)
        #Average contacts over simulation and append to list
        seq_contacts.append(resi_contacts)

        print(f'Calculated contacts for resi {resi}')

    return np.array(seq_contacts)

traj = md.load('prodrun_joinedreplicas.xtc', top='free_rep1/PRO_CG.gro')

backbone_atom_indeces = traj.top.select('name BB')
traj = traj.atom_slice(backbone_atom_indeces)

#A/B contacts on LBD-hinge-DBD sequence
DBD_hinge_LBD_seq_AB_contacts = get_seq_contacts(DBD_hinge_LBD_resis, AB_resis, contacts_cutoff, traj)
save_pickle('pickles/DBD_hinge_LBD_seq_AB_contacts_CA.pkl', DBD_hinge_LBD_seq_AB_contacts)

#A/B (shortened) contacts on LBD-hinge-DBD sequence
DBD_hinge_LBD_seq_AB_shorter_contacts = get_seq_contacts(DBD_hinge_LBD_resis, AB_shorter_resis, contacts_cutoff, traj)
save_pickle('pickles/DBD_hinge_LBD_seq_AB_shorter_contacts_CA.pkl', DBD_hinge_LBD_seq_AB_shorter_contacts)

#DBD contacts on A/B sequence
AB_seq_DBD_contacts = get_seq_contacts(AB_resis, DBD_resis, contacts_cutoff, traj)
save_pickle('pickles/AB_seq_DBD_contacts_CA.pkl', AB_seq_DBD_contacts)

#LBD contacts on A/B sequence
AB_seq_LBD_contacts = get_seq_contacts(AB_resis, LBD_resis, contacts_cutoff, traj)
save_pickle('pickles/AB_seq_LBD_contacts_CA.pkl', AB_seq_LBD_contacts)

#DBD-hinge-LBD contacts on A/B sequence
AB_seq_DBD_hinge_LBD_contacts = get_seq_contacts(AB_resis, DBD_hinge_LBD_resis, contacts_cutoff, traj)
save_pickle('pickles/AB_seq_DBD_hinge_LBD_contacts_CA.pkl', AB_seq_DBD_hinge_LBD_contacts)

#D33 contacts to LBD sequence
LBD_seq_D33_contacts = get_seq_contacts(LBD_resis, D33, contacts_cutoff, traj)
save_pickle('pickles/LBD_seq_D33_contacts_CA.pkl', LBD_seq_D33_contacts)
