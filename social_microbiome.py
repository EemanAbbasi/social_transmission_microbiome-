#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import pickle


def calculate_pn_pr(k, N):
    if N <= k + 1:
        raise ValueError("Number of nodes N must be greater than mean degree k + 1 for valid network.")
    bmin = (N - k - 1) / (k * (N - 2))
    pn = 1 - bmin
    pr = ((k * (N - 1 - pn * (N - 2))) - N + 1) / (N - 2) * (N - 1 - k)
    pr = max(0, min(1, pr))
    
    return pn, pr

    
def pbpnprDyn(ADJM, PN, PR, NOD, DB=None):
    if DB is None:
        DB = np.random.choice(NOD, 2, replace=False)
    inherit = ADJM[:, DB[1]] * (1 - (np.random.rand(NOD) - PN > 0))
    randconn = (1 - ADJM[:, DB[1]]) * (1 - (np.random.rand(NOD) - PR > 0))
    newconn = inherit + randconn
    newconn[DB[1]] = 1
    ADJM[:, DB[0]] = newconn
    ADJM[DB[0], :] = newconn
    ADJM[DB[0], DB[0]] = 0
    return ADJM

def initGraph(NOD, DEG, PN, PR):
    P = (DEG * NOD) / (NOD ** 2)
    ADJM = np.zeros((NOD, NOD))
    lower_tri_indices = np.tril_indices(NOD, -1)
    npairs = len(lower_tri_indices[1])
    ADJM[lower_tri_indices] = (np.random.rand(npairs) <= P).astype(int)
    ADJM = ADJM + ADJM.T
    for i in range(10 * NOD):
        ADJM = pbpnprDyn(ADJM, PN, PR, NOD)
    return ADJM

def stable_mat(adjm,rounds,sampling_interval,nod,PN,PR):
    adjm_list = []
    all_hosts = np.arange(nod)
    for times in range(rounds):
        newBorn = np.random.choice(all_hosts, 1)[0]
        parents = np.random.choice(all_hosts[all_hosts != newBorn], 1)[0]
        adjm = pbpnprDyn(adjm, PN, PR, nod, [newBorn, parents])
    for times in range(1200):
        if times % sampling_interval == 0:
            adjm_list.append(adjm)
    return adjm_list

def initial_microbe_assignment(nod, microbial_pool_size): #capped at smax 
    microbiome_community = np.zeros((nod, microbial_pool_size))
    for i in range(nod):
        num_microbes = np.random.choice(np.arange(1, microbial_pool_size + 1), 1)
        microbes_indices = np.random.choice(np.arange(microbial_pool_size), num_microbes, replace=False)
        microbiome_community[i, microbes_indices] = 1
    return microbiome_community


def main():
    
    # Exclude the script name
    print("sys.argv:", sys.argv)
    
    script_name, NOD, DEG, PN, PR, epsilon, MICROBIAL_POOL_SIZE, path_to_save = sys.argv

    # Convert arguments to the correct type

    nod = int(NOD)
    DEG = int(DEG)
    PN = float(PN)
    PR = float (PR)
    epsilon = float(epsilon)
    microbial_pool_size = int(MICROBIAL_POOL_SIZE)
    path_to_save = path_to_save
    
    #PN,PR = calculate_pn_pr(DEG,nod)

    all_hosts = np.arange(nod) 
    extinction_rate = np.zeros(nod)
    colonization_rate = np.full(nod, 0.5)
    p = 0.1
    smax = 50
    #Indivdual carrying capacity = 50 (fraction of the microbial pool size)
    learning_attempts = 1000
    PE = 0.3
    timesteps = 100
    num_gen = 20 
    sampling_interval = 60
    mutation_probability = 0.001
    #epsilon = 0.1


   # filebase = args.file
    
    adjm = initGraph(nod, DEG, PN, PR)
    adjm_list = stable_mat(adjm, num_gen*timesteps,sampling_interval, nod,PN,PR)
    
    print(len(adjm_list))
    
    net_microbiome_final = []
    extinction_adjm_final = []
    colnization_adjm_final = []
    species_richness_adjm_final = []
    
    

    adjm_counter = 0
    for adjm in adjm_list:
        adjm_counter = adjm_counter + 1
    
        extinction_rate_history = []
        colonization_rate_history = []
        species_richness = []
        
    
        net_microbiome = initial_microbe_assignment(nod, smax)
        for host in range(nod):
            #Update the colonization and extinction rate based on initial microbiome composition
            colonization_rate[host] = (1 - p) ** (np.sum(net_microbiome[host, :]) - 1)
            extinction_rate[host] = ((1 - p) ** (smax - np.sum(net_microbiome[host, :])))
        
    
        for times in range(nod * 200):
            extinction_rate_history.append(np.copy(extinction_rate))
            colonization_rate_history.append(np.copy(colonization_rate))
            species_richness.append(np.sum(net_microbiome, axis=0))
            
            num_hosts = adjm.shape[0]
            host_degrees = np.sum(adjm, axis=1)
            host_probabilities = host_degrees / np.sum(host_degrees)
            
            host = np.random.choice(num_hosts, p=host_probabilities)
            
            if np.sum(net_microbiome[host, :]) == 0:
                colonization_rate[host] = 1
                extinction_rate[host] = 0
                
            else:
                colonization_rate[host] = (1 - p) ** (np.sum(net_microbiome[host, :]) - 1)
                extinction_rate[host] = ((1 - p) ** (smax - np.sum(net_microbiome[host, :])))
                
            
            # social transmission 
            
            network_hosts = np.where(adjm[host] == 1)[0]
            
            if len(network_hosts) > 0:
                nh = np.random.choice(network_hosts)
                present_microbes = np.where(net_microbiome[nh] == 1)[0]
                
                if len(present_microbes) > 0:
                    microbe_to_colonize = np.random.choice(present_microbes)
                    
                    if np.random.rand() < colonization_rate[host]:
                        net_microbiome[host, microbe_to_colonize] = 1
            
            
            present_microbes = np.where(net_microbiome[host] == 1)[0]
            if len(present_microbes) > 0:
                microbe_to_extinct = np.random.choice(present_microbes)
                if np.random.rand() < extinction_rate[host]:
                    net_microbiome[host, microbe_to_extinct] = 0
                    
                    
                    
            if np.sum(net_microbiome[host, :]) == 0:
                colonization_rate[host] = 1
                extinction_rate[host] = 0
                
            else:
                colonization_rate[host] = (1 - p) ** (np.sum(net_microbiome[host, :]) - 1)
                extinction_rate[host] = ((1 - p) ** (smax - np.sum(net_microbiome[host, :])))
                
                    
            
            host = np.random.choice(num_hosts)
            
            if np.random.rand() < PE:
                
                #sum up the abundance of the microbes - weighted aquisition 
                #environmnetal aqusition episoln * extinction and colonization rate 
                #1/microbal_pool_size add that to everything - abundance list   
                num_microbes = net_microbiome.shape[1]
                microbial_abundance = np.sum(net_microbiome, axis=0)
                microbial_abundance = microbial_abundance + 0.001
                microbe_prob = microbial_abundance / np.sum(microbial_abundance)
                env_microbe_index = np.random.choice(num_microbes, p=microbe_prob)
                
                if np.sum(net_microbiome[host, :]) == 0:
                    colonization_rate[host] = 1
                    extinction_rate[host] = 0
                    
                else:
                    colonization_rate[host] = (1 - p) ** (np.sum(net_microbiome[host, :]) - 1) * epsilon
                    extinction_rate[host] = ((1 - p) ** (smax - np.sum(net_microbiome[host, :]))) * epsilon
                
                
                if np.random.rand() < colonization_rate[host]:
                    net_microbiome[host, env_microbe_index] = 1
                
                present_microbes = np.where(net_microbiome[host] == 1)[0]
                if len(present_microbes) > 0:
                    microbe_to_extinct = np.random.choice(present_microbes)
                    if np.random.rand() < extinction_rate[host]:
                        net_microbiome[host, microbe_to_extinct] = 0
                        
                        
        
        #save the final runs
        net_microbiome_final.append(net_microbiome)
        extinction_adjm_final.append(extinction_rate_history)
        colnization_adjm_final.append(colonization_rate_history)
        species_richness_adjm_final.append(species_richness)
    

            
        
    save_path = f'{path_to_save}/simulation_output_{nod}_{DEG}_{PN}_{PR}_{MICROBIAL_POOL_SIZE}.npz'
    params = [NOD, DEG, PN, PR, MICROBIAL_POOL_SIZE]
    np.savez(save_path, params= params, adjm=adjm_list,net_microbiome=net_microbiome_final,
            extinction_adjm=extinction_adjm_final,colonization_adjm=colnization_adjm_final,
            species_richness=species_richness_adjm_final)
    
    
if __name__ == "__main__":
     main()   
 
 




import sys
sys.argv = [
    'social_microbiome.py',  # Script name (or any placeholder name)
    '50',    # NOD
    '4',     # DEG
    '0.2',   # PN
    '0.01',  # PR
    '0.1',   # epsilon
    '100',   # MICROBIAL_POOL_SIZE
    '/Users/eemanabbasi/Desktop/Paper_3/output/'  # path_to_save
]

