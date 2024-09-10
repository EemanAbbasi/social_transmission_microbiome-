#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:15:28 2024

@author: eemanabbasi
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 07:46:01 2024

@author: eemanabbasi
"""
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:40:55 2024

@author: eemanabbasi
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.spatial.distance import jaccard
import pickle
import os

import re 

path_to_save = '/Users/eemanabbasi/Desktop/Paper_3/output/'   


#Things to do:
#correlation between degree of host and the within host diversity (species richness) of that host  
#correlation between degree of host and the average beta diversity between that host and its connections 
#correlation between degree of host and the average beta diversity between that host and random individuals in the population 
#gamma diversity
#we want to calculate these across differen pn and pr values
#and also changing pn and pr while keeping mean degree constant
#20 replicate networks per pn and pr value
#network size = 50, 100, 200
#average degree = 3, 5, 10, 20, 50



def extract_parameters_from_filename(filename):
    # Regex pattern to match the filename format
    pattern = r'output_(\d+)_(\d+)_([\d\.eE\-]+)_([\d\.eE\-]+)_(\d+)\.npz'
    match = re.search(pattern, filename)
    if match:
        return {
            'NOD': int(match.group(1)),
            'DEG': int(match.group(2)),
            'PN': float(match.group(3)),
            'PR': float(match.group(4)),
            'MICROBIAL_POOL_SIZE': int(match.group(5))
        }
    else:
        raise ValueError(f"Filename does not match expected format: {filename}")



def calculate_gamma_diversity(microbiome_community):
    return np.sum(microbiome_community, axis=1)

def calculate_beta_diversity(microbiome_community):
    n_hosts = microbiome_community.shape[0]
    beta_diversity_matrix = np.zeros((n_hosts, n_hosts))
    for i in range(n_hosts):
        for j in range(i + 1, n_hosts):
            set_i = set(np.where(microbiome_community[i] == 1)[0])
            set_j = set(np.where(microbiome_community[j] == 1)[0])
            intersection = len(set_i.intersection(set_j))
            union = len(set_i.union(set_j))
            if union == 0:
                jaccard_distance = 1.0
            else:
                jaccard_distance = 1 - intersection / union
            beta_diversity_matrix[i, j] = jaccard_distance
            beta_diversity_matrix[j, i] = jaccard_distance
    return beta_diversity_matrix

def calculate_mean_pairwise_beta_diversity(beta_diversity_matrix):
    n_hosts = beta_diversity_matrix.shape[0]
    mean_pairwise_beta_diversity = np.zeros(n_hosts)
    for i in range(n_hosts):
        mean_pairwise_beta_diversity[i] = np.mean(beta_diversity_matrix[i, :])
    return mean_pairwise_beta_diversity


def calculate_beta_diversity_connected_hosts(adjacency_matrix, microbiome_matrix):
    n_hosts = adjacency_matrix.shape[0]
    beta_diversity_matrix = np.zeros((n_hosts, n_hosts))
    
    # Iterate over pairs of hosts
    for i in range(n_hosts):
        for j in range(i + 1, n_hosts):
            if adjacency_matrix[i, j] == 1:  # Check if hosts i and j are connected
                jaccard_distance = 1 - jaccard(microbiome_matrix[i], microbiome_matrix[j])
                beta_diversity_matrix[i, j] = jaccard_distance
                beta_diversity_matrix[j, i] = jaccard_distance
    
    return beta_diversity_matrix


def calculate_degree_centrality(adjacency_matrix):
    degree = np.sum(adjacency_matrix, axis=1)
    return degree


def calculate_averages(adjm_list, net_microbiome_final, extinction_adjm_final, colonization_adjm_final, species_richness_adjm_final):
    cumulative_extinction_rates = []
    cumulative_colonization_rates = []
    cumulative_species_richness = []
    cumulative_gamma_diversity = []
    cumulative_beta_diversity = []
    cumulative_beta_diversity_connected_hosts = []

    num_adjm = len(adjm_list)

    for i in range(num_adjm):
        cumulative_extinction_rates.append(extinction_adjm_final[i])
        cumulative_colonization_rates.append(colonization_adjm_final[i])
        cumulative_species_richness.append(species_richness_adjm_final[i])
        gamma_diversity = calculate_gamma_diversity(net_microbiome_final[i])
        cumulative_gamma_diversity.append(gamma_diversity)
        beta_diversity_matrix = calculate_beta_diversity(net_microbiome_final[i])
        mean_pairwise_beta_diversity = calculate_mean_pairwise_beta_diversity(beta_diversity_matrix)
        cumulative_beta_diversity.append(mean_pairwise_beta_diversity)
        beta_diversity_connected_hosts = calculate_beta_diversity_connected_hosts(adjm_list[i], net_microbiome_final[i])
        mean_pairwise_beta_diversity_connected_hosts = calculate_mean_pairwise_beta_diversity(beta_diversity_connected_hosts)
        cumulative_beta_diversity_connected_hosts.append(mean_pairwise_beta_diversity_connected_hosts)

    cumulative_extinction_rates = np.array(cumulative_extinction_rates)
    cumulative_colonization_rates = np.array(cumulative_colonization_rates)
    cumulative_species_richness = np.array(cumulative_species_richness)
    cumulative_gamma_diversity = np.array(cumulative_gamma_diversity)
    cumulative_beta_diversity = np.array(cumulative_beta_diversity)
    cumulative_beta_diversity_connected_hosts = np.array(cumulative_beta_diversity_connected_hosts)

    averages = {
        'extinction_rates': np.mean(cumulative_extinction_rates, axis=0),
        'colonization_rates': np.mean(cumulative_colonization_rates, axis=0),
        'species_richness': np.mean(cumulative_species_richness, axis=0),
        'gamma_diversity': np.mean(cumulative_gamma_diversity, axis=0),
        'beta_diversity': np.mean(cumulative_beta_diversity, axis=0),
        'beta_diversity_connected_hosts': np.mean(cumulative_beta_diversity_connected_hosts, axis=0),
        'adjm_list':adjm_list
    }

    return averages


def process_data(file_path):
    data = np.load(file_path)
    avg = calculate_averages(data['adjm'], data['net_microbiome'], data['extinction_adjm'], data['colonization_adjm'], data['species_richness'])
    return avg




def plot_results(varying_param, param_values, results, ylabel, title):

    plt.figure(figsize=(10, 6))
    
    # Check if results is a dictionary
    if isinstance(results, dict):
        for label, result in results.items():
            plt.plot(param_values, result, label=label)
    else:
        # Assuming results is a list of results if not a dictionary
        for result in results:
            plt.plot(param_values, result, label=f'{ylabel} vs {varying_param}')
    
    plt.xlabel(varying_param)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()
    plt.show()


averages_dict = {}

for filename in os.listdir(path_to_save):
    print(filename)
    
    if filename.endswith(".npz"):
        print(filename)
        file_path = os.path.join(path_to_save, filename)
        try:
            params = extract_parameters_from_filename(filename)
            averages = process_data(file_path)
            key = tuple(params.values())
            averages_dict[key] = averages
        except ValueError as e:
            print("error now")
            print(e)
            
import pandas as pd      
            
NOD_vals = []
DEG_vals = []
PN_vals = []
PR_vals = []
MICROBIAL_POOL_SIZE_vals = []
extinction_rates = []
colonization_rates = []
species_richness = []
gamma_diversity = []
beta_diversity = []
beta_diversity_connected_hosts = []

# Iterate over the dictionary to populate the lists
for key, avg in averages_dict.items():
    NOD, DEG, PN, PR, MICROBIAL_POOL_SIZE = key
    NOD_vals.append(NOD)
    DEG_vals.append(DEG)
    PN_vals.append(PN)
    PR_vals.append(PR)
    MICROBIAL_POOL_SIZE_vals.append(MICROBIAL_POOL_SIZE)
    
    extinction_rates.append(np.mean(avg['extinction_rates']))
    colonization_rates.append(np.mean(avg['colonization_rates']))
    species_richness.append(np.mean(avg['species_richness']))
    gamma_diversity.append(np.mean(avg['gamma_diversity']))
    beta_diversity.append(np.mean(avg['beta_diversity']))
    beta_diversity_connected_hosts.append(np.mean(avg['beta_diversity']))  # Assuming this is the correct data source

# Create a DataFrame with the unpacked parameters and metrics
df = pd.DataFrame({
    'NOD': NOD_vals,
    'DEG': DEG_vals,
    'PN': PN_vals,
    'PR': PR_vals,
    'MICROBIAL_POOL_SIZE': MICROBIAL_POOL_SIZE_vals,
    'extinction_rates': extinction_rates,
    'colonization_rates': colonization_rates,
    'gamma_diversity': species_richness,
    'alpha_diversity': gamma_diversity,
    'beta_diversity': beta_diversity,
    'beta_diversity_connected_hosts': beta_diversity_connected_hosts
})

    
deg_values = [3, 5, 10, 20, 50]
nod_vals = [30, 50, 100]
pr_vals = [0.01, 0.05, 0.1]
pn_vals = [0.2, 0.4, 0.8]



#import seaborn as sns
#import matplotlib.pyplot as plt


#for i in pn_vals:
   # for j in pr_vals:
        # Create a grid of heatmaps for varying DEG and NOD values
       # heatmap_data = df[(df['PN'] == i) & (df['PR'] == j)]
       # 
        # Check if the filtered DataFrame has data
        #if heatmap_data.empty:
        #    print("No data available for the specified PN and PR values.")
        #else:
            # Create a pivot table for the heatmap
           # pivot_table = heatmap_data.pivot_table(index='DEG', columns='NOD', values='gamma_diversity', aggfunc='mean')
        
            # Plot the heatmap
           # plt.figure(figsize=(10, 8))
           # sns.heatmap(pivot_table, annot=True, cmap='viridis', fmt='.2f', linewidths=0.5)
           # plt.title(f'Heatmap of alpha diversity (PN={i}, PR={j})')
            #plt.xlabel('NOD')
            #plt.ylabel('DEG')
            #plt.show()
    
    
    
import seaborn as sns
import matplotlib.pyplot as plt

# Create a FacetGrid with line plots
g = sns.FacetGrid(df, col='PN', row='PR', margin_titles=True)
g.map(sns.lineplot, 'NOD', 'gamma_diversity', 'DEG', marker='o')

g.set_axis_labels('NOD', 'Gamma Diversity')
g.add_legend()
plt.show()


# Create a FacetGrid with line plots
g = sns.FacetGrid(df, col='PN', row='PR', margin_titles=True)
g.map(sns.lineplot, 'NOD', 'alpha_diversity', 'DEG', marker='o')

g.set_axis_labels('NOD', 'Alpha Diversity')
g.add_legend()
plt.show()


# Create a FacetGrid with line plots
g = sns.FacetGrid(df, col='PN', row='PR', margin_titles=True)
g.map(sns.lineplot, 'NOD', 'extinction_rates', 'DEG', marker='o')

g.set_axis_labels('NOD', 'Extinction Rates')
g.add_legend()
plt.show()



# Create a FacetGrid with line plots
g = sns.FacetGrid(df, col='PN', row='PR', margin_titles=True)
g.map(sns.lineplot, 'NOD', 'colonization_rates', 'DEG', marker='o')

g.set_axis_labels('NOD', 'Colonization Rates')
g.add_legend()
plt.show()

g = sns.FacetGrid(df, col='PN', row='PR', margin_titles=True)
g.map(sns.lineplot, 'NOD', 'beta_diversity', 'DEG', marker='o')

g.set_axis_labels('NOD', 'Beta Diversity')
g.add_legend()
plt.show()

g = sns.FacetGrid(df, col='PN', row='PR', margin_titles=True)
g.map(sns.lineplot, 'NOD', 'beta_diversity_connected_hosts', 'DEG', marker='o')

g.set_axis_labels('NOD', 'Beta Diversity Connected Hosts')
g.add_legend()
plt.show()
