# Social Microbiome Simulation

This Python script simulates the dynamics of microbial colonization and extinction within a social network, where nodes represent hosts and edges represent social connections. It models social transmission of microbes, environmental acquisition, and microbial community changes over time.

## Features

- Initializes a random social network based on input parameters like the number of nodes and mean degree.
- Simulates microbial colonization, extinction, and social transmission among hosts in the network.
- Tracks changes in microbial composition and social interactions over time.
- Outputs the simulation results, including network structures, microbial communities, extinction/colonization rates, and species richness.

## Prerequisites

Ensure you have the following Python libraries installed:

- `numpy`
- `matplotlib`
- `argparse`
- `pickle`

You can install these dependencies using pip:

```bash
pip install numpy matplotlib
```

## Usage

The script can be run from the command line or within a Python environment. Below is an example of how to run it with arguments:

```bash
python social_microbiome.py <NOD> <DEG> <PN> <PR> <epsilon> <MICROBIAL_POOL_SIZE> <path_to_save>
```

### Parameters:

- `NOD`: Number of nodes (hosts) in the network.
- `DEG`: Mean degree of the network.
- `PN`: Probability of a inherited connections.
- `PR`: Probability of random connections.
- `epsilon`: Environmental acquisition rate.
- `MICROBIAL_POOL_SIZE`: The number of different microbes that can be assigned to hosts.
- `path_to_save`: Path to save the simulation results as `.npz` files.

### Example:

```bash
python social_microbiome.py 50 4 0.2 0.01 0.1 100 /path/to/output
```

### Output:

The simulation results are saved in `.npz` format and include:

- **Network adjacency matrices** over time.
- **Microbiome composition** for each host.
- **Extinction rates** of microbes for each host.
- **Colonization rates** of microbes for each host.
- **Species richness** across the network.

The output file is saved as:

```
simulation_output_<NOD>_<DEG>_<PN>_<PR>_<MICROBIAL_POOL_SIZE>.npz
```

## How It Works

1. **Social Network Initialization:**  
   A network of hosts is initialized, where edges represent potential microbial transmission pathways between hosts.
   
2. **Microbial Colonization and Extinction:**  
   Each host's microbiome is simulated based on colonization and extinction rates. Hosts can acquire microbes from their social contacts or the environment.
   
3. **Network Dynamics:**  
   Over time, the network evolves as hosts change their microbial communities through social transmission and environmental acquisition.

4. **Saving Results:**  
   The simulation outputs are saved in `.npz` format for later analysis or visualization.

## Plotting

Visualizations:
- Generates line plots to visualize the relationships between Number of Hosts, DEG (Degree), and diversity metrics such as host microbial gamma, beta, alpha diversity, and extinction and colonization rates).
- Creates FacetGrid plots to display how diversity metrics change across different PN (probability of inherited connection) and PR (probability of random) values.

---

### Additional Notes

- The values of `PN` and `PR` are crucial in determining the network's structure and microbial transmission patterns.
- The script allows you to control the frequency of environmental acquisition using the `epsilon` parameter.
  
## Contact

For any questions or issues, feel free to reach out to [Eeman Abbasi](eabbasi@sas.upenn.edu).
