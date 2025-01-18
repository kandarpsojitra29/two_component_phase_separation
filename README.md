# two_component_phase_separation

This repository contains scripts for generating initial configuration and running simulations of minimialistic polymer model for two component phase separation

![amov-ezgif com-video-to-gif-converter](https://github.com/user-attachments/assets/2fb83f84-16f8-4350-8dc4-b2b38827edd7)

## Generating initial configuration
### 1) Create the initial configuration using generating_lattice.py

''' 
python generating_lattice.py
'''

Required inputs: seq.seq (contains sequence of two homopolymer); stats_module.dat (contains mass, charge, sigma, and lambda for each monomer)

Output: temp.xyz (Coordinates file); start.gsd (Initial configuration before resizing) restart_tmp.gsd (output trajectory from resizing)

This script will first create a lattice and then assign monomer A and B to the lattice based on the sequence provided in seq.seq. The code also defines the chain parameters (from stats_module.dat) and resize the simulation box to a cubic box of size 10nm.
