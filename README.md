# two_component_phase_separation

This repository contains scripts for generating the initial configuration and running simulations of a minimalistic polymer model for two-component phase separation.

![amov-ezgif com-video-to-gif-converter](https://github.com/user-attachments/assets/2fb83f84-16f8-4350-8dc4-b2b38827edd7)

## Generating initial configuration
#### 1) Create the initial configuration using `generating_lattice.py`
```
python generating_lattice.py
```
<ins>Required inputs:</ins>
- `seq.seq`: Contains the sequence of two homopolymers.
- `stats_module.dat`: Contains mass, charge, sigma, and lambda for each monomer.

<ins>Output:</ins>
- `coord.xyz`: Coordinates file.
- `start.gsd`: Initial configuration before resizing.
- `resized_config.gsd`: Output trajectory from resizing.

<ins>Description:</ins>
This script will first create a lattice and assign monomers A and B to the lattice based on the sequence provided in `seq.seq`. The code also defines the chain parameters (from `stats_module.dat`) and resizes the simulation box to a cubic box of size 10 nm.

#### 2) Obtain the last frame of the resized configuration using `get_lastframe_newgsd.py`
```
python get_lastframe_newgsd.py resized_config.gsd
```
<ins>Required input:</ins>
- `resized_config.gsd`: Output from the previous step.

<ins>Output:</ins>
- `10.094831x10.094831x10.2002_box.gsd`: Last frame of the resized configuration.

<ins>Description:</ins>
This step extracts the last frame from the resized box obtained in the previous step.

#### 3) Extend the box length in the z-direction using `box2slab_extend_newHOOMD.py`
```
python box2slab_extend_newHOOMD.py 10.094831x10.094831x10.2002_box.gsd
```
<ins>Required inputs:</ins>
- `10.094831x10.094831x10.2002_box.gsd`: Last frame GSD file obtained in the previous step.

<ins>Output:</ins>
- `box2slab_extend.gsd`: Final configuration with the extended z-dimension.

<ins>Description:</ins>
The code will ask for user inputs, including:
- Box length to be extended in the z-direction (e.g., 75 nm).
- Number of components (e.g., 2).
- Length of each component (e.g., 20).
- Number of chains of each component (e.g., 250).

The final output GSD file consists of a rectangular slab configuration.

## Running the simulation
```
python runfile.py box2slab_extend.gsd
```

<ins>Required inputs:</ins> 
- box2slab_extend.gsd: Rectangular slab configuration.

<ins>Output:</ins> 
- test_BB_0.000.gsd: Trajectory file.
- restart_BB_0.000.gsd: Restart file.
- run_BB_0.000.log: Log file.
- last_0.0.gsd: Last frame of the simulation.
