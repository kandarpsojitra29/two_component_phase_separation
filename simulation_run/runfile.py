
# coding: utf-8

"""
For scenatio where interaction between A and B monomer is 0.85 and between A monomers is 1. Interaction between B monomers is varied from 0 to 1.

User inputs in order 
1) HOOMD input file gsd
"""
import sys,os,numpy as np
import hoomd, hoomd.md as md
from hoomd import azplugins
import gsd, gsd.hoomd, gsd.pygsd

prefactor={}

eps= np.arange(0,1.01,0.02)                  ##############lambda values over which scanning is done

for j in range(len(eps)):
    
    hoomd.context.initialize()
    
    if j==0:
        system=hoomd.init.read_gsd(sys.argv[1])
    else:
        system=hoomd.init.read_gsd('last_%s.gsd'%eps[j-1])             ##########Load previous configuration
    
    nl=hoomd.md.nlist.cell()

    long_period = 10000
    long_length = 1e7

    harmonic=hoomd.md.bond.harmonic()
    harmonic.bond_coeff.set('AA_bond',k=8368,r0=0.38)

    nl.reset_exclusions(exclusions=['1-2','body'])
    nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)

    yukawa=hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
    yukawa.pair_coeff.set(['A','B'],['A','B'],epsilon=0.0,kappa=1.0)

    all=hoomd.group.all()

    hoomd.md.integrate.mode_standard(dt=0.01)
    kTinput = 300*8.3144598/1000
    integrator=hoomd.md.integrate.langevin(group=all,kT=kTinput,seed=63535)
    integrator.set_gamma('A', gamma=0.1)
    integrator.set_gamma('B', gamma=0.1)

    nb.pair_coeff.set('A','A',lam=1.00,epsilon=0.8368,sigma=0.5,r_cut=2.0)
    nb.pair_coeff.set('A','B',lam=0.85,epsilon=0.8368,sigma=0.5,r_cut=2.0)
    nb.pair_coeff.set('B','B',lam=eps[j],epsilon=0.8368,sigma=0.5,r_cut=2.0)
    
	######### Equilibrate the system for 500ns at the first data point############
    if j ==0:
        hoomd.run(5e7)
    
    logger= hoomd.analyze.log("run_BB_%.3f.log"%(eps[j]),quantities=['potential_energy','pair_ashbaugh_energy','pair_yukawa_energy','bond_harmonic_energy','pressure_xx','pressure_yy','pressure_zz','temperature','lx','ly','lz'],period=1000, overwrite=True,header_prefix='#')
    trj= hoomd.dump.gsd("test_BB_%.3f.gsd"%(eps[j]),period=50000,dynamic=['property','momentum'],group=all,overwrite=True)
    writegsd=hoomd.dump.gsd('restart_BB_%.3f.gsd'%(eps[j]),period=10000000, group=all,overwrite=True,truncate=True,dynamic=['property', 'momentum'])
    hoomd.run(long_length)
    writegsd.disable()
    trj.disable()
    logger.disable()

	##### Get last frame configuration to use it for next run
    lastframe = gsd.hoomd.HOOMDTrajectory(gsd.pygsd.GSDFile(open("test_BB_%.3f.gsd"%(eps[j]),'rb')))[-1]
    lastframe.configuration.step = 0
    newgsdfile=gsd.hoomd.open('last_%s.gsd'%(eps[j]),'wb')
    newgsdfile.append(lastframe)