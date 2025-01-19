import sys,os,numpy as np
import hoomd, hoomd.md as md
import hoomd.deprecated as old
from hoomd import azplugins
import gsd, gsd.hoomd

#####Creating lattice###############
a1=[0.38,0,0]                         #0.38 corresponds to bond length in nm
a2=[0,1,0]
a3=[0,0,1]
a1=np.array(a1)
a2=np.array(a2)
a3=np.array(a3)

outputFile= open('coord.xyz','w')

#########Chain1
x1= 20								  #number of  monomer of chains1 in x direction
y1= 25								  #repeats of chains1 in y direction
z1= 10								  #repeats of chain1 in z direction
for k in range(z1):   											
    for j in range(y1):											
        for i in range(x1):										
            R = i*a1 + j*a2 + k*a3
            outputFile.write(str(R[0]) + ' ' + str(R[1]) + ' ' + str(R[2]) + '\n')

##########Chain2
x2= 20                                #number of  monomer of chains2 in x direction
y2= 25								  #repeats of chains2 in y direction
z2= 20   							  #repeats of chain2 in z direction                                                       
for k in range(z1,z2):   										
    for j in range(y2):											
        for i in range(x2):										
            R = i*a1 + j*a2 + k*a3
            outputFile.write(str(R[0]) + ' ' + str(R[1]) + ' ' + str(R[2]) + '\n')
outputFile.close()

####################################

#####Generating sequence for the entire system from initial sequence of chain1 and chain2
n1= 250                               #number of chain1 in the system
n2= 250                               #number of chain2 in the system
seq={'A':'A','B':'B'}
count=0
nline=1
fout=open('seq_seq3.dat','w')         #seq file for entire system
f= open('seq.seq','r')                #Load the sequence file with sequence of two polymers
fr= f.readlines()
fid1 = fr[0][0:x1]			          #Defining chain1 in the seq file
fid2=  fr[0][x1:-1]					  #Defining chain2 in the seq file
fid1 = fid1*n1													
fid2 = fid2*n2  												
fid = fid1 + fid2                                              

for i in fid:
    if i[0]!='#':
        for j in i:
            if j in seq:
                fout.write(' %s'%seq[j])
                count+=1
                if count==nline:
                    fout.write('\n')
                    count=0
fout.close()

# Input parameters of all the monomers
ff_para = 'stats_module.dat'
aalist={}
with open(ff_para,'r') as fid:
    for i in fid:
        if i[0]!='#':
            tmp=i.rsplit()
            aalist[tmp[0]]=np.loadtxt(tmp[1:],dtype=float)
aakeys=list(aalist.keys())

## Mass
aamass=[]
aacharge=[]
aaradius=[]
aahps=[]
for i in aakeys:
    aamass.append(aalist[i][0])
    aacharge.append(aalist[i][1])
    aaradius.append(aalist[i][2])
    aahps.append(aalist[i][3])

chain_id=[]
chain_mass=[]
chain_charge=[]
with open('seq_seq3.dat','r') as fid:
    for i in fid:
        iname=i.rsplit()[0]
        chain_id.append(aakeys.index(iname))
        chain_mass.append(aalist[iname][0])
        chain_charge.append(aalist[iname][1])
print('This is the sequence coded in numbers:',chain_id)

bond_length=0.38
num_total=len(chain_id)
box_length=100                                      

# Now we can build HOOMD structure
# data structure for one single frame
s=gsd.hoomd.Snapshot()
s.particles.N = num_total
# Claim the name of monomers
s.particles.types = aakeys
s.particles.typeid = chain_id
s.particles.mass = chain_mass
s.particles.charge = chain_charge
print(s.particles.typeid)
print(s.particles.types)
# Build initial position as a linear chain

pos= np.loadtxt('coord.xyz',comments=['#'])
pos=np.array(pos)
s.particles.position= pos
# Initialize bond
nbonds=(num_total-1)                                 
s.bonds.N = nbonds-(n1+n2-1)								   #No. of atoms for forming bond pairs
print(s.bonds.N)
s.bonds.types = ['AA_bond']
s.bonds.typeid = [0]*(nbonds-(n1+n2-1))							 
print(len(s.bonds.typeid))
bond_pairs=np.zeros((nbonds,2),dtype=int)

for i in range(nbonds):
    count+=1
    bond_pairs[i,:]=np.array([i,i+1])
print(bond_pairs)

for i in range(x1*n1):										   #Setting unwanted bond pairs as 0 for chain1
    if bond_pairs[i][1] % x1 == 0:                             
        bond_pairs[i][0]=0
        bond_pairs[i][1]=0

for i in range((x1*n1),(num_total-1)):						   #Setting unwanted bond pairs as 0 for chain2 
    if bond_pairs[i][1] % x2 == 0:                             
        bond_pairs[i][0]=0
        bond_pairs[i][1]=0


#Final bond pair list	
new_bond_pairs=[]
for i in range(len(bond_pairs)):
        if bond_pairs[i][1]!=0:
            new_bond_pairs.append(bond_pairs[i])
for i in range(len(new_bond_pairs)):
    new_bond_pairs[i]= list(new_bond_pairs[i])

new_bond_pairs = np.array(new_bond_pairs)
s.bonds.group = new_bond_pairs
print(s.bonds.group)
print(len(s.bonds.group))
print(s.particles.typeid)
print(len(s.particles.typeid))

# Box size
s.configuration.dimensions=3
s.configuration.box=[box_length,box_length,box_length,0,0,0]
s.configuration.step=0


# Write gsd file (before resizing)
f = gsd.hoomd.open(name='start.gsd', mode='wb')
f.append(s)
f.close()

# Replicate single to enough number of chains for a slab 
hoomd.context.initialize("--mode=gpu")
system = hoomd.init.read_gsd('start.gsd')


snap = system.take_snapshot(all=True)
print (snap.particles.N)


#Bonds
harmonic=hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('AA_bond',k=8368,r0=bond_length)

# Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', 'body'])


nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)
nb.pair_coeff.set('A','A', lam=1.00, epsilon=0.8368, sigma=0.5, r_cut=2.0)
nb.pair_coeff.set('A','B', lam=0.00, epsilon=0.8368, sigma=0.5, r_cut=2.0)
nb.pair_coeff.set('B','B', lam=0.00, epsilon=0.8368, sigma=0.5, r_cut=2.0)
      
# Electrostatics
# 8.98755e9*1e9*1.6e-19**2*6.02e23/1000./80.=1.73136
yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
yukawa.pair_coeff.set(['A','B'],['A','B'], epsilon=0.0,kappa=1.0)

## Group Particles
all = hoomd.group.all()

## Set up integrator
hoomd.md.integrate.mode_standard(dt=0.01) # Time units in ps

kTinput=300 * 8.3144598/1000.
integrator = hoomd.md.integrate.langevin(group=all, kT=kTinput, seed=63535)
hoomd.run(tsteps=1e6)

##########
#Change the dimensions you need for final resized box

hoomd.update.box_resize(Lx=hoomd.variant.linear_interp([(5e4,100),(1e6,10)]),Ly=hoomd.variant.linear_interp([(5e4,100),(1e6,10)]),Lz=hoomd.variant.linear_interp([(5e4,200),(1e6,10)]),scale_particles=True)
#########

gamma=0.01
for i in aakeys:
    integrator.set_gamma(seq, gamma=0.1)

## Outputs
hoomd.analyze.log(filename='seq.log', quantities=['bond_harmonic_energy','pair_yukawa_energy','potential_energy','kinetic_energy','temperature','lx','ly','lz'], period=10, overwrite=True, header_prefix='#')
restart=hoomd.dump.gsd('resized_config.gsd', period=1000, dynamic=['protperty','momentum'], group=all, overwrite=True)


## Run simulation
hoomd.run(tsteps=1e6)

print('This is the sequence coded in numbers:',chain_id)
count=0
for i in chain_id:
    if i ==1:
        count+=1
print(count)


