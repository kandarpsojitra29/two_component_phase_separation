import gsd, gsd.hoomd, gsd.pygsd 
import sys 
import numpy as np 
def extend(s):
    boxdim = s.configuration.box[:3]
    zmin,zmax,dz = -boxdim[2]/2., boxdim[2]/2., boxdim[2]
    pos1 =  s.particles.position
    pos = pos1.copy()
    skip=0
    ncomp=int(input("No of components: "))
    for k in range(ncomp):
        chain_length = int(input("Length of comp %i:"%(k+1)))
        nchain = int(input("No of chains of comp %i:"%(k+1)))
        nres = chain_length
        for i in range(nchain):
            mol_coord = pos[i*nres+skip:(i+1)*nres+skip,2]
            for j in range(1,nres):
                dist2 = (mol_coord[j] - mol_coord[j-1])**2
                if dist2 > 8:
                    excess = np.sign(mol_coord[j] - mol_coord[j-1])*dz
                    mol_coord[j] = mol_coord[j] - excess 
                    #mol_coord[j] -= np.sign(mol_coord[j] - mol_coord[j-1])*dz # checking whether a bond is getting cutoff 
                com = np.mean(mol_coord)
                if com < zmin:
                    mol_coord += dz
                elif com > zmax:
                    mol_coord -= dz
            pos[i*nres+skip:(i+1)*nres+skip,2] = mol_coord
        skip += nchain*nres
    return pos
slab_z_length = float(input("Enter slab length: "))
f = gsd.pygsd.GSDFile(open(sys.argv[1],'rb'))
t = gsd.hoomd.HOOMDTrajectory(f)
s1 = t[0]
s = gsd.hoomd.Snapshot()
s.particles.N = s1.particles.N
s.particles.types = s1.particles.types 
s.particles.typeid = s1.particles.typeid 
s.particles.mass = s1.particles.mass
s.particles.charge = s1.particles.charge
s.particles.position = extend(s1)
s.bonds.N = s1.bonds.N
s.bonds.types = s1.bonds.types
s.bonds.typeid = s1.bonds.typeid
s.bonds.group = s1.bonds.group
s.configuration.box = s1.configuration.box
s.configuration.dimensions=3
s.configuration.box = [s1.configuration.box[0],s1.configuration.box[1],slab_z_length,0,0,0] 
s.configuration.step = 0
outfile = gsd.hoomd.open('box2slab_extend.gsd','wb')
outfile.append(s)
outfile.close()