import sys 
import gsd, gsd.hoomd, gsd.pygsd 

lastframe = gsd.hoomd.HOOMDTrajectory(gsd.pygsd.GSDFile(open(sys.argv[1],'rb')))[-1]

lastframe.configuration.step = 0

#print ('Box size is %f x %f x %f')%(lastframe.configuration.box[0],lastframe.configuration.box[1],lastframe.configuration.box[2])
newgsdfile=gsd.hoomd.open('%sx%sx%s_box.gsd'%(lastframe.configuration.box[0],lastframe.configuration.box[1],lastframe.configuration.box[2]),'wb')
newgsdfile.append(lastframe)
newgsdfile.close()
#gsd.hoomd.create('lastframe%f.gsd'%(lastframe.configuration.box[0]),snapshot=lastframe)