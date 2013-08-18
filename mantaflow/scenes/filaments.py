#
# Filament test
# don't use yet --- work in progress

from manta import *
from math import *
import random
import math

# solver params
res = 64
gs = vec3(res,2*res,res)
s = Solver(name='main', gridSize = gs,is_debug=1)
s.timestep = 0.05

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)
filaments = s.create(VortexFilamentSystem)
tracer = s.create(TracerParticleSystem)

# scene setup
flags.initDomain(boundaryWidth=1)
flags.fillGrid()

#filaments.addRing(position=gs*vec3(0.5,0.1,0.5), circulation=10, radius=0.2*res, normal=(0,1,0), number=40)
#filaments.addRing(position=gs*vec3(0.5,0.14,0.5), circulation=10, radius=0.2*res, normal=(0,1,0), number=1000)
#filaments.remesh(maxLen=3, minLen=1)
    
#filaments.setPos(10,filaments.getPos(10)+vec3(0,1.5,0))
source = s.create(Cylinder, center=gs*vec3(0.5,0.09,0.5), radius=res*0.17, z=gs*vec3(0, 0.03, 0))
fixedRegion = s.create(Box, center=gs*vec3(0.5,0.07,0.5), size=gs*vec3(0.4,0.03,0.4))
mesh.fromShape(source)

if (GUI):
    gui = Gui()
    gui.show()
    gui.pause()

r = 0.17
scale = 1
number_seg=41
min_num_edge=number_seg/5
max_num_edge=4*number_seg
cosine_threshold=0.9
ring_radius=res*0.13
theta = random.gauss(0,0.01*pi)
phi = random.gauss(0,0.01*pi)
n = vec3(sin(theta)*cos(phi), cos(theta), sin(theta)*sin(phi))
#for k in range(20):
 #   filaments.addRing(position=gs*vec3(0.5,0.12,0.5),circulation=random.uniform(100,100),radius=res*random.gauss(r,0.01),normal=n,number=number_seg)
#main loop
for t in range(1000):
    markAsFixed(mesh=mesh, shape=fixedRegion)

    # seed rings
    #if t< 20 and random.random()<0.05:
    #if random.random()<0.05:
    #theta = random.gauss(0,0.01*pi)
    #phi = random.gauss(0,0.01*pi)
    #n = vec3(sin(theta)*cos(phi), cos(theta), sin(theta)*sin(phi))
    #filaments.addRing(position=gs*vec3(0.5,0.14,0.5),circulation=random.uniform(60,60),radius=res*random.gauss(r,0.01),normal=n,number=3)
    
    if t %2==0:
        theta = random.gauss(0,0.01*pi)
        phi = random.gauss(0,0.01*pi)
        n = vec3(sin(theta)*cos(phi), cos(theta), sin(theta)*sin(phi))
        filaments.addRing(position=gs*vec3(0.5,0.09,0.5),circulation=random.uniform(14,45),radius=ring_radius,normal=n,number=number_seg)


    # seed tracer particles
    for i in range(80):
        x=100
        z=100
        #while x*x+z*z > r*r :
        x=random.uniform(-r/5,r/5)
        y=random.uniform(0,0.1)
        z=random.uniform(-r/5,r/5)
        tracer.addParticle(vec3(x+0.5,y+0.1,z+0.5)*gs)
        #tracer.addParticle(vec3(x+0.5,y+0.1,z+0.5)*gs)
    #if t%1==0:
        #tracer.add_Particle(gs*vec3(0.5,0.09,0.5),radius=res*0.17,normal=n,num=number_seg)
    
    # mesh cosmetics
    
        

    #gui.pause()
    #filaments.ddTest(d=0.25,phi=0.1)
    edge_length=math.sin(math.pi/number_seg)*ring_radius*2;

    #filaments.advectMesh(mesh=mesh, scale=1, regularization=2, integrationMode=IntRK4)




    filaments.advectParticles(sys=tracer, scale=1, regularization=2, integrationMode=IntRK4)
    if t>0:
        filaments.advectSelf(scale=1, regularization=2, integrationMode=IntRK4)
    #if t%2==0:
    #    filaments.Direct_motion(scale=1,regularization=2,integrationMode=IntRK4)


    #filaments.doublyDiscreteUpdate(regularization=2)
    #filaments.Decimate_ring(min_num_edge,min_num_edge*edge_length,max_num_edge*edge_length)
    #filaments.split_ring(-cosine_threshold,edge_length/2)
    max_edge_length=3*edge_length;
    #filaments.merge_adj_edge(cosine_threshold,max_edge_length)
    #filaments.divide_ring(edge_length*2,max_edge_length);
    filaments.resample_ring(max_edge_length,edge_length,number_seg);

    #filaments.reconnect_ring(-cosine_threshold,edge_length/3)
    filaments.reset_dirty()
    #filaments.revise_circulation();

    #filaments.remesh(maxLen=2*edge_length, minLen=edge_length/2)
    #smoothMesh(mesh=mesh, strength=1e-4, steps=1)
    #subdivideMesh(mesh=mesh, minAngle=0.01, minLength=scale, maxLength=2*scale, cutTubes=True)
    #killSmallComponents(mesh=mesh, elements=10)

    #filaments.Debug_fun()
    
    s.printTimings()
    s.step()
    
