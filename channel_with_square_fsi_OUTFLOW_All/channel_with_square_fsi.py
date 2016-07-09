from __future__ import print_function, absolute_import, division 
import sys
import timeit

sys.path.append("..")
print(sys.path)

from numpy import *
from pyKratos import *

#this i a modified copy of the stokes testcase to test the first prototype of the navier stokes element (steady) 

#messure runtime
start = timeit.default_timer()

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    VELOCITY_X,
    VELOCITY_Y,
    ACCELERATION_X,
    ACCELERATION_Y,
    PRESSURE,
    IS_LAGRANGIAN,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y,
    IS_STRUCTURE
]

property_list = {
    0: {VISCOSITY: 1.5e-5,
        DENSITY: 1.0,
        BODY_FORCE_X: 0.0,
        BODY_FORCE_Y: 0.0,
        }
}

#create buffer size and model part
buffer_size = 2  # store current step and 1 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
model_part.AddProperties(property_list)

#read in geometry from input - nodes and elements
input_file = "02_channel_w_squareFluid.mdpa"
#input_file = "triangle_testFluid.mdpa"
output_file = "02_channel_w_squareFluid.OUT"
dragLiftFile= "dragLift.txt"

#initialize GiDIO object
from pyKratos import gid_io_navier_stokes
gid_io_input = gid_io_navier_stokes.GidIONS(input_file, output_file)
#read model from mdpa file
gid_io_input.ReadModelPart(model_part)


import bossak_scheme
alphaBossak = -0.3;
time_scheme = bossak_scheme.BossakScheme(model_part, alphaBossak)


builder_and_solver = builder_and_solver.BuilderAndSolver(
    model_part, time_scheme)

from pyKratos import newton_raphson_strategy
strategy = newton_raphson_strategy.NewtonRaphsonStrategy(
    model_part, time_scheme, builder_and_solver, 100, 1.0e-5, 1.0e-3)

strategy.Initialize()
#open file for drag and lift computation
strategy.OpenFile(dragLiftFile)

mesh_name = "test"
gid_io_input.WriteMesh(model_part, mesh_name)
dt = 0.01
nsteps = 50 #2000
outputStep = 10
step = 1

# WEI BEGIN =========================================
# define the variables for structural part
# vertical
M = 10.0
K = 200.0
pInf = 1.0

u0 = 0.0
v0 = 0.0
a0 = 0.0
# rotational
I = M/12.0    # I = 1/12*bh^3
Kr = 0.0

phi0 = 0.0
omega0 = 0.0
alpha0 = 0.0

# structural solver
from oscillatorGeneralizedAlpha1 import StructureSDOF
# vertical
structure1 = StructureSDOF(dt, M, K, pInf, u0, v0, a0)
# rotational
structure2 = StructureSDOF(dt, I, Kr, pInf, phi0, omega0, alpha0)
u = []
phi = []
t = []
# WEI END ===========================================

#to be consistend with kratos results
model_part.CloneTimeStep(dt)
model_part.CloneTimeStep(2*dt)

for i in range(3,nsteps):
    time = i*dt
    model_part.CloneTimeStep(time)
    print("time = ", time)
    
    strategy.Solve()
    
    #check if this step results are written
    if(step >= outputStep):
        gid_io_input.WriteNodalResults(PRESSURE, model_part.NodeIterators(), time)
        gid_io_input.WriteNodalResults(VELOCITY, model_part.NodeIterators(), time)
        #gid_io_input.WriteNodalResults(ACCELERATION, model_part.NodeIterators(), time)
        step = 0
    else:
        step = step + 1
        
    # WEI BEGIN =========================================
    # fluid force
    dragLift = [0.0, 0.0]
    moment = builder_and_solver.ComputeReactions(strategy.A, strategy.dx, strategy.b, dragLift)
    print("Fy = ",-dragLift[1])
    print("Moment = ", -moment)
    # solve structure problem --- vertical movement
    structure1.updateStructureTimeStep()
    structure1.solveStructure(-dragLift[1])
    structure1.printValuesAtCurrentStep(step)
    # solve structure problem --- vertical movement
    structure2.updateStructureTimeStep()
    structure2.solveStructure(-moment)
    structure2.printValuesAtCurrentStep(step)
    # record the displacement
    u.append(structure1.getDisplacement())
    phi.append(structure2.getDisplacement())
    t.append(time)
    
    # apply Dirichlet B.C. for fluid by setting solution value at boundary nodes
    for fluid_node in builder_and_solver.model_part.NodeIterators():
            if(fluid_node.GetSolutionStepValue(IS_STRUCTURE,0)):
                # if rotation is considered, we need to compute the velocity at each node
                vel = structure1.v1 + structure2.v1 * fluid_node.coordinates[0]
                fluid_node.SetSolutionStepValue(VELOCITY_Y, 0, vel) 
    # WEI END ===========================================

strategy.CloseFile()

#stop timer
stop = timeit.default_timer()

LogFile = open("LOG.txt","w")
resultTime = "Program runntime: " + str(stop - start)
LogFile.write(resultTime)
LogFile.close()
print("Program runntime: ",stop - start)

# WEI BEGIN ===========================================
# plot the result for structure
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.subplot(111)
ax.set_ylim(min(u+phi), max(u+phi))
ax.set_xlim(0.0, nsteps*dt)
ax.grid()

ax.plot(t, u, label='Vertical displacement')
ax.plot(t, phi, label='Rotational displacement')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()
# WEI END =============================================

# End of channel_with_square.py
