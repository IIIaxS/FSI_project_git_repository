from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
print(sys.path)

from numpy import *
from pyKratos import *

#this i a modified copy of the stokes testcase to test the first prototype of the navier stokes element (steady) 

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
    0: {VISCOSITY: 1.0e-6,
        DENSITY: 1000.0,
        BODY_FORCE_X: 0.0,
        BODY_FORCE_Y: 0.0,
        }
}

#create buffer size and model part
buffer_size = 2  # store current step and 1 in the past
model_part = ModelPart(buffer_size, solution_step_variables)
model_part.AddProperties(property_list)

#read in geometry from input - nodes and elements
input_file = "NS_Bossak_channel.mdpa"
#input_file = "triangle_testFluid.mdpa"
output_file = "NS_Bossak_channel_MDPA.OUT"
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
dt = 0.1
nsteps = 20
for i in range(1,nsteps):
    time = i*dt
    model_part.CloneTimeStep(time)
    print("time = ", time)
    strategy.Solve()
    gid_io_input.WriteNodalResults(PRESSURE, model_part.NodeIterators(), time)
    gid_io_input.WriteNodalResults(VELOCITY, model_part.NodeIterators(), time)
    gid_io_input.WriteNodalResults(ACCELERATION, model_part.NodeIterators(), time)


strategy.CloseFile()