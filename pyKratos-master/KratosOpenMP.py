from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters
import define_output

# BEGIN IVAN ===================================================================
from numpy import *
from oscillatorGeneralizedAlpha1 import StructureSDOF
from ProgressBar import ProgressBar
import shutil
# END IVAN ==================================================


#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

#
#
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

# BEGIN IVAN ===================================================================
from KratosMultiphysics.ALEApplication import *
from KratosStructuralApplication import *
# END IVAN ===================

# defining variables to be used
# GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
# CODE IS NOT USING IT

variables_dictionary = {"PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "REACTION": REACTION,
                        "DISTANCE": DISTANCE, }

# defining a model part for the fluid
fluid_model_part = ModelPart("FluidPart")

if "REACTION" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#
#
# importing the solvers needed
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_module = import_solver(SolverSettings)

#
#
# importing variables
solver_module.AddVariables(fluid_model_part, SolverSettings)

# BEGIN IVAN ===================================================================
import mesh_solver
mesh_solver.AddVariables(fluid_model_part)
# adding unused Kratos variable to store the mesh accelerations (used in SetMeshVelocityBossak)
fluid_model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL)
# END IVAN ============

# introducing input file name
input_file_name = ProjectParameters.problem_name


# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# Check tetrahedral mesh for wrong orientation
throw_errors = False
orientation_check = TetrahedralMeshOrientationCheck(fluid_model_part,throw_errors)
orientation_check.Execute()

solver_module.AddDofs(fluid_model_part, SolverSettings)

# BEGIN IVAN ===================================================================
mesh_solver.AddDofs(fluid_model_part)
# END IVAN ====

# If Lalplacian form = 2, free all pressure Dofs
# laplacian_form = ProjectParameters.laplacian_form
# if(laplacian_form >= 2):
    # for node in fluid_model_part.Nodes:
        # node.Free(PRESSURE)

# copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL, 0)
    node.SetValue(Y_WALL, y)

#
#
# Creating the fluid solver
fluid_solver = solver_module.CreateSolver(
    fluid_model_part, SolverSettings)

# activate turbulence model
if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
    # apply the initial turbulent viscosity on all of the nodes
    turb_visc = SolverSettings.TurbulentViscosity
    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
        visc = node.GetSolutionStepValue(VISCOSITY)
        node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
        if (node.IsFixed(VELOCITY_X) and node.GetSolutionStepValue(VELOCITY_X, 0) != 0.0) or \
           (node.IsFixed(VELOCITY_Y) and node.GetSolutionStepValue(VELOCITY_Y, 0) != 0.0) or \
           (node.IsFixed(VELOCITY_Z) and node.GetSolutionStepValue(VELOCITY_Z, 0) != 0.0):
            node.Fix(TURBULENT_VISCOSITY)

    # select nodes on the wall
    fluid_solver.wall_nodes = []
    for i in SolverSettings.SA_wall_group_ids:
        # get the nodes of the wall for SA.
        nodes = fluid_model_part.GetNodes(i)
        for node in nodes:
            fluid_solver.wall_nodes.append(node)
            node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
            node.Fix(TURBULENT_VISCOSITY)


fluid_solver.Initialize()
print("fluid solver created")
#
#

# initialize GiD  I/O
from gid_output import GiDOutput
gid_io = GiDOutput(input_file_name,
                   ProjectParameters.VolumeOutput,
                   ProjectParameters.GiDPostMode,
                   ProjectParameters.GiDMultiFileFlag,
                   ProjectParameters.GiDWriteMeshFlag,
                   ProjectParameters.GiDWriteConditionsFlag)

if not ProjectParameters.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(fluid_model_part, cut_list)

gid_io.initialize_results(fluid_model_part)

# 33
# 33
# define the drag computation list
drag_list = define_output.DefineDragList()
drag_file_output_list = []
for it in drag_list:
    f = open(it[1], 'w')
    drag_file_output_list.append(f)
    tmp = "#Drag for group " + it[1] + "\n"
    f.write(tmp)
    tmp = "#time RX RY RZ\n"
    f.write(tmp)
    f.flush()

print(drag_file_output_list)


def PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time):
    i = 0
    for it in drag_list:
        print(it[0])
        nodes = fluid_model_part.GetNodes(it[0])
        drag = Vector(3)
        drag[0] = 0.0
        drag[1] = 0.0
        drag[2] = 0.0
        for node in nodes:
            reaction = node.GetSolutionStepValue(REACTION, 0)
            drag[0] += reaction[0]
            drag[1] += reaction[1]
            drag[2] += reaction[2]

        output = str(time) + " " + str(drag[0]) + " " + str(
            drag[1]) + " " + str(drag[2]) + "\n"
        # print drag_file_output_list[i]
        # print output
        drag_file_output_list[i].write(output)
        drag_file_output_list[i].flush()
        i = i + 1


# 33
# preparing output of point graphs
import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(
    output_nodes_list,
    fluid_model_part,
    variables_dictionary,
    domain_size)


# Stepping and time settings
Dt = ProjectParameters.Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

print("Stepping and time settings\n") 
print("Number of steps = ", Nsteps)
print("Time step = ", Dt)
print("Output time step = ", output_time)
print("Final time = ", final_time)

# BEGIN IVAN ===================================================================
# outlet
for node in fluid_model_part.GetMesh(2).Nodes:
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)    
# object
for node in fluid_model_part.GetMesh(4).Nodes:
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Z)
# walls
for node in fluid_model_part.GetMesh(3).Nodes:
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)
# inlet
for node in fluid_model_part.GetMesh(5).Nodes:
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)        


fluid_mesh_solver = mesh_solver.MeshSolver(fluid_model_part, domain_size, True)
fluid_mesh_solver.time_order = 2
fluid_mesh_solver.Initialize()

# move the FSI mesh, prescribe displacement of nodes at the interface
def MoveFSI(mesh, u):
    disp = Vector(3)
    disp[0] = 0.0
    disp[2] = 0.0
    for node in mesh.Nodes:
        disp[1] = u
        node.SetSolutionStepValue(DISPLACEMENT, 0, disp)

# computes the mesh velocity (WBZ-alpha scheme)        
def SetMeshVelocityBossak(fluid_model_part, dt):
    # Bossak parameter (in Kratos = -0.3)
    alphaB = -0.3
    # Newmark parameters
    gamma = 0.5 - alphaB
    beta = 0.25 * (1.0 - alphaB) ** 2
    # coefficients for velocity update
    a1v = gamma / (beta * dt)
    a2v = 1.0 - gamma / beta
    a3v = (1.0 - gamma / (2 * beta)) * dt
    
    for node in fluid_model_part.Nodes:
        u0 = node.GetSolutionStepValue(DISPLACEMENT, 1)
        u1 = node.GetSolutionStepValue(DISPLACEMENT, 0)
        
        v0 = node.GetSolutionStepValue(MESH_VELOCITY, 1)
        a0 = node.GetSolutionStepValue(ACCELERATION_NULL, 1)
        
        v1 = a1v * (u1 - u0) + a2v * v0 + a3v * a0
        node.SetSolutionStepValue(MESH_VELOCITY, 0, v1)
        a1 = 1.0 / (gamma * dt) * (v1 - v0) - (1.0 - gamma) / gamma * a0
        node.SetSolutionStepValue(ACCELERATION_NULL, 0, a1)

# lists to store the results
times = []
disp = []

# initialize the structural part
Ms = 0.183
Ks = 652.01
pInf = 1.0

# relaxation parameter for FSI coupling
eta = 1.0

# initial conditions
u0 = 0.0
v0 = 0.0
a0 = 0.0

# instantiate an object: structure
structure = StructureSDOF(Dt, Ms, Ks, pInf, u0, v0, a0)
# structure.printSetUp()

# initialize fluid
MoveFSI(fluid_model_part.GetMesh(4), u0)

print("Relaxation (eta): ", eta)
print()

# Mate Ramp Begin
class InletVelocityRamp:

    """ Mate's custom inlet velocity. """

    def __init__(self, inletNodes, inletVelocity, maxTime):
        self.inletVelocity = inletVelocity
        self.maxTime = maxTime
        self.inletNodes = inletNodes

    def ApplyInletVelocity(self, time):
        if time < self.maxTime:
            U = (time / self.maxTime) * self.inletVelocity
        else:
            U = self.inletVelocity
        for node in self.inletNodes:
            node.SetSolutionStepValue(VELOCITY_X, U)
#Mesh(5) is inlet
inletVelocity = InletVelocityRamp(fluid_model_part.GetNodes(5), 5.3, Dt * 100)
v = Vector(3)
v[0] = 0.0
v[1] = 0.0
v[2] = 0.0
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY, 0, v)
    node.SetSolutionStepValue(VELOCITY, 1, v)
    node.SetSolutionStepValue(VELOCITY, 2, v)
print("RampUp initialized")
# Mate Ramp End

# clean up the currentResults.dat
open("currentResults.dat", "w").close()

#input("Press Enter to Start Simulation ...")

# initialize and start Progress bar window
# progressBar = ProgressBar(0, Nsteps)
# progressBar.startProgressBar()

while(time <= final_time):

    time = time + Dt
    step = step + 1
    # update fluid and structure in time
    fluid_model_part.CloneTimeStep(time)
    structure.updateStructureTimeStep()

    print()
    print("STEP = ", step)
    print("TIME = ", time)
    
    # update progress bar
#     progressBar.updateProgressBar(step)
    
    # Mate modification begin
    inletVelocity.ApplyInletVelocity(time)
    # Mate modification end

    if(step >= 3):
              
        oldDisp = structure.getDisplacement()
        structure.printValuesAtCurrentStep(step)
        
        # !!! for smaller tolerances the resulting displacement 
        # amplitudes amplify in time      
        tol = 1e-5     
        iMax = 10      
        # inner loop for iterations of FSI coupling
        for i in range(0, iMax):
            print("\nFSI iteration ", i)      
                   
            MoveFSI(fluid_model_part.GetMesh(4), oldDisp)

            SetMeshVelocityBossak(fluid_model_part, Dt)

            fluid_mesh_solver.Solve()
            fluid_solver.Solve()
            
            # summation over the reaction forces at the interface
            reaction = 0.0
            for node in fluid_model_part.GetMesh(4).Nodes:
                reaction += node.GetSolutionStepValue(REACTION_Y)
                print("reaction = ", reaction)
                         
            print("Superposed reaction = ", reaction)
             
            # solve the structure for the input force from fluid reaction
            structure.solveStructure(-1.0 * reaction)

            # without relaxation
#             newDisp = structure.getDisplacement()

            # apply relaxation
            newDisp = eta * structure.getDisplacement() + (1.0 - eta) * oldDisp
            
            # !!! for smaller tolerances the resulting displacement amplitudes don't match
            if (structure.checkConvergence(oldDisp, newDisp, tol) == 1):
                print("FSI strong coupling converged at iteration = ", i)
                break
            
            if (i == iMax - 1):
                print("Max number of FSI iterations was reached")
        
            oldDisp = newDisp
        # end of inner loop

        # append the values to the result lists
        times.append(time)
        disp.append(structure.getDisplacement())
        
        # append the results to the currentresult.dat
        file1 = open("currentResults.dat", "a")
        file1.write(str(time) + "\t")
        file1.write(str(structure.getDisplacement()) + "\t")
        file1.write("\n")
        file1.close()
        
        # copy to dropbox
        shutil.copy("currentResults.dat", "/home/mate/Documents/Dropbox/MScThesisOwn/" )
        shutil.copy("ResonanceObject1.dat", "/home/mate/Documents/Dropbox/MScThesisOwn/" )
     
# END IVAN ======================================================================
   
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if(output_time <= out):
        gid_io.write_results(
            time,
            fluid_model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)
        out = 0

    out = out + Dt

gid_io.finalize_results()

for i in drag_file_output_list:
    i.close()
    
# BEGIN IVAN ===================================================================
# write the results to the result.dat
file = open("results.dat", "w")
file.write("#Strong FSI coupling")
file.write("\n#time\t")
file.write("#fsi displacement\t")

for l in range(0, len(times)):
    file.write("\n")
    file.write(str(times[l]))
    file.write("\t" + str(disp[l]))
file.close()
# END IVAN ======================================================================

# exit progress bar
# progressBar.exitProgressBar()
