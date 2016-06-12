from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters
import define_output
from numpy import *
import h5py


#
#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

#
#
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ALEApplication import *

# defining variables to be used

variables_dictionary = {"DISPLACEMENT": DISPLACEMENT,
                        "PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "REACTION": REACTION,
                        "DISTANCE": DISTANCE,
                        "Q_VALUE": Q_VALUE,
                        "VORTICITY_MAGNITUDE": VORTICITY_MAGNITUDE}

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

import trilinos_mesh_solver_structural_similarity as fluid_mesh_solver

#
#
# importing variables
solver_module.AddVariables(fluid_model_part, SolverSettings)

# mesh_velocity, displacement
fluid_mesh_solver.AddVariables(fluid_model_part)

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)

# do parallel reading ######################
number_of_partitions = mpi.size  # we set it equal to the number of processors
if mpi.rank == 0:
    partitioner = MetisDivideHeterogeneousInputProcess(
        model_part_io_fluid,
        number_of_partitions,
        domain_size,
        1,
        True)
    partitioner.Execute()

mpi.world.barrier()

MPICommSetup = SetMPICommunicatorProcess(fluid_model_part)
MPICommSetup.Execute()

my_input_filename = input_file_name + "_" + str(mpi.rank)
model_part_io_fluid = ModelPartIO(my_input_filename)
model_part_io_fluid.ReadModelPart(fluid_model_part)
#

Comm = CreateCommunicator()


# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)


# adding dofs
solver_module.AddDofs(fluid_model_part, SolverSettings)
fluid_mesh_solver.AddDofs(fluid_model_part)

#
#
# ALE boundary conditions
for node in fluid_model_part.GetNodes(2): # outlet
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)
for node in fluid_model_part.GetNodes(5): # ground
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)
for node in fluid_model_part.GetNodes(4): # mirror_back
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)
for node in fluid_model_part.GetNodes(3): # mirror_front
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)
for node in fluid_model_part.GetNodes(6): # side
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)
for node in fluid_model_part.GetNodes(7): # top
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)
for node in fluid_model_part.GetNodes(8): # inlet
    node.Fix(DISPLACEMENT_X)
    node.Fix(DISPLACEMENT_Y)
    node.Fix(DISPLACEMENT_Z)

#
#
# setting FSI interface
mirrorBackIds = set()
for node in fluid_model_part.GetNodes(4): # mirror_back
    node.Set(INTERFACE, True)
    mirrorBackIds.add(node.Id)
mirrorFrontIds = set()
for node in fluid_model_part.GetNodes(3): # mirror_front
    node.Set(INTERFACE, True)
    mirrorFrontIds.add(node.Id)
for cond in fluid_model_part.Conditions:
    if cond.GetNode(0).Is(INTERFACE) and cond.GetNode(1).Is(INTERFACE) and cond.GetNode(2).Is(INTERFACE):
        cond.Set(INTERFACE, True)
        if cond.GetNode(0).Id in mirrorFrontIds and cond.GetNode(1).Id in mirrorFrontIds \
                and cond.GetNode(2).Id in mirrorFrontIds:
            fluid_model_part.AddCondition(cond, 3)
        elif cond.GetNode(0).Id in mirrorBackIds and cond.GetNode(1).Id in mirrorBackIds \
                and cond.GetNode(2).Id in mirrorBackIds:
            fluid_model_part.AddCondition(cond, 4)

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
if(SolverSettings.TurbulenceModel == "Spalart-Allmaras"):
    # apply the initial turbulent viscosity on all of the nodes
    turb_visc = SolverSettings.TurbulentViscosity
    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
        visc = node.GetSolutionStepValue(VISCOSITY)
        node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
        if node.IsFixed(VELOCITY_X):
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
# fix outlet BCs to prevent instabilities
for node in fluid_model_part.GetNodes(2): # outlet
    if node.X > 163.83 and node.Z < 0.01:
        node.SetValue(IS_STRUCTURE,0.0)
        node.Fix(VELOCITY_Z)
        node.Fix(VELOCITY_Y)

# creating the mesh motion solver
mesh_solver = fluid_mesh_solver.TrilinosMeshSolverStructuralSimilarity(fluid_model_part, domain_size, False)
mesh_solver.time_order = ProjectParameters.FluidSolverConfiguration.time_order
mesh_solver.Initialize()

# initialize GiD  I/O
from trilinos_gid_output import TrilinosGiDOutput
gid_io = TrilinosGiDOutput(input_file_name,
                           ProjectParameters.VolumeOutput,
                           ProjectParameters.GiDPostMode,
                           ProjectParameters.GiDMultiFileFlag,
                           ProjectParameters.GiDWriteMeshFlag,
                           ProjectParameters.GiDWriteConditionsFlag)

if not ProjectParameters.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(fluid_model_part, cut_list)

gid_io.initialize_results(fluid_model_part)

#
#
# define the drag computation list
drag_list = define_output.DefineDragList()
drag_file_output_list = []

if(mpi.rank == 0):
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
        nodes = fluid_model_part.GetNodes(it[0])
        rx = 0.0
        ry = 0.0
        rz = 0.0
        mx = 0.0
        my = 0.0
        mz = 0.0
        for node in nodes:
            if(node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
                reaction = node.GetSolutionStepValue(REACTION, 0)
                rx -= reaction[0]
                ry -= reaction[1]
                rz -= reaction[2]
                x = node.X - 54.891
                y = node.Y - 65.190
                z = node.Z -  2.840
                mx -= y * reaction[2] - z * reaction[1]
                my -= z * reaction[0] - x * reaction[2]
                mz -= x * reaction[1] - y * reaction[0]
        auxrx = mpi.gather(mpi.world,rx,0)
        auxry = mpi.gather(mpi.world,ry,0)
        auxrz = mpi.gather(mpi.world,rz,0)
        auxmx = mpi.gather(mpi.world,mx,0)
        auxmy = mpi.gather(mpi.world,my,0)
        auxmz = mpi.gather(mpi.world,mz,0)
        rx = 0.0
        ry = 0.0
        rz = 0.0
        mx = 0.0
        my = 0.0
        mz = 0.0
        for k in auxrx:
            rx += k
        for k in auxry:
            ry += k
        for k in auxrz:
            rz += k
        for k in auxmx:
            mx += k
        for k in auxmy:
            my += k
        for k in auxmz:
            mz += k
        if(mpi.rank == 0):
            output = str(time) + " " + str(rx) + " " + str(ry) + " " + str(rz) + " " + str(mx) + " " + str(my) + " " + str(mz) + "\n"
            drag_file_output_list[i].write(output)
            drag_file_output_list[i].flush()
        i = i + 1

#
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

time = 2.0 * Dt
out = 0
step = 2

normal_util = NormalCalculationUtils()
normal_util.CalculateOnSimplex(fluid_model_part,domain_size)

class StructureSDOF:

    """ Direct time integration of linear SDOF (Generalized-alpha method). """

    def __init__(self, Dt, J, zeta, K, pInf, u0, v0, a0):
        self.Dt = Dt
        # structure moment of inertia, damping and spring stiffness
        self.J = J
        self.B = 2.0 * sqrt(J * K) * zeta # zeta is the damping ratio
        self.K = K
        # generalized alpha parameters (to ensure unconditional stability, 2nd order accuracy)
        self.alphaM = (2.0 * pInf - 1.0) / (pInf + 1.0)
        self.alphaF = pInf / (pInf + 1.0)
        self.beta = 0.25 * (1 - self.alphaM + self.alphaF)**2
        self.gamma = 0.5 - self.alphaM + self.alphaF
        # coefficients for LHS
        self.a1h = (1.0 - self.alphaM) / (self.beta * self.Dt**2)
        self.a2h = (1.0 - self.alphaF) * self.gamma / (self.beta * self.Dt)
        self.a3h = 1.0 - self.alphaF
        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.Dt
        self.a3m = (1.0 - self.alphaM - 2.0 * self.beta) / (2.0 * self.beta)
        # coefficients for damping
        self.a1b = (1.0 - self.alphaF) * self.gamma / (self.beta * self.Dt)
        self.a2b = (1.0 - self.alphaF) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alphaF) * (0.5 * self.gamma / self.beta - 1.0) * self.Dt
        # coefficient for stiffness
        self.a1k = -1.0 * self.alphaF
        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.Dt)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.Dt
        # coefficients for acceleration update
        self.a1a = self.a1v / (self.Dt * self.gamma)
        self.a2a = -1.0 / (self.beta * self.Dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)
        # initial displacement, velocity and acceleration
        self.u0 = u0
        self.v0 = v0
        self.a0 = a0
        #
        self.u1 = u0
        self.v1 = v0
        self.a1 = a0
        # moment from a previous time step (initial moment)
        self.m0 = J * a0 + self.B * v0 + K * u0
        self.m1 = J * a0 + self.B * v0 + K * u0
        # output
        if mpi.rank == 0:
            self.support_output = open("mirror_structure.dat", 'w')
            self.support_output.write("# (1): time [s] (2): rotation [degrees] (3): support moment [N m]\n")
            self.support_output.flush()

    def getRotation(self):
        return self.u1

    def predictRotation(self):
        return 2.0 * self.u1 - self.u0

    def printSupportOutput(self, time):
        if mpi.rank == 0:
            self.support_output.write(str(time) + " " + str(self.u1 * 180.0 / pi) + " " + str(self.K * self.u1) + "\n")
            self.support_output.flush()
    
    def solveStructure(self, m1):
        # sys of eq reads: LHS * u1 = RHS
        # m1 is the externally applied moment at the current time step
        F = (1.0 - self.alphaF) * m1 + self.alphaF * self.m0
        LHS = self.a1h * self.J + self.a2h * self.B + self.a3h * self.K
        RHS = self.J * (self.a1m * self.u0 + self.a2m * self.v0 + self.a3m * self.a0)
        RHS += self.B * (self.a1b * self.u0 + self.a2b * self.v0 + self.a3b * self.a0)
        RHS += self.a1k * self.K * self.u0 + F
        #
        # update self.m1
        self.m1 = m1
        # updates self.u1,v1,a1
        self.u1 = RHS / LHS
        self.v1 = self.a1v * (self.u1 - self.u0) + self.a2v * self.v0 + self.a3v * self.a0
        self.a1 = self.a1a * (self.u1 - self.u0) + self.a2a * self.v0 + self.a3a * self.a0

    def incrementTimeStep(self):
        # update angular displacement, velocity and acceleration
        self.u0 = self.u1
        self.v0 = self.v1
        self.a0 = self.a1
        # update the moment
        self.m0 = self.m1

def RotateMirror(theta, nodes):
    for node in nodes:
        rx0 = node.X0 - 54.89128
        rz0 = node.Z0 - 2.84
        rx  = cos(theta) * rx0 + sin(theta) * rz0
        rz  =-sin(theta) * rx0 + cos(theta) * rz0
        dx  = rx - rx0
        dz  = rz - rz0
        node.SetSolutionStepValue(DISPLACEMENT_X,dx)
        node.SetSolutionStepValue(DISPLACEMENT_Z,dz)

def ExtractMoment(nodes):
    my = 0.0
    for node in nodes:
        if(node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
            reaction = node.GetSolutionStepValue(REACTION, 0)
            fx = -reaction[0]
            fz = -reaction[2]
            rx = node.X - 54.891
            rz = node.Z -  2.840
            my += rz * fx - rx * fz
    auxmy = mpi.allgather(mpi.world,my)
    my = 0.0
    for k in auxmy:
        my += k
    return my

#structure_solver = StructureSDOF(Dt, 3321.0, 0.05, 16060.7, 0.16, 0.0, 0.0, 0.0) # f = 0.35
structure_solver = StructureSDOF(Dt, 3321.0, 0.05, 27742.4, 0.16, 0.0, 0.0, 0.0) # f = 0.46
#structure_solver = StructureSDOF(Dt, 3321.0, 0.05, 41115.4, 0.16, 0.0, 0.0, 0.0) # f = 0.56
#structure_solver = StructureSDOF(Dt, 3321.0, 0.05, 57110.5, 0.16, 0.0, 0.0, 0.0) # f = 0.66
#structure_solver = StructureSDOF(Dt, 3321.0, 0.05, 77733.8, 0.16, 0.0, 0.0, 0.0) # f = 0.77

class H5UnsteadyBC:

    """ Add unsteady turbulent fluctuations in HDF5 format to BC. """

    def __init__(self, inlet_nodes, y0, z0):
        self.H5_filename = "../../../ptsc10mps.h5"
        self.H5_file     = h5py.File(self.H5_filename, 'r')
        self.lx          = self.H5_file.get('lx')[0]
        self.ly          = self.H5_file.get('ly')[0]
        self.lz          = self.H5_file.get('lz')[0]
        self.u           = self.H5_file.get('u')
        self.v           = self.H5_file.get('v')
        self.w           = self.H5_file.get('w')
        z                = self.H5_file.get('z')[0]
        self.log_z0      = self.H5_file.get('log_z0')[0]
        umean            = self.H5_file.get('umean')[0]
        self.ubulk       = umean * (log(73.8/self.log_z0) - 1.0) / log(z/self.log_z0) # average over CFD domain height         
        self.utau        = umean * 0.41 / log(z / self.log_z0)
        self.nx          = self.u.shape[0]
        self.ny          = self.u.shape[1]
        self.nz          = self.u.shape[2]
        self.dx          = self.lx / (self.nx - 1)
        self.dy          = self.ly / (self.ny - 1)
        self.dz          = self.lz / (self.nz - 1)
        self.y0          = y0
        self.z0          = z0
        self.inlet_nodes = inlet_nodes
        if len(self.inlet_nodes) > 0:
            # y-direction
            self.ymin = min([node.Y for node in self.inlet_nodes])
            self.ymax = max([node.Y for node in self.inlet_nodes])
            if self.ymax - self.ymin > self.ly:
                raise RuntimeError('H5InletBC: Found node.Y out of range')
            self.ja = int(floor((self.ymin - self.y0) / self.dy))
            self.je = int( ceil((self.ymax - self.y0) / self.dy)) + 1
            self.jj = self.je - self.ja
            # z-direction
            self.zmin = min([node.Z for node in self.inlet_nodes])
            self.zmax = max([node.Z for node in self.inlet_nodes])
            if self.zmax - self.zmin > self.lz:
                raise RuntimeError('H5InletBC: Found node.Z out of range')
            self.ka = int(floor((self.zmin - self.z0) / self.dz))
            self.ke = int( ceil((self.zmax - self.z0) / self.dz)) + 1
            self.kk = self.ke - self.ka
            # x-slice
            self.local_inlet_component = zeros((self.jj, self.kk))
            self.x  = 0.0 # start
    def ScaleBC(self, scal):
        for node in self.inlet_nodes:
            vel = node.GetSolutionStepValue(VELOCITY_X)
            node.SetSolutionStepValue(VELOCITY_X, scal * vel)
            vel = node.GetSolutionStepValue(VELOCITY_Y)
            node.SetSolutionStepValue(VELOCITY_Y, scal * vel)
            vel = node.GetSolutionStepValue(VELOCITY_Z)
            node.SetSolutionStepValue(VELOCITY_Z, scal * vel)
    def Update(self, dt):
        if len(self.inlet_nodes) > 0:
            for node in self.inlet_nodes:
                vel = 2.439 * self.utau * log((node.Z + self.log_z0) / self.log_z0)
                node.SetSolutionStepValue(VELOCITY_X, vel)
                node.SetSolutionStepValue(VELOCITY_Y, 0.0)
                node.SetSolutionStepValue(VELOCITY_Z, 0.0)
            self.x = self.x + dt * self.ubulk
            if self.x >= self.lx:
                self.x = self.x - self.lx
            i = int(floor(self.x / self.dx))
            tx = (self.x - self.dx * i) / self.dx
            self.local_inlet_component = (1.0 - tx) *            \
                self.u[i,self.ja:self.je,self.ka:self.ke] +      \
                tx * self.u[i+1,self.ja:self.je,self.ka:self.ke]
            self.AddComponent(VELOCITY_X)
            self.local_inlet_component = (1.0 - tx) *            \
                self.v[i,self.ja:self.je,self.ka:self.ke] +      \
                tx * self.v[i+1,self.ja:self.je,self.ka:self.ke]
            self.AddComponent(VELOCITY_Y)
            self.local_inlet_component = (1.0 - tx) *            \
                self.w[i,self.ja:self.je,self.ka:self.ke] +      \
                tx * self.w[i+1,self.ja:self.je,self.ka:self.ke]
            self.AddComponent(VELOCITY_Z)
    def AddComponent(self, component):
        for node in self.inlet_nodes:
            y = node.Y - self.y0
            z = node.Z - self.z0
            ty = (y % self.dy) / self.dy # [0,1]
            tz = (z % self.dz) / self.dz
            ty = 2.0 * (ty - 0.5)        # [-1,1]
            tz = 2.0 * (tz - 0.5)
            wt = [0.25*(1.0-ty)*(1.0-tz), 0.25*(1.0+ty)*(1.0-tz), \
                  0.25*(1.0+ty)*(1.0+tz), 0.25*(1.0-ty)*(1.0+tz)]
            j = int(floor(y / self.dy)) - self.ja
            k = int(floor(z / self.dz)) - self.ka
            vel = node.GetSolutionStepValue(component)
            vel = wt[0] * self.local_inlet_component[j,k]     \
                + wt[1] * self.local_inlet_component[j+1,k]   \
                + wt[2] * self.local_inlet_component[j+1,k+1] \
                + wt[3] * self.local_inlet_component[j,k+1]   \
                + vel
            node.SetSolutionStepValue(component, vel)
            
turbul_inlet = H5UnsteadyBC(fluid_model_part.GetNodes(8), 0.0, 0.0)            

#
#
# IAC parameter
fluid_model_part.ProcessInfo[DENSITY] = 8000.0 * 0.004

#
#
# initialize simulation
fluid_model_part.CloneTimeStep(0.0 * Dt)
fluid_model_part.CloneTimeStep(1.0 * Dt)
fluid_model_part.CloneTimeStep(2.0 * Dt)
sys.stdout.flush()
mpi.world.barrier()

while(time <= final_time):

    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    turbul_inlet.Update(Dt)
    if time < 30.0:
        turbul_inlet.ScaleBC(float(time) / 30.0)

    if mpi.rank == 0:
        print("STEP = ", step)
        print("TIME = ", time)
        sys.stdout.flush()

    k = 0
    # predict structure rotation
    pitching_angle = structure_solver.predictRotation()
    # solve the fluid problem
    RotateMirror(pitching_angle, fluid_model_part.GetNodes(3))
    RotateMirror(pitching_angle, fluid_model_part.GetNodes(4))
    mesh_solver.Solve()
    normal_util.CalculateOnSimplex(fluid_model_part,domain_size)
    fluid_solver.Solve()
    # solve the structure problem
    pitching_moment = 0.0
    pitching_moment += ExtractMoment(fluid_model_part.GetNodes(3))
    pitching_moment += ExtractMoment(fluid_model_part.GetNodes(4))
    structure_solver.solveStructure(pitching_moment)
    # compute the residual
    initial_residual = structure_solver.getRotation() - pitching_angle
    old_residual = initial_residual
    residual = old_residual
    
    while not (abs(residual) < 1.0e-09 or abs(residual) < abs(initial_residual) * 1.0e-02):
        relax_coef = 1.0
        # update structure rotation
        pitching_angle += relax_coef * residual
        k = k + 1
        # solve the fluid problem
        RotateMirror(pitching_angle, fluid_model_part.GetNodes(3))
        RotateMirror(pitching_angle, fluid_model_part.GetNodes(4))
        mesh_solver.Solve()
        normal_util.CalculateOnSimplex(fluid_model_part,domain_size)
        fluid_solver.Solve()
        # solve the structure problem
        pitching_moment = 0.0
        pitching_moment += ExtractMoment(fluid_model_part.GetNodes(3))
        pitching_moment += ExtractMoment(fluid_model_part.GetNodes(4))
        structure_solver.solveStructure(pitching_moment)
        # compute the residual
        old_residual = residual
        residual = structure_solver.getRotation() - pitching_angle

        if mpi.rank == 0:
            print('Iteration[', k, ']: relax_coef = ', relax_coef, ' rel. residual = ', 
                  abs(residual / initial_residual), ' abs. residual = ', abs(residual))
            sys.stdout.flush()
        if k >= 15:
            break
    
    structure_solver.printSupportOutput(time)

    graph_printer.PrintGraphs(time)
    PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if(output_time <= out):
        gid_io.write_results(
            time,
            fluid_model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)
        out = 0
    elif (time >= 100.0 and time <= 110.0):
        gid_io.write_results(
            time,
            fluid_model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)
        out = 0

    out = out + Dt
    structure_solver.incrementTimeStep()

gid_io.finalize_results()

for i in drag_file_output_list:
    i.close()

if mpi.rank == 0:
    structure_solver.support_output.close()
