from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = line2d.Line2D(list_of_nodes)
    return RobustOutflowCondition(Id, prop, geom)


class RobustOutflowCondition:
    #this elements constructs a stiffness matrix which mixes velocities and pressures
    #for each pair of nodes I and J this matrix has a 3*3 subblock ordered as
    #  | Kvv Kvp | 
    #  | Kpv Kpp |3x3
    
##    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type
##    include_dynamics = True

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 2

    def GetVectorValueOnGauss(self, var_x, var_y, N,step=0):
        value = zeros(2)
        for i in range(0, self.geometry.GetNumberOfNodes() ):
            value[0] += N[i] * self.geometry[i].GetSolutionStepValue(var_x, step)
            value[1] += N[i] * self.geometry[i].GetSolutionStepValue(var_y, step)
        return value
                
    def CalculateLocalSystem(self, im_useless):
        nnodes = self.geometry.GetNumberOfNodes()
        dofs_per_node = self.GetDofsPerNode()
        mat_size = nnodes*dofs_per_node # mat_size = 4
        
        [Ns, derivatives, weights] = self.geometry.ShapeFunctions()
        number_of_gauss = len(Ns)
        
        RHS = zeros(mat_size)  # no external forces so far
        LHS = zeros((mat_size,mat_size))

        # define variables of the method
        U_0 = 1 # TO BE CHANGED: computed as local characteristic velocity
        delta = 0.01 # a small positive number

        x1 = self.geometry[0].coordinates[0]
        y1 = self.geometry[0].coordinates[1]
        x2 = self.geometry[1].coordinates[0]
        y2 = self.geometry[1].coordinates[1]

        A= sqrt((x1 - x2)**2 + (y1 - y2)**2)
        
        #integrate external forces to the RHS
        for gauss in range(0, number_of_gauss):
            weight = weights[gauss]
            N = Ns[gauss]

            # last step: in k-step of Newton-Raphson method, the k-1 values are
            # needed (so the last ones since the k-step values are to be
            # computed).
            step = -1 # -1 index => last value
            
            # k-1-step values of velocity at the current Gauss point
            [u, v] = self.GetVectorValueOnGauss(VELOCITY_X, VELOCITY_Y, N, step)

            # squared module of velocity at the current Gauss point
            squared_velocity_module = u**2 +v**2

            # normal vector
            n = [0, 1]

            # u_vector * n (dot product)
            projected_vel = u * n[0] + v * n[1] 

            # the S_0 function 
            S=0.5*(1-tanh( projected_vel / (U_0*delta) ) )

            
            RHS[0] += 0.5*squared_velocity_module*S*n[0]*N[0]*A # u1
            RHS[1] += 0.5*squared_velocity_module*S*n[1]*N[0]*A # v1
            RHS[2] += 0.5*squared_velocity_module*S*n[0]*N[1]*A # u2
            RHS[3] += 0.5*squared_velocity_module*S*n[1]*N[1]*A # v2

        return [LHS, RHS]
##        ???
        return C

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], VELOCITY_X))
        unknowns.append(Dof(self.geometry[0], VELOCITY_Y))
        unknowns.append(Dof(self.geometry[1], VELOCITY_X))
        unknowns.append(Dof(self.geometry[1], VELOCITY_Y))
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[0].EquationId(VELOCITY_Y))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_X))
        equation_ids.append(self.geometry[1].EquationId(VELOCITY_Y))
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(VELOCITY_X, step)
        values[1] = self.geometry[0].GetSolutionStepValue(VELOCITY_Y, step)
        values[2] = self.geometry[1].GetSolutionStepValue(VELOCITY_X, step)
        values[3] = self.geometry[1].GetSolutionStepValue(VELOCITY_Y, step)
        return values


