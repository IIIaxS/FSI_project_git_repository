# ----------------------------------------------------------------------
# author  : Martin Ruchti
# contact : martin.ruchti@tum.de
# ----------------------------------------------------------------------
from __future__ import print_function, absolute_import, division 
from numpy import *
from .variables import *

'''
    this scheme is an implementation of the newton raphson scheme
    
    
    the scheme is created as a modified version of 'static_scheme.py' and follows the implementation of the
    residual based newton raphson strategy in KratosMultiphysics
'''
class NewtonRaphsonStrategy:

    def __init__(self, model_part, scheme, builder_and_solver, max_iteration, epsAbs, ratioRel = 1e-10):
        self.model_part             = model_part
        self.scheme                 = scheme
        self.builder_and_solver     = builder_and_solver
        #self.adjoint_builder_and_solver = adjointbuilder_and_solv
        self.max_iteration_number   = max_iteration
        self.epsilon                = epsAbs
        self.relativeRatio          = ratioRel
        
        #allocate matrices
        self.A                      = zeros((0, 0))
        self.b                      = zeros((0))
        self.dx                     = zeros((0))
        
        #file
        self.file                   = 0

    def Initialize(self):
        # find list of dofs
        self.builder_and_solver.SetupDofSet()

        # allocate memory for the system
        self.A, self.x, self.b = self.builder_and_solver.SetupSystem(self.A, self.dx, self.b)

    def Solve(self):
        print("=================================================================")
        print("start fluid solving process ...")
                
        self.Initialize()
        
        # do prediction once
        self.scheme.Predict()
        
        # initialize parameters for NR - strategy
        iteration_number = 1
        is_converged = False
        error_L2_norm = 0.0
        init_L2_norm = 0.0
        ratio = 1.0
        dragLift = [0, 0]
        # Dennis begin
        moment = 0.0
        # Dennis end

        # first solve
        self.A, self.dx, self.b = self.builder_and_solver.BuildAndSolve(self.A, self.x, self.b)
                
        # check for convergence
        error_L2_norm = 0.0
        for i in range(0,len(self.dx)):
            error_L2_norm += (self.dx[i])**2
                
        #scale error with number of nodes
        error_L2_norm = sqrt(error_L2_norm)/len(self.dx)
        init_L2_norm = error_L2_norm
        
        if( error_L2_norm <= self.epsilon ):
            is_converged = True
            print("coverged after step: ",iteration_number)
            print("error is: ",error_L2_norm)
            moment = self.builder_and_solver.ComputeReactions(self.A, self.x, self.b, dragLift)
        else:
            print("not converged, error is: ",error_L2_norm)
            print("ratio is: ", 1)
            print("-----------------------------------------------------------------")
            
        # call scheme to do update
        self.scheme.Update(self.builder_and_solver.dofset, self.dx)
        
        #iterate if not converged
        while(not is_converged and iteration_number < self.max_iteration_number):

            # do build and solve
            self.A, self.dx, self.b = self.builder_and_solver.BuildAndSolve(self.A, self.x, self.b)
            
            # call scheme to do update
            self.scheme.Update(self.builder_and_solver.dofset, self.dx)
            
            #check for convergence
            error_L2_norm = 0.0
            for i in range(0,len(self.dx)):
                error_L2_norm += (self.dx[i])**2
                
            #scale error with number of nodes
            error_L2_norm   = sqrt(error_L2_norm)/len(self.dx)
            
            #compute relative error
            ratio = error_L2_norm / init_L2_norm
            
            if( error_L2_norm <= self.epsilon or ratio <= self.relativeRatio ):
                is_converged = True
            else:
                print("not converged, error is: ",error_L2_norm)
                print("ratio is: ",ratio)
                print("-----------------------------------------------------------------")
            
            iteration_number += 1
        if(iteration_number == self.max_iteration_number):
            print("*********maximum iterations reached*********")
            print("error is: ",error_L2_norm)
            print("ratio is: ",ratio)
            moment = self.builder_and_solver.ComputeReactions(self.A, self.x, self.b, dragLift)
            print("solving fluid process done!")
            print("=================================================================")
        elif(iteration_number > 1):
            print("coverged after step: ",iteration_number)
            print("error is: ",error_L2_norm)
            print("ratio is: ",ratio)
            moment = self.builder_and_solver.ComputeReactions(self.A, self.x, self.b, dragLift)
            print("solving fluid process done!")
            print("=================================================================")
        
        if(self.file != 0):
            self.WriteDragforceToFile(dragLift, moment)

        
    def SpyMatrix(self):
        try:
            import matplotlib.pylab as pl
            pl.spy(self.A)
            pl.show()
        except:
            raise Exception(
                "error in function Spy. Probably matplotlib not installed")

    
    def WriteDragforceToFile(self, dragLift, moment):
        # dennis begin
        print('WriteDragforceToFile has been called')
        print('moment = ', moment)
        print('drag lift =', dragLift)
        # dennis end
        time = self.model_part.ProcessInfo[TIME]
        
        #output 20 digits
        output = str(time) + " " + str.format("{0:.20f}", dragLift[0]) + " " + str.format("{0:.20f}", dragLift[1]) + " " + str.format("{0:.20f}", moment) + "\n"
            
        self.WriteToFile(output)

    #file operations
    def OpenFile(self, filename):
        self.file   = open(filename, "w")
        
    def WriteToFile(self, data):
        self.file.write(data)
        self.file.flush()
        
    def CloseFile(self):
        self.file.close()