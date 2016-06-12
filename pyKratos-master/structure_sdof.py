#===============================================================================
# this is SDOF system solved by direct time integration - Generalized alpha scheme
# implemented following the: Andre.M (2012), Generalized alpha metod for LAGRANGE
# 
# this is a currently implemented and used version that stores u0,v0,a0,f0 and u1,v1,a1,f1
#
# Ivan Hanzlicek, 2014
#===============================================================================
 
# import
from pylab import *
from math import sqrt, fabs

# class to define the structure and Generalized alpha scheme parameter (fixed time step)
class StructureSDOF:


    """ Direct time integration scheme of linear damped SDOF system (Generalized alpha with fixed dt)"""
    
    def __init__(self, dt, M, B, K, pInf, u0, v0, a0):
        self.dt = dt

        # structure mass and spring stiffness
        self.M = M
        self.B = B
        self.K = K

        # generalized alpha parameters (to ensure unconditional stability, 2nd order accuracy)
        self.alphaM = (2.0 * pInf - 1.0) / (pInf + 1.0)
        self.alphaF = pInf / (pInf + 1.0)
        self.beta = 0.25 * (1 - self.alphaM + self.alphaF) ** 2
        self.gamma = 0.5 - self.alphaM + self.alphaF

        # coefficients for LHS
        self.a1h = (1.0 - self.alphaM) / (self.beta * self.dt ** 2)
        self.a2h = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a3h = 1.0 - self.alphaF

        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.dt
        self.a3m = (1.0 - self.alphaM - 2.0 * self.beta) / (2.0 * self.beta)

        # coefficients for damping
        self.a1b = self.a2h
        self.a2b = ((1.0 - self.alphaM) * self.gamma - self.beta) / self.beta 
        self.a3b = (self.gamma - 2.0 * self.beta) * (1.0 - self.alphaF) * self.dt / (2.0 * self.beta)

        # coefficient for stiffness
        self.a1k = -1.0 * self.alphaF

        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.dt)
        self.a2v = 1.0 - self.gamma / self.beta

        #self.a3v = (1.0 - self.gamma / (2 * self.beta))
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.dt # CHECK

        # coefficients for acceleration update
        # self.a1a = self.a1v / self.dt
        self.a1a = self.a1v / (self.dt * self.gamma) # CHECK
        self.a2a = -1.0 / (self.beta * self.dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)   
        
        # initial displacement, velocity and acceleration
        self.u0 = u0
        self.v0 = v0
        self.a0 = a0
        self.u1 = u0
        self.v1 = v0
        self.a1 = a0
        
        # force from a previous time step (initial force)
        self.f0 = 0.0
        self.f1 = 0.0
        
    def printSetUp(self):
        print("Printing integration scheme set up:")
        print("dt: ", self.dt)
        print("alphaM: ", self.alphaF)
        print("alphaF: ", self.alphaM)
        print("gamma: ", self.gamma)
        print("betta: ", self.beta)
        print("Printing structural set up:")
        print("Structural Mass: ", self.M)
        print("Structural Stiffness: ", self.K)
        print("Structural Damping: ", self.B)
        print("Printing initial conditions:")
        print("u0: ", self.u0)
        print("v0: ", self.v0)
        print("a0: ", self.a0)
        
        
    def printValuesAtCurrentStep(self, n):
        print("Printing values at step no: ", n)
        print("u0: ", self.u1)
        print("v0: ", self.v1)
        print("a0: ", self.a1)
        print("f0: ", self.f1)

    def getDisplacement(self):
        return self.u1

    def getVelocity(self):
        return self.v1
    
    def getAcceleration(self):
        return self.a1
                   
    def solveStructure(self, f1):
        # sys of eq reads: LHS * u1 = RHS
        F = (1.0 - self.alphaF) * f1 + self.alphaF * self.f0 
        LHS = self.a1h * self.M + self.a2h * self.B + self.a3h * self.K
        RHS = self.M * (self.a1m * self.u0 + self.a2m * self.v0 + self.a3m * self.a0)
        RHS += self.B * (self.a1b * self.u0 + self.a2b * self.v0 + self.a3b * self.a0)
        RHS += self.a1k * self.K * self.u0 + F
        
        # update self.f1
        self.f1 = f1
        
        # updates self.u1,v1,a1
        self.u1 = RHS / LHS
        self.v1 = self.a1v * (self.u1 - self.u0) + self.a2v * self.v0 + self.a3v * self.a0
        self.a1 = self.a1a * (self.u1 - self.u0) + self.a2a * self.v0 + self.a3a * self.a0
        
    def updateStructureTimeStep(self):    
        # update displacement, velocity and acceleration 
        self.u0 = self.u1
        self.v0 = self.v1
        self.a0 = self.a1
       
        # update the force   
        self.f0 = self.f1
        
    def checkConvergence(self, a, b, tol):
        print("a = ", a)
        print("b = ", b)
        if fabs(a-b) < tol:
            print("e1 norm |a-b| = ", fabs(a-b), " is smaller than tolerance ", tol)
            return 1
        else:
            return 0   
