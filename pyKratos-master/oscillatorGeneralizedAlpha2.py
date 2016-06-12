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
    
    """ Direct time integration of linear undamped SDOF (Generalized alpha with fixed dt)"""
    
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
        #
        #
        # coefficients for LHS
        self.a1h = (1.0 - self.alphaM) / (self.beta * self.dt ** 2)
        self.a2h = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a3h = 1.0 - self.alphaF
        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.dt
        self.a3m = (1.0 - self.alphaM - 2.0 * self.beta) / (2.0 * self.beta)
        #coefficients for damping
        self.a1b = self.a2h
        self.a2b = ((1.0 - self.alphaM) * self.gamma - self.beta) / self.beta 
        self.a3b = (self.gamma - 2.0 * self.beta) * (1.0 - self.alphaF) * self.dt / (2.0 * self.beta)
        # coefficient for stiffness
        self.a1k = -1.0 * self.alphaF
        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.dt)
        self.a2v = 1.0 - self.gamma / self.beta
#         self.a3v = (1.0 - self.gamma / (2 * self.beta))
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.dt # CHECK
        # coefficients for acceleration update
#         self.a1a = self.a1v / self.dt
        self.a1a = self.a1v / (self.dt * self.gamma) # CHECK
        self.a2a = -1.0 / (self.beta * self.dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)   
        #
        # initial displacement, velocity and acceleration
        self.u0 = u0
        self.v0 = v0
        self.a0 = a0
        #
        self.u1 = u0
        self.v1 = v0
        self.a1 = a0
        # force from a previous time step (initial force)
        self.f0 = 0.0
        self.f1 = 0.0
        
    def printSetUp(self):
        print("Printing integration scheme set up:")
        print()
        print("dt: ", self.dt)
        print("alphaM: ", self.alphaF)
        print("alphaF: ", self.alphaM)
        print("gamma: ", self.gamma)
        print("betta: ", self.beta)
        print()
        print("Printing structural set up:")
        print()
        print("Structural Stiffness: ", self.K)
        print("Structural Mass: ", self.M)
        
    def printValuesAtCurrentStep(self, n):
        print("Printing values at step no: ", n)
        print("u0: ", self.u1)
        print("v0: ", self.v1)
        print("a0: ", self.a1)
        print("f0: ", self.f1)
        print()

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
        #
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
                
# END of class StructureSDOF:

#===============================================================================
# demo        
if __name__ == "__main__":

    tEnd = 10.0

    # steps = sampling frequency
    n = 10000
    dt = tEnd / n

#     # 10 Hz excitation force - for cosine function
#     freq = 10
#     Amp = 1

    currentTime = 0.0
 
    M = 0.183
#     #for resonance try
#     M = 1
#     
    K = 652.01
#     #for resonance try
#     K = 4.0 * math.pi ** 2 * 100.0
    
    xi = 0.015
    B = xi * 2.0 * sqrt(M * K)
#     #for resonance try
#     B = 0.0
    
    # 1.0 introduces numerical error - numerical oscillations
    pInf = 0.9 #1.0

    # initial displacement from measured results 
    u0 = 0.004
#     #for resonance try
#     u0 = 0.0
    
    v0 = 0.0
    a0 = 0.0
 
    # instantiate an object: structure
    structure = StructureSDOF(dt, M, B, K, pInf, u0, v0, a0)
    structure.printSetUp()
 
    # lists to store the results
    time = []
    disp = []
    acc  = []
    veloc  = []

    force = []
    
    # clean up the result file
    open("results.dat", "w").close()

    # file header
    file1 = open("results.dat", "a")
    file1.write("Time" + "\t")
    file1.write("Disp" + "\t")
    file1.write("Vel" + "\t")
    file1.write("Acc" + "\t")
    file1.write("Force" + "\t")
    file1.write("\n")
    file1.close()

    # initial values
    time.append(currentTime)
    disp.append(u0)
    veloc.append(v0)
    acc.append(a0)

    f = 0.0
#     #for resonance try
#     f = Amp * cos(2* math.pi * freq * currentTime)
    force.append(f)

    file1 = open("results.dat", "a")
    file1.write(str(currentTime) + "\t")
    file1.write(str(u0) + "\t")
    file1.write(str(v0) + "\t")
    file1.write(str(a0) + "\t")
    file1.write(str(f) + "\t")
    file1.write("\n")
    file1.close()
    
    # computation loop
    for n in range(1, n):
        currentTime = currentTime + dt
        
        f = 0.0
#         # for resonance try
#         # cosine excitation of 10 Hz
#         f = Amp * cos(2* math.pi * freq * currentTime)
        
        force.append(f)
        
        # solve the problem    
        structure.solveStructure(f)
        #  update results
        structure.updateStructureTimeStep()

        time.append(currentTime)
        disp.append(structure.getDisplacement())
        veloc.append(structure.getVelocity())
        acc.append(structure.getAcceleration())

        # append the results to the currentresult.dat
        file1 = open("results.dat", "a")
        file1.write(str(currentTime) + "\t")
        file1.write(str(structure.getDisplacement()) + "\t")
        file1.write(str(structure.getVelocity()) + "\t")
        file1.write(str(structure.getAcceleration()) + "\t")
        file1.write(str(f) + "\t")
        file1.write("\n")
        file1.close()
    
    # compute analytical solution of undamped SDOF with initial displacement
    omega = sqrt( K / M)
    print("Analytical: circular natural frequency: ", omega)
    print("Analytical: natural period: ", 2.0 * math.pi / omega)
     
    # plot out the results in the file
    plt.figure(1)
    plt.plot(time, disp, "-k", label="$\mathrm{REF}$", lw=0.5)
    plt.xlim(0.0, 10.0 )#time[-1])
    #plt.ylim(-2.0, 2.0)
    plt.xlabel("time")
    plt.ylabel("displacement")
    savefig("displacement.pdf", format="pdf")

    plt.figure(2)
    plt.plot(time, veloc, "-k", label="$\mathrm{REF}$", lw=0.5)
    plt.xlim(0.0, 10.0 )#time[-1])
    #plt.ylim(-2.0, 2.0)
    plt.xlabel("time")
    plt.ylabel("velocity")
    savefig("velocity.pdf", format="pdf")

    plt.figure(3)
    plt.plot(time, acc, "-k", label="$\mathrm{REF}$", lw=0.5)
    plt.xlim(0.0, 10.0 )#time[-1])
    #plt.ylim(-2.0, 2.0)
    plt.xlabel("time")
    plt.ylabel("acceleration")
    savefig("acceleration.pdf", format="pdf")

    plt.figure(4)
    plt.plot(time, force, "-k", label="$\mathrm{REF}$", lw=0.5)
    plt.xlim(0.0, 10.0 )#time[-1])
    #plt.ylim(-2.0, 2.0)
    plt.xlabel("time")
    plt.ylabel("force")
    savefig("force.pdf", format="pdf")
    
    print("File created")
     
    show()
#===============================================================================        
