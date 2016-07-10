# import
from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
#print(sys.path)
from numpy import *
from pyKratos import *
from pylab import *
from math import sqrt, fabs
from structure_sdof import *


# begin testing
#tEnd = 10.0

# steps = sampling frequency
#n = 10000
#dt = tEnd/n

# DENNIS BEGIN
# redefinition of time integration parameter to be consistent with the fluid solver
dt = 0.01
n = 1000
tEnd = n*dt

# name the result file
nameResultFileRes = "resultFileResDampTest.dat"
# DENNIS END

currentTime = 0.0

M = 1

# frequency
freq = 2.0

K = (2*pi*freq)**2 * M
#for resonance try
#K = 4.0 * math.pi ** 2 * 100.0

xi = 0.1 # x% damping
B = xi * 2.0 * sqrt(M * K)
#for resonance try
#B = 0.0

#omega0 = sqrt(K/M)
#print('omega_0 = ', omega0)
#omega = omega0*sqrt(1-2*xi**2)
#print('omega = ', omega)
#resFreq = omega0/(2*pi)
#print('resonance frequency = ', resFreq)

# excitation force - for cosine function
#freq = resFreq
Amp = 1.0

# 1.0 introduces numerical error - numerical oscillations
pInf = 1.0 #1.0

# initial displacement from measured results 
#u0 = 0.004
#for resonance try
u0 = 0.0
v0 = 0.0
a0 = 0.0
 
# instantiate an object: sdofSystem
sDofSystem = StructureSDOF(dt, M, B, K, pInf, u0, v0, a0)
#sDofSystem.printSetUp()
 
# lists to store the results
time = []
disp = []
acc  = []
veloc  = []
force = []


# clean up the result file
open(nameResultFileRes, "w").close()

# file header
file1 = open(nameResultFileRes, "a")
file1.write("dt" + "\t")
file1.write("M" + "\t")
file1.write("K" + "\t")
file1.write("B" + "\t")
file1.write("Amp" + "\t")
file1.write("freq" + "\t")
file1.write("pInf" + "\t")
file1.write("\n")
file1.write(str(dt) + "\t")
file1.write(str(M) + "\t")
file1.write(str(K) + "\t")
file1.write(str(B) + "\t")
file1.write(str(Amp) + "\t")
file1.write(str(freq) + "\t")
file1.write(str(pInf) + "\t")
file1.write("\n")
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

#for resonance try
f = Amp * sin(2 * pi * freq * currentTime)
force.append(f)

file1 = open(nameResultFileRes, "a")
file1.write(str(currentTime) + "\t")
file1.write(str(u0) + "\t")
file1.write(str(v0) + "\t")
file1.write(str(a0) + "\t")
file1.write(str(f) + "\t")
file1.write("\n")
file1.close()

# computation loop
for n in range(1, n):
    currentTime += dt
    
    #f = 0.0
    # for resonance try
    # sin excitation of 10 Hz
    f = Amp*sin(2*pi*freq*currentTime)
    
    force.append(f)
    
    # solve the problem    
    sDofSystem.solveStructure(f)
    #  update results
    sDofSystem.updateStructureTimeStep()

    time.append(currentTime)
    disp.append(sDofSystem.getDisplacement())
    veloc.append(sDofSystem.getVelocity())
    acc.append(sDofSystem.getAcceleration())

    # append the results to the currentresult.dat
    file1 = open(nameResultFileRes, "a")
    file1.write(str(currentTime) + "\t")
    file1.write(str(sDofSystem.getDisplacement()) + "\t")
    file1.write(str(sDofSystem.getVelocity()) + "\t")
    file1.write(str(sDofSystem.getAcceleration()) + "\t")
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
plt.xlim(0.0, tEnd )#time[-1])
#plt.ylim(-2.0, 2.0)
plt.title('$dt=%s, M=%s, K=%s, B=%s, xi=%s, u_0=%s, v_0=%s, a_0=%s, Amp=%s, f=%s, pInf=%s $ \n $ f(t) = Amp*sin(2*pi*f*t) $' % (dt, round(M,1), round(K,1), round(B,1), round(xi,2), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3), round(pInf,3)))
plt.xlabel("time $ t $ in $ s $")
plt.ylabel("displacement $u$")
#savefig("displDt=%stEnd=%sM=%sK=%sB=%su0=%sv0=%sa0=%sAmp=%sf=%s.pdf" % (dt, round(tEnd,2), round(M,1), round(K,1), round(B,1), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3)), format="pdf")

# plt.figure(2)
# plt.plot(time, veloc, "-k", label="$\mathrm{REF}$", lw=0.5)
# plt.xlim(0.0, tEnd)#time[-1])
# #plt.ylim(-2.0, 2.0)
# plt.title('$dt=%s, M=%s, K=%s, B=%s, u_0=%s, v_0=%s, a_0=%s, Amp=%s, f=%s $\n $ f(t) = Amp*sin(2*pi*f*t) $' % (dt, round(M,1), round(K,1), round(B,1), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3)))
# plt.xlabel("time t in s")
# plt.ylabel("velocity")
# #savefig("velDt=%stEnd=%sM=%sK=%sB=%su0=%sv0=%sa0=%sAmp=%sf=%s.pdf" % (dt, round(tEnd,2), round(M,1), round(K,1), round(B,1), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3)), format="pdf")
#
# plt.figure(3)
# plt.plot(time, acc, "-k", label="$\mathrm{REF}$", lw=0.5)
# plt.xlim(0.0, tEnd)#time[-1])
# #plt.ylim(-2.0, 2.0)
# plt.title('$dt=%s, M=%s, K=%s, B=%s, u_0=%s, v_0=%s, a_0=%s, Amp=%s, f=%s $\n $ f(t) = Amp*sin(2*pi*f*t) $' % (dt, round(M,1), round(K,1), round(B,1), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3)))
# plt.xlabel("time t in s")
# plt.ylabel("acceleration")
# #savefig("accDt=%stEnd=%sM=%sK=%sB=%su0=%sv0=%sa0=%sAmp=%sf=%s.pdf" % (dt, round(tEnd,2), round(M,1), round(K,1), round(B,1), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3)), format="pdf")
#
# plt.figure(4)
# plt.plot(time, force, "-k", label="$\mathrm{REF}$", lw=0.5)
# plt.xlim(0.0, tEnd)#time[-1])
# #plt.ylim(-2.0, 2.0)
# plt.title('$dt=%s, M=%s, K=%s, B=%s, u_0=%s, v_0=%s, a_0=%s, Amp=%s, f=%s $\n $ f(t) = Amp*sin(2*pi*f*t) $' % (dt, round(M,1), round(K,1), round(B,1), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3)))
# plt.xlabel(" t in s")
# plt.ylabel("force")
#savefig("forceDt=%stEnd=%sM=%sK=%sB=%su0=%sv0=%sa0=%sAmp=%sf=%s.pdf" % (dt, round(tEnd,2), round(M,1), round(K,1), round(B,1), round(u0,1), round(v0,1) , round(a0,1), round(Amp,1), round(freq,3)), format="pdf")

#print("File created")
 
show()
