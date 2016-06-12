#===============================================================================
# aerodynamic evaluation class
# 
# Mate Pentek, 2015 
#===============================================================================
 
# import
from pylab import *
from math import sqrt, fabs
#from sympy import coef
import sympy
from numpy import fft
from pylab import xticks

def cm2inch(value):#(*tupl):
    inch = 2.54
    #if type(tupl[0]) == tuple:
    #    return tuple(i/inch for i in tupl[0])
    #else:
    #    return tuple(i/inch for i in tupl)
    return value/inch

#from file 1
simulTime = loadtxt("displacement.dat", skiprows=1, usecols = (0,))
#change sign of forces - due to definition of results from Kratos
displacement = loadtxt("displacement.dat", skiprows=0, usecols = (1,)) #drag - stramwise -x
struct_react = loadtxt("displacement.dat", skiprows=0, usecols = (2,)) #drag - stramwise -x
fluid_react = multiply(loadtxt("structure_drag.dat", skiprows=0, usecols = (2,)),-1) #drag - stramwise -x
rotation = loadtxt("rotation.dat", skiprows=0, usecols = (1,)) #drag - stramwise -x
struct_mom = loadtxt("rotation.dat", skiprows=0, usecols = (2,)) #drag - stramwise -x
fluid_mom = multiply(loadtxt("structure_drag.dat", skiprows=0, usecols = (6,)),-1) #drag - stramwise -x

#plot forces
# row and column sharing
f1, ((ax1, ax2)) = plt.subplots(2, 1)#, sharex='col', sharey='row')
f1.set_size_inches(cm2inch(16),cm2inch(8))
title = "Deck movement"
f1.suptitle(title, fontsize=12)


title = "Displacement"
ax1.set_title(title, fontsize="12")
ax1.set_ylabel("Displacement $[m]$", fontsize="12")
ax1.plot(simulTime, displacement, "-r", lw=1.0)
ax1.legend(loc='best', fontsize="12")
ax1.set_xlim(0, 6)

title = "Rotation"
ax2.set_title(title, fontsize="12")
ax2.set_ylabel("Rotation $[rad]$", fontsize="12")
ax2.plot(simulTime, rotation, "-r", lw=1.0)
ax2.set_xlim(0, 6)

# save data (optionally you can change .pdf to .svg)
outputName = "Wt " + "Movement" + ".pdf"
savefig(outputName + "1", format="pdf", dpi=200)
print("File created: " + outputName)    

#plot forces
# row and column sharing
f2, ((ax1, ax2)) = plt.subplots(2, 1)#, sharex='col', sharey='row')
f1.set_size_inches(cm2inch(16),cm2inch(8))
title = "Deck forces"
f1.suptitle(title, fontsize=12)


title = "Ry  - action (fluid) and reaction (structure)"
ax1.set_title(title, fontsize="12")
ax1.set_ylabel("Force $[N]$", fontsize="12")
ax1.plot(simulTime, struct_react,"-r", label= "struct", lw=1.0)
ax1.plot(simulTime, fluid_react, "-k", label= "fluid", lw=1.0)
ax1.legend(loc='best', fontsize="12")
ax1.set_xlim(0, 6)

title = "Mz - action (fluid) and reaction (structure)"
ax2.set_title(title, fontsize="12")
ax2.set_ylabel("Moment $[Nm]$", fontsize="12")
ax2.plot(simulTime, struct_mom,"-r", label= "struct", lw=1.0)
ax2.plot(simulTime, fluid_mom, "-k", label= "fluid", lw=1.0)
ax2.legend(loc='best', fontsize="12")
ax2.set_xlim(0, 6)

# save data (optionally you can change .pdf to .svg)
outputName = "Wt " + "Forces" + ".pdf"
savefig(outputName + "2", format="pdf", dpi=200)
print("File created: " + outputName)    

