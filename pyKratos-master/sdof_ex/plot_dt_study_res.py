# ----------------------------------------------------------------------
# author  : Dennis Kasper
# contact : dennis.kasper@tum.de
#
# code for plotting results
# ----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# read in the file header parameters
resultFileRes01 = open("resultFileRes501.dat", "r")
para01 = resultFileRes01.readline().split()
paraValues01 = resultFileRes01.readline().split()
resultFileRes01.close()
print "Parameter File 01"
print para01
print paraValues01

resultFileRes02 = open("resultFileRes502.dat", "r")
para02 = resultFileRes02.readline().split()
paraValues02 = resultFileRes02.readline().split()
resultFileRes02.close()
print "Parameter File 02"
print para02
print paraValues02

resultFileRes03 = open("resultFileRes503.dat", "r")
para03 = resultFileRes03.readline().split()
paraValues03 = resultFileRes03.readline().split()
resultFileRes03.close()
print "Parameter File 03"
print para03
print paraValues03

resultFileRes04 = open("resultFileRes504.dat", "r")
para04 = resultFileRes04.readline().split()
paraValues04 = resultFileRes04.readline().split()
resultFileRes04.close()
print "Parameter File 04"
print para04
print paraValues04

# resultFileRes05 = open("resultFileRes105.dat", "r")
# para05 = resultFileRes05.readline().split()
# paraValues05 = resultFileRes05.readline().split()
# resultFileRes05.close()
# print "Parameter File 05"
# print para05
# print paraValues05

# resultFileRes06 = open("resultFileRes06.dat", "r")
# para06 = resultFileRes06.readline().split()
# paraValues06 = resultFileRes06.readline().split()
# resultFileRes06.close()
# print "Parameter File 06"
# print para06
# print paraValues06


#read in the data
time01, disp01 = np.genfromtxt("resultFileRes501.dat", skip_header = 3, usecols = (0, 1),unpack = True)
time02, disp02 = np.genfromtxt("resultFileRes502.dat", skip_header = 3, usecols = (0, 1),unpack = True)
time03, disp03 = np.genfromtxt("resultFileRes503.dat", skip_header = 3, usecols = (0, 1),unpack = True)
time04, disp04 = np.genfromtxt("resultFileRes504.dat", skip_header = 3, usecols = (0, 1),unpack = True)
#time05, disp05 = np.genfromtxt("resultFileRes105.dat", skiprows = 3, usecols = (0, 1),unpack = True)
#time06, disp06 = np.genfromtxt("resultFileRes06.dat", skiprows = 3, usecols = (0, 1),unpack = True)

# compute analytic solution
import math as math
dt = 0.0001
tEnd = 10
nsteps = int(tEnd/dt)
dispAna = []
timeAna = []
time = 0
timeAna.append(time)
dispAna.append(0)
F0 = 1.0
m = float(paraValues01[1])
print m
freq = float(paraValues01[5])
print freq
#k = (2*np.pi*f)**2 * m
k = float(paraValues01[2])
print k
Omega = math.sqrt(k/m)

for i in range(0, nsteps):
    time += dt
    disp = F0/(2*k)*(math.sin(Omega*time) - Omega*time*math.cos(Omega*time))
    timeAna.append(time)
    dispAna.append(disp)


plt.figure(1)
plt.title("$ F(t) = 1 * sin(2 \pi f t) $ \n Structural solver $ \Delta t $ study, M=%s, K=%s, B=%s, f=%s \n" % (round(float(paraValues01[1]),1), round(float(paraValues01[2]),1), round(float(paraValues01[3]),1), round(float(paraValues01[5]),3)))
plt.plot(timeAna, dispAna, 'k-', label = " ana sol")
plt.plot(time01, disp01, 'r--', label = " $ \Delta t = 0.1 $")
plt.plot(time02, disp02, 'g:', label = " $ \Delta t = 0.05 $")
plt.plot(time03, disp03, 'm-.', label = " $ \Delta t = 0.01 $")
plt.plot(time04, disp04, 'c', label = " $ \Delta t = 0.005 $")
#plt.plot(time05, disp05, 'r:', label = " $ \Delta t = 0.001 $")
#plt.plot(time04, disp04, 'g', label = " $ \Delta t = 0.0001 $")
plt.legend(loc = "upper left")
plt.grid(True)
plt.xlabel("time t in s")
plt.ylabel("displacement v")
plt.show()