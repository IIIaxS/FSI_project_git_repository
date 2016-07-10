# ----------------------------------------------------------------------
# author  : Dennis Kasper
# contact : dennis.kasper@tum.de
#
# code for plotting results
# ----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# read in the file header parameters
resultFile01 = open("resultFile101.dat", "r")
para01 = resultFile01.readline().split()
paraValues01 = resultFile01.readline().split()
resultFile01.close()

resultFile02 = open("resultFile102.dat", "r")
para02 = resultFile02.readline().split()
paraValues02 = resultFile02.readline().split()
resultFile02.close()

resultFile03 = open("resultFile103.dat", "r")
para03 = resultFile03.readline().split()
paraValues03 = resultFile03.readline().split()
resultFile03.close()

resultFile04 = open("resultFile104.dat", "r")
para04 = resultFile04.readline().split()
paraValues04 = resultFile04.readline().split()
resultFile04.close()


#read in the data
time01, disp01 = np.genfromtxt("resultFile101.dat", skiprows = 3, usecols = (0, 1),unpack = True)
time02, disp02 = np.genfromtxt("resultFile102.dat", skiprows = 3, usecols = (0, 1),unpack = True)
time03, disp03 = np.genfromtxt("resultFile103.dat", skiprows = 3, usecols = (0, 1),unpack = True)
time04, disp04 = np.genfromtxt("resultFile104.dat", skiprows = 3, usecols = (0, 1),unpack = True)


import matplotlib.pyplot as plt

plt.figure(1)
plt.title("Structural solver $ \Delta t $ study, M=%s, K=%s, B=%s, f=%s \n $ F(t) = 10 cos(2 \pi f t) $" % (round(float(paraValues01[1]),1), round(float(paraValues01[2]),1), round(float(paraValues01[3]),1), round(float(paraValues01[5]),3)))
plt.plot(time01, disp01, label = " $ \Delta t = 0.1 $")
plt.plot(time02, disp02, label = " $ \Delta t = 0.01 $")
plt.plot(time03, disp03, label = " $ \Delta t = 0.001 $")
plt.plot(time04, disp04, label = " $ \Delta t = 0.0001 $")
plt.legend(loc = "upper left")
plt.xlabel("time t in s")
plt.ylabel("displacement v")
plt.show()