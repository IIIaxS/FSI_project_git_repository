import numpy as np
from scipy.fftpack import fft
import matplotlib.pyplot as plt

# import data
data =  "forceAndDisplacement0_01.dat"
# read the columns
time, disp, lift, rot, mom = np.genfromtxt(data, skip_header = 0 , unpack = True)

# number of sample points
N = len(time)
# sample spacing
T = time[1] - time[0]
x = np.linspace(0.0, N*T, N)
fftlift = fft(lift)
fftmom = fft(mom)
ffttime = np.linspace(0.0, 1.0/(2.0*T), N/2)


# plot in time domain
plt.subplot(4, 1, 1)
plt.plot(time, lift)
plt.title("signal in time and frequency domain \n $ u_{in} = 30 $ m/s")
plt.ylabel("lift")
plt.xlabel("time in $ s $")
plt.grid(True)

plt.subplot(4, 1, 2)
plt.plot(time, mom)
plt.ylabel("moment")
plt.xlabel("time in $ s $")
plt.grid(True)

# plot in frequency domain
plt.subplot(4 , 1, 3)
plt.plot(ffttime, 2.0/N * np.abs(fftlift[0:N/2]))
#plt.title("signal in frequency domain")
plt.xlabel("frequency in Hz")
plt.ylabel("lift")
plt.grid(True)

plt.subplot(4, 1, 4)
plt.plot(ffttime, 2.0/N * np.abs(fftmom[0:N/2]))
plt.xlabel("frequency in Hz")
plt.ylabel("moment")
plt.grid(True)

plt.show()



