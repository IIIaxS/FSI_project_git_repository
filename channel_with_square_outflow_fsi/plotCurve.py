# plotCurve.py
# plot the displacement and force history

import matplotlib.pyplot as plt
import numpy as np

# read data from "forceAndDisplacement.dat"
data = "forceAndDisplacement0_01.dat"
t, disp, lift, rot, mom = np.genfromtxt(data, skip_header = 0, unpack = True)

# plot vertical motion
plt.subplot(2, 1, 1)             # the first subplot in the first figure
plt.plot(t, lift)
#plt.title('Vertical motion')
#plt.xlabel('time in $ s $')
plt.ylabel('lift force')
plt.grid(True)

plt.subplot(2, 1, 2)             # the second subplot in the first figure
plt.plot(t, mom)
plt.xlabel('time in $ s $')
plt.ylabel('moment')
plt.grid(True)

plt.show()
#plt.savefig('Translation.pdf')