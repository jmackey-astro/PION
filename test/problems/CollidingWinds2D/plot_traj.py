

import numpy as np
import matplotlib.pyplot as plt
#data=np.loadtxt("trajectory.txt")
data=np.genfromtxt("trajectory.txt", skip_footer=1)
x1 = data[:,2]
y1 = data[:,3]

x2 = data[:,6]
y2 = data[:,7]

plt.figure()
plt.gca().set_aspect('equal', adjustable='box')

plt.plot(x1,y1,label="star 1")
plt.plot(x2,y2,":",label="star 2")

#x = np.linspace(-2e12,2e12,100)
#y = np.ones(100)*1.597399e+12
#plt.plot(x,y)

plt.legend()

plt.show()
