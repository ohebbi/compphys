import numpy as np
import random
from math import exp, log
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D


f = open("values3.txt", "r")
x = []
y = []
z = []
rx = []
ry = []
rz = []
jrx = []
jry = []
jrz = []


f1 = f.readlines()
for i in f1:
    a = i.split(" ")
    x.append(float(a[0]))
    y.append(float(a[1]))
    z.append(float(a[2]))

f.close()

for i in range(0, len(x)):
    if(i%2==0):
        rx.append(x[i])
        ry.append(y[i])
        rz.append(z[i])
    else:
        jrx.append(x[i])
        jry.append(y[i])
        jrz.append(z[i])
        
        
fig = plt.figure()
ax = fig.gca(projection = '3d')

u = np.linspace(0, 2 * np.pi, 1000)
v = np.linspace(0, np.pi, 1000)
x = 0.8*np.outer(np.cos(u), np.sin(v))
y = 0.8*np.outer(np.sin(u), np.sin(v))
z = 0.00008*np.outer(np.ones(np.size(u))*0.04, np.cos(v))

# Plot the surface
ax.plot_surface(0.04*x, 0.04*y, z, color='r')

ax.plot(rx,ry, rz)
ax.plot(jrx, jry, jrz)


ax.set_xlabel("rx")
ax.set_ylabel("ry")
ax.set_zlabel("rz")
plt.title("n=100000 for the general matrix")
plt.show();

#plt.savefig("Banejorda.jpg")


