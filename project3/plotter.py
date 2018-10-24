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
    rx.append(x[i])
    ry.append(y[i])
    rz.append(z[i])
   
        
fig = plt.figure()
ax = fig.gca(projection = '3d')


# Plot the surface


ax.plot(rx,ry, rz)

ax.set_xlabel("rx")
ax.set_ylabel("ry")
ax.set_zlabel("rz")
plt.title("Ten years of the earths orbit with n = 1e9 - With Verlet")


plt.savefig("Verlet9.eps")
plt.show();

