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
n = []




f1 = f.readlines()
for i in f1:
    a = i.split(" ")
    x.append(float(a[0]))
    y.append(float(a[1]))
    z.append(float(a[2]))
    n.append(int(a[3]))

f.close()


fig = plt.figure()

ax = fig.gca(projection = '3d')

u = np.linspace(0, 2 * np.pi, 1000)
v = np.linspace(0, np.pi, 1000)
xu = 0.8*np.outer(np.cos(u), np.sin(v))
yu = 0.8*np.outer(np.sin(u), np.sin(v))
zu = 0.00008*np.outer(np.ones(np.size(u))*0.04, np.cos(v))


ax.plot_surface(0.04*xu, 0.04*yu, zu, color='r')






for j in range (0, (max(n)+1)):
    rx = []
    ry = []
    rz = []
    for i in range(0, len(x)):
    
        if(n[i] == j):
            rx.append(x[i])
            ry.append(y[i])
            rz.append(z[i])
    
    ax.plot(rx,ry, rz) 
    
     
ax.set_xlabel("rx[AU]")
ax.set_ylabel("ry[AU]")
ax.set_zlabel("rz[AU]")


plt.show();
