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


plt.savefig("Relativity.eps")
plt.show();
