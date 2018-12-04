
import numpy as np
import matplotlib.pyplot as plt
import sys

T = []
Ek = []
Ep = []
Et = []

f = open("statistics.txt", "r")
f.readline()
for line in f:
    a = line.split(" ")
    T.append(float(a[0])*119.735) 
    Ek.append(float(a[1]))
    Ep.append(float(a[2]))
    Et.append(float(a[3]))
x = np.linspace(0, T[-1], len(T))*1e-13

plt.subplot(2,2,1)
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.plot(x, T)

plt.subplot(2,2,2)
plt.xlabel("Time (s)")
plt.ylabel("Kinetic energy (eV)")
plt.plot(x,Ek)

plt.subplot(2,2,3)
plt.xlabel("Time (s)")
plt.ylabel("Potential energy (eV)")
plt.plot(x,Ep)

plt.subplot(2,2,4)
plt.plot(x,Et)
plt.xlabel("Time (s)")
plt.ylabel("Total energy (eV)")
plt.tight_layout()
plt.show()
