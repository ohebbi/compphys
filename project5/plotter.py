    
import numpy as np
import matplotlib.pyplot as plt
import sys
t=[]
T = []
Ek = []
Ep = []
Et = []

f = open("statistics.txt", "r")
f.readline()
for line in f:
    a = line.split(" ")
    t.append(float(a[0])*1.00224e-13)
    T.append(float(a[1])*119.735) 
    Ek.append(float(a[2]))
    Ep.append(float(a[3]))
    Et.append(float(a[4]))
siste = 8000
T0 = T[0]

for i in range(0,len(T)):
    T[i] = T[i]/T0
 

plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.plot(t[0:siste], T[0:siste], ',')    
    
"""
plt.subplot(2,2,1)
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.plot(t[0:siste], T[0:siste], '.')

plt.subplot(2,2,2)
plt.xlabel("Time (s)")
plt.ylabel("Kinetic energy (eV)")
plt.plot(t[0:siste],Ek[0:siste])

plt.subplot(2,2,3)
plt.xlabel("Time (s)")
plt.ylabel("Potential energy (eV)")
plt.plot(t[0:siste],Ep[0:siste])

plt.subplot(2,2,4)
plt.plot(t[0:siste],Et[0:siste])
plt.xlabel("Time (s)")
plt.ylabel("Total energy (eV)")
plt.tight_layout()
"""
plt.show()
