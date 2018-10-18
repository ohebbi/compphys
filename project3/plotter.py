import numpy as np
import random
from math import exp, log
import matplotlib.pyplot as plt
import sys

f = open("values3.txt", "r")
x = []
y = []
rx = []
ry = []
jrx = []
jry = []

f1 = f.readlines()
for i in f1:
    a = i.split(" ")
    x.append(float(a[0]))
    y.append(float(a[1]))

f.close()



for i in range(0, len(x)):
    if(i%2==0):
        rx.append(x[i])
        ry.append(y[i])
    else:
        jrx.append(x[i])
        jry.append(y[i])
        
    
plt.plot(rx,ry)
plt.plot(jrx,jry)


plt.xlabel("rx")
plt.ylabel("ry")
plt.legend(['rx', 'ry'])
plt.title("Our solar system. What a beauty. :)")
plt.show()
