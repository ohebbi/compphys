import numpy as np
import random
from math import exp, log
import matplotlib.pyplot as plt
import sys

f = open("task.txt", "r")
rx = []
ry = []

f1 = f.readlines()
for i in f1:
    a = i.split(" ")
    rx.append(float(a[0]))
    ry.append(float(a[1]))

f.close()

plt.plot(rx,ry)
plt.xlabel("rx")
plt.ylabel("ry")
plt.legend(['rx', 'ry'])
plt.title("n=1000 for the general matrix")
plt.show()
