# -*- coding: utf-8 -*-
"""
Created on Wed Sep 05 11:16:00 2018

@author: Erlendo
"""

import numpy as np
import random
from math import exp, log
import matplotlib.pyplot as plt
import sys

f = open("task.txt", "r")
x = []
v = []
u = []
f1 = f.readlines()
for i in f1:
    a = i.split(" ")
    x.append(float(a[0]))
    v.append(float(a[1]))
    u.append(float(a[2]))
f.close()

plt.plot(x,v)
plt.plot(x, u)
plt.xlabel("x")
plt.ylabel("u")
plt.legend(['Numerical solution', 'Exact solution'])
plt.title("n=1000 for the general matrix")
plt.show()

h = []
eps = []

for i in range (1,8):
    fy = open("taskdn"+str(i)+".txt", "r")
    f1 = fy.readlines()
    for j in f1:
        a = j.split(" ")
        eps.append(np.log10(float(a[0])))
        h.append(np.log10(float(a[1])))
    fy.close()
plt.plot(h, eps)
plt.title("Logarithmic Plot The Error as a Function of h")
plt.legend(['log(eps)'])
plt.xlabel("log(h)")
plt.ylabel("log(eps)")
plt.show()
