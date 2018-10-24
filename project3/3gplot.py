# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:35:05 2018

@author: Erlendo
"""

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
    x.append(float(a[1]))
    y.append(float(a[0]))
        

f.close()

plt.scatter(x,y)

