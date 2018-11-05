import numpy as np
import matplotlib.pyplot as plt

e=[]
t=[]
f = open("plot1.txt", "r")
f1 = f.readlines()
del f1[-1]
for i in f1:
    a = i.split(" ")
    t.append(float(a[0]))
    e.append(float(a[2]))

f.close()


plt.plot(t,e)
plt.show()

