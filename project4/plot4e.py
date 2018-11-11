import numpy as np
import matplotlib.pyplot as plt

e=[]
t=[]
cv=[]
m=[]
x=[]

f1 = open("L40-1.txt", "r")
f11 = f1.readlines()
f2 = open("L40-2.txt", "r")
f22 = f2.readlines()

for i in f11:
    a = i.split(" ")
    t.append(float(a[0]))
    e.append(float(a[2]))
    cv.append(float(a[3]))
    m.append(float(a[4]))
    x.append(float(a[5]))
for i in f22:
    a = i.split(" ")
    t.append(float(a[0]))
    e.append(float(a[2]))
    cv.append(float(a[3]))
    m.append(float(a[4]))
    x.append(float(a[5]))

f1.close()
f2.close()

plt.subplot(2, 2, 1)
plt.plot(t,e)
plt.legend(["$<E>$"])
plt.xlabel("Temperature")
plt.ylabel("Energy")

plt.subplot(2, 2, 2)
plt.plot(t,cv)
plt.legend(["$C_V$"])
plt.xlabel("Temperature")
plt.ylabel("Heat capacity")

plt.subplot(2, 2, 3)
plt.plot(t,m)
plt.legend(["<M>"])
plt.xlabel("Temperature")
plt.ylabel("Magnetization")

plt.subplot(2, 2, 4)
plt.plot(t,x)
plt.legend(["$\chi$"])
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")

plt.savefig("4eL40.pdf")
plt.show()

