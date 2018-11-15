import numpy as np
import matplotlib.pyplot as plt

e=[]
t=[]
cv=[]
m=[]
x=[]
L=60

f1 = open("L%s.txt" %L, "r")
f11 = f1.readlines()

for i in f11:
    a = i.split(" ")
    t.append(float(a[0]))
    e.append(float(a[1]))
    cv.append(float(a[2]))
    m.append(float(a[3]))
    x.append(float(a[4]))

f1.close()


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

plt.savefig("4eL%i.pdf" %L)
plt.show()
