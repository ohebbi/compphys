import numpy as np
import matplotlib.pyplot as plt

f = open("plot.txt", "r")
a = []
b = []
c = []
d = []
n = []

E=-7.9836/4
Cv=0.1283/4
X=15.973/4

while len(a)<0.8e5:
    i = f.readline()
    allah = i.split(" ")
    n.append(float(allah[0]))
    a.append(float(allah[1]))
    b.append(float(allah[2]))
    c.append(float(allah[3]))
    d.append(float(allah[4]))


f.close()

plt.subplot(2, 2, 1)
plt.plot(n[:1000],a[:1000])
plt.legend(["<E>"])
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean energy")


plt.subplot(2, 2, 2)
plt.plot(n[:1000],b[:1000])
#plt.legend(["$C_V$"])
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Heat capacity")

plt.subplot(2, 2, 3)
plt.plot(n[:1000],c[:1000])
plt.legend(["<M>"])
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean magnetization")

plt.subplot(2, 2, 4)
plt.plot(n[:1000],d[:1000])
plt.legend(["$\chi$"])
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Susceptibility")


plt.tight_layout()

"""
plt.xlabel("Number of Monte-Carlo cycles")
plt.ylabel("")
plt.title("4b")
#plt.plot(x,y)
plt.plot(n[-50000:-1],a[-50000:-1])
#plt.plot(n,b)
#plt.plot(n,c)
#plt.plot(n,d)
plt.plot()
plt.legend(["Mean energy", "Heat capacity", "Magnetization", "Susceptibility"])

print b[-1]
plt.hist(a, normed=True, bins=65)
plt.xlabel("Energy (J)")
plt.ylabel("Normalized probability")
plt.savefig("paral.pdf")
"""
#plt.savefig("allup.pdf")
plt.show()
