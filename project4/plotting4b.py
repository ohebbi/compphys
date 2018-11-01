import numpy as np
import matplotlib.pyplot as plt



f = open("plot.txt", "r")
a = []
b = []
c = []
d = []
n = []




f1 = f.readlines()
for i in f1:
    allah = i.split(" ")
    n.append(float(allah[0]))
    a.append(float(allah[1]))
    b.append(float(allah[2]))
    c.append(float(allah[3]))
    d.append(float(allah[4]))

f.close()


plt.xlabel("Number of Monte-Carlo cycles")
plt.ylabel("Energy")
plt.title("4b")
#plt.plot(x,y)
plt.plot(n,a)
plt.plot(n,b)
plt.plot(n,c)
plt.plot(n,d)
plt.plot()
plt.legend(["Mean energy", "Heat capacity", "Magnetization", "Susceptibility"])

#plt.savefig("FixedSolarSystem.pdf")
plt.show()
