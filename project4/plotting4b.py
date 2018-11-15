import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import stats  

f = open("plot.txt", "r")
a = []
b = []
c = []
d = []
n = []

E=-7.9836/4
Cv=0.1283/4
X=15.973/4

while len(n)<1.0e6:
    i = f.readline()
   
    
    n.append(float(i))
    #a.append(float(allah[1]))
    #b.append(float(allah[2]))
    # c.append(float(allah[3]))
    #d.append(float(allah[4]))


x = np.linspace(min(n)-4,max(n)+4, 40)
y = np.exp(-((x-(-1.24))**2)/(2*8.11126))*1.0/(np.sqrt(2*np.pi*8.11126))

plt.hist(n, normed = 1, bins = 23)
plt.plot(x,y)
plt.xlabel("Energy (1/J)")
plt.ylabel("Proabability (%)")
plt.legend(["probability density","nummerical values"])
f.close()
"""
plt.subplot(1, 2, 1)
plt.plot(n[2:],b[2:])
plt.legend(["$\sigma^2$"])
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Variance in energy")



plt.subplot(1, 2, 2)
plt.plot(n[-100:],b[-100:])
plt.legend(["$\sigma^2$"])
plt.xlabel("Monte Carlo cycles")
plt.ylabel("Variance in energy")

plt.subplot(1, 2, 2)
plt.plot(n,c)
plt.legend(["<|M|>"])
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Mean magnetization")

plt.subplot(2, 2, 4)
plt.plot(n,d)
plt.legend(["$\chi$"])
plt.xlabel("Number of Monte Carlo cycles")
plt.ylabel("Susceptibility")

"""



plt.tight_layout(pad=0.4, w_pad=1.5, h_pad=1.0)

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
plt.savefig("L20prob24.pdf")
plt.show()

