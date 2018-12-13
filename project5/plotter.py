
import numpy as np
import matplotlib.pyplot as plt
import sys



start = 500
slutt = 700
n = 50
#f = open("50diffusion.txt", "r")
L = np.linspace(start,slutt,n)

j=0
Diffusjon = []
Temperatur = []
tInit = []
tSlutt = []
while start <= slutt:
    
	
    f = open("%idiffusion.txt" %start, "r")
    f.readline()
    t=[]
    
    T = []
    Ek = []
    Ep = []
    Et = []
    Diff = []
    rho = []

    differseoi = 0.0
    differseoi1 = 0.0
    T0=0.0
    T1=0.0

    for line in f:
        a = line.split(" ")
        t.append(float(a[0]))
        T.append(float(a[1]))
        Ek.append(float(a[2]))
        Ep.append(float(a[3]))
        Et.append(float(a[4]))
        Diff.append(float(a[5]))
        rho.append(float(a[6]))
    Diffusjon.append(np.mean(Diff[-1000:-1]))    
    tInit.append(T[0])
    tSlutt.append(np.mean(T[-1000:-1]))
    #Average values
    start += n
tForhold = []
for i in range(len(tInit)):
    tForhold.append(tSlutt[i]/tInit[i])
    

plt.plot(tSlutt, tForhold, 'x-')
plt.show()

plt.plot(tSlutt, Diffusjon)
plt.show()
