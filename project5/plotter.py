
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

for i in L:

    f = open("%idiffusion.txt" %i, "r")
    f.readline()
    t=[]
    T = []
    Ek = []
    Ep = []
    Et = []
    Diff = []

    differseoi = 0.0
    differseoi1 = 0.0
    T0=0.0
    T1=0.0

    for line in f:
        a = line.split(" ")
        t.append(float(a[0]))
        T.append(float(a[1])*119.735)
        Ek.append(float(a[2]))
        Ep.append(float(a[3]))
        Et.append(float(a[4]))
        Diff.append(float(a[5]))

    #Average values
    for i in range(len(T)):
        T0 += T[i]
        differseoi += Diff[i]
    T1 = T0 / len(T)
    differseoi1 = float(differseoi / len(T))


    #Pretty plotting
    if j>-1:
        plt.plot(t,np.array(Diff),label='$T=%i K$'%(T1))
    else:
        plt.plot(t,np.array(Diff),',',alpha=0.02)

    Temperatur.append(T1)
    Diffusjon.append(differseoi1)

    j +=1


for i in range(len(Temperatur)):
    print(Temperatur[i], " ", Diffusjon[i])



plt.legend(loc='upper center')
plt.xlabel("Time [MD-units]")
plt.ylabel("Diffusion Coeffisient [MD-units]")
plt.axis([-2,100, 0,0.035])
plt.tight_layout()
#plt.savefig("Time_vs_Diff.pdf")
plt.show()

plt.plot(Temperatur,Diffusjon)
plt.xlabel("Temperatur [K]")
plt.ylabel("Diffusjonskoeffisient [MD-units]")
#plt.savefig("Temp_vs_Diff.pdf")
plt.tight_layout()
plt.show()

plt.subplot(2,1,1)
plt.plot(Temperatur,L)
plt.xlabel("$T$ [K]")
plt.ylabel("$T_0$ [K]")

L = np.linspace(start, slutt, len(Diffusjon))
plt.subplot(2,1,2)
plt.plot(Temperatur,Temperatur/L)
plt.xlabel("$T$ [K] ")
plt.ylabel("$T/T_0$")
plt.tight_layout()
#plt.savefig("rateTemp.pdf")
plt.show()


"""
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

plt.subplot(2,1,1)
plt.plot(Temperatur,L)
plt.xlabel("$T$ [K]")
plt.ylabel("$T_0$ [K]")

L = np.linspace(start, slutt, len(Diffusjon))
plt.subplot(2,1,2)
plt.plot(Temperatur,Temperatur/L)
plt.xlabel("$T$ [K] ")
plt.ylabel("$T/T_0$")
plt.tight_layout()
#plt.savefig("rateTemp.pdf")
plt.show()

"""
"""
plt.subplot(2,2,1)
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.plot(t[0:siste], T[0:siste], '.')

plt.subplot(2,2,2)
plt.xlabel("Time (s)")
plt.ylabel("Kinetic energy (eV)")
plt.plot(t[0:siste],Ek[0:siste])

plt.subplot(2,2,3)
plt.xlabel("Time (s)")
plt.ylabel("Potential energy (eV)")
plt.plot(t[0:siste],Ep[0:siste])

plt.subplot(2,2,4)
plt.plot(t[0:siste],Et[0:siste])
plt.xlabel("Time (s)")
plt.ylabel("Total energy (eV)")
plt.tight_layout()
"""
