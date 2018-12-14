
import numpy as np
import matplotlib.pyplot as plt
import sys




start = 400
slutt = 800
n = 10	
#f = open("50dif.txt", "r")

L = np.linspace(start,slutt,n)

j=0
Diffusjon = []
Temperatur = []
totalForhold = np.zeros(9999)
j=0
Diffusjon = []
Temperatur = []
tInit = []
tSlutt = []
tForhold= []
while start <= slutt:


    f = open("%iboundary_5e.txt" %start, "r")
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
        T.append(float(a[1])*119.8)
        Ek.append(float(a[2]))
        Ep.append(float(a[3]))
        Et.append(float(a[4]))
        Diff.append(float(a[5]))
        rho.append(float(a[6]))
    Diffusjon.append(np.mean(Diff[-1000:-1])*8.556e-2)
    tInit.append(T[0])
    tSlutt.append(np.mean(T[-1000:-1]))
    #Average values
    
    start += n
    tForhold.append(tSlutt[-1]/T[0])
    #for i in range(len(T)):
        
    
        #totalForhold[i]+=(tForhold[i])
    j += 1
    #if j%1==0:
        #plt.plot(t,tForhold, 'b,')
"""    
for i in range(len(totalForhold)):
    totalForhold[i]=totalForhold[i]/((slutt-400)/n)
plt.plot(t,totalForhold, 'crimson')
plt.xlabel("Time [MD-units]")
plt.ylabel("T/$T_0$")
#plt.savefig("TemperatureofTime.pdf")
plt.show()

plt.scatter(tSlutt, Diffusjon, c=tInit, s = 50)
plt.xlabel("T[K]")
plt.ylabel("D[$cm^2$ /s]")
plt.axis([0, 400, -0.00005, 0.00042])
plt.set_cmap('cool')
plt.grid()

plt.colorbar().ax.set_ylabel('$T_0$ [K]')
plt.savefig("diffconst.pdf")
plt.show()
"""
print len(tSlutt), len(tForhold)

plt.subplot(1, 2, 1)
plt.scatter(tInit, tSlutt, c=tInit, s=50)
plt.xlabel("$T_0$ [K]")
plt.ylabel("T [K]")
plt.set_cmap('cool')
plt.grid()

plt.subplot(1, 2, 2)
plt.scatter(tSlutt, tForhold, c=tInit, s=50)
plt.xlabel("T [K]")
plt.ylabel("$T/T_0$")
plt.set_cmap('cool')
plt.grid()



plt.colorbar().ax.set_ylabel("$T_0$ [K]")
plt.tight_layout()
plt.savefig('TempOnTime.pdf')
plt.show()


"""
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
