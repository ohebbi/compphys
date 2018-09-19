import numpy as np
import math
import random

n = 10

a = -np.ones(n-1)
d = 2*np.ones(n)
d[0] += 1E-15

"""
d[0]=random.randint(0,10)
for i in range(n):
    a[i] = random.randint(0,10)
    d[i+1] = random.randint(0,10)
"""
#print d, a

A=1
eps = 1E-15
while A>eps:
    print A
    m = np.argmax(abs(a))
    #print a[m]
    theta = 0.5*np.arctan(2*a[m]/(d[m+1]-d[m]))
    s=np.sin(theta)
    c=np.cos(theta)
    x=d[m]*c**2   - 2.0*a[m]*c*s + d[m+1]*s**2
    y=d[m+1]*c**2 + 2.0*a[m]*c*s + d[m]*s**2
    a[m]=(d[m] - d[m+1])*c*s + a[m]*(c**2 - s**2)
    d[m] = x
    d[m+1] = y
    z=0.
    for i in range(len(a)):
        z+=a[i]**2
    A=np.sqrt(z)



m1 = np.reshape(np.zeros(n**2),(n,n))

for i in range(n):
    for j in range(n):
        if i==j:
            m1[i,i] = 2.
        elif i+1 == j:
            m1[i,j] = -1.
        elif i-1 == j:
            m1[i,j] = -1.

"""
A=1.
while A>eps:
    m = np.argmax(abs(m1))
    print m


    theta = (np.pi/4.0)-0.5*np.arctan((m1[m][m]-m1[m+1][m+1])/(2.0*a[m]))
    s=np.sin(theta)
    c=np.cos(theta)
    x=d[m+1]*c**2   - 2.0*a[m]*c*s + d[m]*s**2
    y=d[m]*c**2 + 2.0*a[m]*c*s + d[m+1]*s**2
    a[m]=(d[m+1] - d[m])*c*s + a[m]*(c**2 - s**2)
    d[m+1] = x
    d[m] = y

    z=0.
    for i in range(len()):
        z+=a[i]**2
    A=np.sqrt(z)
"""
for i in range(len(d)):
    print d[i], "           ", np.linalg.eigvalsh(m1)[i]
