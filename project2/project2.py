import numpy as np
import math
import random

n = 10
h = 1./n #
A=1.

"""
Making a matrix and finding eigenvalues with lib func.
"""
m1 = np.reshape(np.zeros(n**2),(n,n))

for i in range(n):
    for j in range(n):
        if i==j:
            m1[i,i] = 2.
        elif i+1 == j:
            m1[i,j] = -1.
        elif i-1 == j:
            m1[i,j] = -1.
m1/=(h**2)
m2 = m1.copy() #controlling eigenvalues at the end

"""
Finding indices with max value
"""
def maksoffdiag(matrise, n, max=0,q=0,p=0):
    for i in range(n):
        for j in range(n):
            if i != j:
                aij = np.abs(matrise[i][j])
                if aij > max:
                    max = aij
                    p = i
                    q = j
    return p, q

"""
Jacobi's method
"""
eps = 1E-4
while A>eps:

    l, k = maksoffdiag(m1,n)
    tau = (m1[l][l]-m1[k][k])/(2.*m1[k][l])

    if tau >= 0: #matinf1100 curriculum before mid-term
        t = ( -tau -np.sqrt(1+tau*tau))
    else:
        t = (-tau +np.sqrt(1+tau*tau))

    c = 1./(np.sqrt(1+t*t))
    s = c*t
    for i in range(n): #
        if i != k and i != l:
                ik = m1[i][k]
                il = m1[i][l]
                m1[i][k] = ik*c - il*s
                m1[i][l] = il*c + ik*s
                m1[k][i] = m1[i][k] #symmetric matrix
                m1[l][i] = m1[i][l]

    x=m1[k][k]*c**2 - 2.0*m1[k][l]*c*s + m1[l][l]*s**2
    y=m1[l][l]*c**2 + 2.0*m1[k][l]*c*s + m1[k][k]*s**2

    m1[k][l]= 0. #(m1[k][k] - m1[l][l])*c*s + m1[k][l]*(c**2 - s**2)
    m1[l][k] = m1[k][l] #symmetric matrix

    m1[k][k] = x
    m1[l][l] = y

    z=0.

    #Frobenius norm
    for u in range(n):
        for v in range(n):
            if u != v:
                z+=(m1[u][v])**2
    A=np.sqrt(z)

"""
printing eigenvalues.
"""

for i in range(n):
    print "computed values = ",m1[i][i]
x = np.linalg.eigvalsh(m2)
for i in range(len(x)):
    print "exact values = ", x[i]
