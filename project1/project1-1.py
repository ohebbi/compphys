import numpy as np
import random
from math import exp, log
import matplotlib.pyplot as plt
import sys

n=int(10**float(sys.argv[1]))
h=1.0/(n+1)
x=np.linspace(0,1,n)


a=np.zeros(n-1)
b=np.zeros(n)
c=np.zeros(n-1)

v=np.zeros(n)
b_tilde=np.zeros(n) #solution

def f(x):
    return 100*exp(-10*x)


for i in range(n):
    b_tilde[i]=h**2*f(x[i])


#b)
def task_b():
    b[0]=random.uniform(0.0, 100.0)
    for i in range(n-1):
        a[i]=random.uniform(0.0, 100.0)
        b[i+1]=random.uniform(0.0, 100.0)
        c[i]=random.uniform(0.0, 100.0)
    return a, b, c


#c)
def task_c():
    b[0]=2
    for i in range(n-1):
        a[i]=-1
        b[i+1]=2
        c[i]=-1
    return a, b, c

task_c()

#forward substitution
for i in range(1,n):
    s=a[i-1]/b[i-1]
    b[i]=b[i]-s*c[i-1]
    b_tilde[i]=b_tilde[i]-s*b_tilde[i-1]


#backward substitution
v[n-1]=b_tilde[n-1]/b[n-1]
for i in range(1,n):
    j=n-1-i
    v[j]=(b_tilde[j]-c[j]*v[j+1])/b[j]

#exact solution
u=np.zeros(n)
for i in range(n):
    u[i]=(1-(1-exp(-10))*x[i]-exp(-10*x[i]))

def rel_error(v, u):
    return log(abs((v-u)/u))

eps=np.zeros(n)
for i in range(n):
    eps[i]=rel_error(v[i], u[i])


    
print u
"""

plt.plot(x, v)
plt.plot(x, u)
plt.show()
"""






"""
#need to read a, b and c values from a file (and b_tilde?)
file=open('values.txt', 'r')
values=file.readlines()
file.close()

for i in range(len(values)):
    values[i]=float(values[i][:-1])

n=int((len(values)+2)/4)


#inserts values form file to a, b, c and b_tilde
b[0]=values[n-1]
b_tilde[0]=values[3*n-1]
for i in range(n-1):
    a[i]=values[i]
    b[i+1]=values[n+i]
    c[i]=values[2*n-1+i]
    b_tilde[i+1]=values[3*n-1+i]
"""
