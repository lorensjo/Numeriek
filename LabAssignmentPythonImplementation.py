# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:46:22 2023

@author: diete
"""

#######exercise (c) 

import numpy as np
import matplotlib.pyplot as plt



dr = 5
r0 = 1
R = 1001
n = int((R-r0)/dr-1)

def re(j,dr):
    return j*dr+1

A0 = np.diag(np.full(n,-2,dtype=float))
A0[0,1]=2*re(0.5,dr)
A0[0,0]=-2*re(0.5,dr)

d=1
e=1
for i in range(n-2):
    k1 = re(d-0.5,dr)/re(d,dr)
    k2 = re(d+0.5,dr)/re(d,dr)
    A0[d,e-1]=k1
    
    A0[d,e+1]=k2
    
    d=d+1
    e=e+1



A = (1/dr**2)*A0 
b = np.zeros(shape=(n,1))
b[0]=(2/dr)
invA = np.linalg.inv(A)
w = invA.dot(-1*b)

ls1 = []
for i in range(n):
    ls1.append(1+dr*i)
    
ls2 = []
for i in range(n):
    ls2.append(w[i])
    
def ph(r,R):
    return -np.log(r) + np.log(R)

ls3 = []
for i in range(n):
    ls3.append(ph(ls1[i],R))
    
plt.plot(ls1,ls2, label = 'Numerical')  
plt.plot(ls1,ls3,'r--',label = 'Analitical')  
plt.legend()
plt.show()


#######exercise (g) 
dt = 10**3


def back_euler(w):
    I = np.identity(199)
    inv = np.linalg.inv(I-dt*A)
    d = w+dt*b
    return inv.dot(d)


w0 = np.zeros((199,1))    
ls = [w0]
for i in range(1,1000):
    w0 = back_euler(w0)
    ls.append(w0)

def function(t,r):
    return ls[t][r][0]*(1+dr*r)



def trap_rule(t):
    a = 0
    b = 198
    h = dr
    s = 0.5*(function(t,a) + function(t,b))
    for i in range(1,n):
        s = s + function(t,a+i)
    r = h*s
    return (2/(1001**2-1))*r


    
    

t0 = 0
lss1 = []
lss2 = []
for i in range(1000):
    lss1.append(t0+i*dt)
    lss2.append(trap_rule(i))
plt.plot(lss1,lss2)
plt.show()


#######exercise(h)

def trap_method(w):
    I = np.identity(199)
    inv = np.linalg.inv(I-(dt*A)/2)
    C = I+(dt*A)/2
    d = C.dot(w)+dt*b
    return inv.dot(d)
    

w0 = np.zeros((199,1))    
ls = [w0]
for i in range(1,1000):
    w0 = trap_method(w0)
    ls.append(w0)

def function2(t,r):
    return lss[t][r][0]*(1+dr*r)


def trap_rule2(t):
    a = 0
    b = 198
    h = dr
    s = 0.5*(function2(t,a) + function2(t,b))
    for i in range(1,n):
        s = s + function2(t,a+i)
    r = h*s
    return (2/(1001**2-1))*r

w0 = np.zeros((199,1))    
lss = [w0]
for i in range(1,1000):
    w0 = trap_method(w0)
    lss.append(w0)


t0 = 0
lsss2 = []
for i in range(1000):
    lsss2.append(trap_rule2(i))
plt.plot(lss1,lsss2)
plt.show()

