# -*- coding: utf-8 -*-
'''
PHY407 Lab04

Done by Seunghyun Park

Q2

'''

import numpy as np
from numpy.linalg import eigvalsh
from numpy.linalg import eigh
from scipy import constants
import matplotlib.pyplot as plt

M = 9.1094e-31 # mass of electron in kg
charge = 1.6022e-19 # charge of electron in C
L = 5e-10 # width of well in m
a = 10*1.602e-19 # in J
hbar = constants.hbar
# 1ev = 1.602 × 10−19 joule.

#(c)
# Define H element value
def Hmatrix(m,n):
    # if m is not equal to n
    if (m !=n) :
        # if they both odd or even
        if ((m%2==0 and n%2 == 0) or (m%2 == 1 and n%2 ==1)):
            H_1 = 0
        # one even, oneodd
        else:
            H_1 = -8*a*m*n/(np.pi**2*(m**2-n**2)**2)
    # m = n
    else:
        H_1 = 1/2*a+np.pi**2*hbar**2*m**2/(2*M*L**2)
    return H_1


# Define matrix 
def H(mmax,nmax):
    # Initialize zero matrix
    H = np.zeros([mmax, nmax], float)
    for m in range(1, mmax+1):
        for n in range(1, nmax+1):
            H[m-1, n-1] = Hmatrix(m, n)
    return H

# Obtain eigenvalue using eigvalsh, 10*10 array
eigenvalue = eigvalsh(H(10,10))/1.602e-19
print()
print("The first ten eigenvales using 10 X 10 array:",eigenvalue,'eV')
print()




#(d)

def H(mmax,nmax):
    H = np.zeros([mmax, nmax], float)
    for m in range(1, mmax+1):
        for n in range(1, nmax+1):
            H[m-1, n-1] = Hmatrix(m, n)
    return H

# Obtain eigenvalue using eigvalsh, 100*100 array
eigenvalue1 = eigvalsh(H(100,100))/1.602e-19
print("The first ten eigenvales using 100 X 100 array",eigenvalue1[:10] ,'eV')
print()


# (e)
def gaussxw(N): 
    a = np.linspace(3, 4*N-1, N)/(4*N+2) 
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))  
    epsilon = 1e-15 
    delta = 1.0 
    while delta>epsilon: 
        p0 = np.ones(N, float) 
        p1 = np.copy(x) 
        for k in range(1, N): 
            p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1) 
        dp = (N+1)*(p0-x*p1)/(1-x*x) 
        dx = p1/dp 
        x -= dx 
        delta= max(abs(dx)) 
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp) 
    return x,w 

#define wavefunction with variable x, state, dimension
def wavefunction(x,state,mmax):
    # find eigenvalue and eigenvector using eigh function
    eva, evc = eigh(H(mmax,mmax))
    #Assign initial value of wavefunction
    psi = 0.0
    for i in range(mmax):
        psi += evc[i,state]*np.sin(np.pi*(i+1)*x/L)
    return psi

#define integrand
def wave2(x,state,mmax):
    return abs(wavefunction(x,state,mmax))**2

#Calculate the integral using gaussxw
#lower bound of integral
lower = 0
#upper bound of integral
upper = L
N = 100
x, w = gaussxw(N)
xp = 0.5 * (upper - lower) * x + 0.5 * (upper + lower)
wp = 0.5 * (upper - lower) * w
s = 0.0
for i in range(N):
    s += wp[i] * wave2(xp[i],0,3)
    
#print the value of integral    
print("A =",s)


x = np.linspace(0,L,100)


plt.figure(figsize =(10,6))
plt.plot(x,wave2(x,0,100)/s,label="Ground state")
plt.plot(x,wave2(x,1,100)/s,label="First excited state")
plt.plot(x,wave2(x,2,100)/s,label="Second excited state")
plt.legend(loc = 'upper left')
plt.xlabel("x (m)",fontsize = "14")
plt.ylabel("$|\psi(x)|^2$",fontsize = "14")
plt.title("Probability Density",fontsize = "14")


