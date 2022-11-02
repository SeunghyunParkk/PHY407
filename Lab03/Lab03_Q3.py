# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 20:34:04 2022

All questions are done by Seunghyun Park
Lab03_Q3
"""

# Q3 (a)
import numpy as np
import matplotlib.pyplot as plt
import math


def gaussxw(N): 
    # Initial approximation to roots of the Legendre polynomial 
    a = np.linspace(3, 4*N-1, N)/(4*N+2) 
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))  
    # Find roots using Newton's method 
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
        
    # Calculate the weights 
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp) 
    
    return x,w 


# Define Hermite polynomial
def H(n,x):
    if n<0:
        return 0
    # H_0(x) = 11 
    if n == 0:
        return 1
    # H_1(x) = 2*x
    elif n == 1:
        return 2*x
    # H_(n)(x) = 2*x*H_(n-1)(x) - 2*(n-1)*H_(n-2)(x) 
    else: 
        return 2*x*H(n-1,x) -2*(n-1)*H(n-2,x)

# Define nth energy level of the one-dimensional quantum harmonic oscillator
def psi(n,x):
    return 1/np.sqrt(float(2**n)*float(math.factorial(n))*np.sqrt(np.pi))*np.exp(-x**2/2)*H(n,x)


x = np.linspace(-4,4,300)
plt.figure()
plt.title("Wave function of the nth evergy level", fontsize = '14')
plt.plot(x,psi(0,x), label= 'n = 0')
plt.plot(x,psi(1,x), label= 'n = 1')
plt.plot(x,psi(2,x), label= 'n = 2')
plt.plot(x,psi(3,x), label= 'n = 3')
plt.legend( fontsize = '14')
plt.xlabel("x", fontsize = '14')
plt.ylabel("\u03A8(x)", fontsize = '14')
plt.grid()

# Q3 (b)
x = np.linspace(-10,10,300)
plt.figure()
plt.title("Wave function of the 30th evergy level", fontsize = '14')
plt.plot(x,psi(30,x))
plt.xlabel("x", fontsize = '14')
plt.ylabel("\u03A8(x)", fontsize = '14')
plt.grid()

# Q3 (c)

# Define derivative of wavefunction dpsi/dx
def der_psi(n,x):
    nume = np.exp(-x**2/2)*(-x*H(n,x)+2*n*H(n-1,x))
    deno = np.sqrt(float(2**n)*float(math.factorial(n))*np.sqrt(np.pi))
    return nume/deno


# Define integrand of <x^2>
# Change variable x to tangent(z)
def x_u(n,z):
    return np.tan(z)**2*abs(psi(n,np.tan(z)))**2/(np.cos(z))**2

# Define integrand of <p^2>
# Change variable x to tangent(z)
def p_u(n,z):
    return abs(der_psi(n,np.tan(z)))**2/np.cos(z)**2


N = 100 # Number of points
a = -np.pi/2 # Lower bound of integral
b = np.pi/2 # Upper bound of integral
x,w = gaussxw(N)
xp = 0.5 * (b - a) * x + 0.5 * (b + a)
wp = 0.5 * (b - a) * w
ns = np.arange(0,16,1)
po = [] # <X^2>
mo = [] # <p^2>
E = [] # total energy

# Calculate <x^2> for a given value of n using for loop
for n in ns:
    s_x = 0
    for i in range(N):
        s_x += wp[i] * x_u(n,xp[i])
    po.append(s_x)
    print('<x^2> for n =',n,' is',np.format_float_scientific(s_x, precision=2))
print()


# Calculate <p^2> for a given value of n using for loop

for n in ns:
    s_p = 0
    for i in range(N):
        s_p += wp[i] * p_u(n,xp[i])
    mo.append(s_p)
    print('<p^2> for n =',n,' is',np.format_float_scientific(s_p, precision=2))
print()

# Calculate the Energy using for loop
for i in range(0,len(po)):
    E.append(0.5*(po[i]+mo[i]))

# Print the total energy of the oscilataor using for loop
for i in range(0,len(ns)):
    print('The total energy of the oscillator when n=',ns[i],'is',np.format_float_scientific(E[i], precision=2))
print()   

u_x = np.sqrt(po) # uncertainty in position
u_m = np.sqrt(mo) # uncertainty in momentum

for i in range(0,len(ns)):
    print('The uncertainty in position when n=',ns[i],'is',np.format_float_scientific(u_x[i], precision=2))
print()       
for i in range(0,len(ns)):
    print('The uncertainty in momentum when n=',ns[i],'is',np.format_float_scientific(u_m[i], precision=2))
print()   






