# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:12:31 2022

@author: Seunghyun Park
Lab02_Q3 
"""
import numpy as np
from scipy import constants

# Q3 (b)
h = constants.Planck # planck constant
c = constants.speed_of_light # speed of light
k = constants.k # boltzmann constnat

# Define the function inside of the integral
def f(x):
    return (x**3) / (np.exp(x)-1.)

def simpson(a,b,N):
    # a = lower bound of integral
    # b = upper bound of integral
    # N = number of slice
    h = (b-a)/N # Width of slice
    c = 0 # Assign initial value of f(a+k*h) odd
    d = 0 # Assign initial value of f(a+k*h) even
    for k in range(1,N,2):
        c += f(a+k*h)
    for k in range(2,N,2):
        d += f(a+k*h)
    tot = 1/3*h*(f(a)+f(b)+4*c+2*d) # Calculate the integral
    return tot

N = 10000
# since f(x=0) = undefined, instead of using a = 0, use very small 0.00000001
a = 0.00000001
# since upper bound cannot be infinite, assign b =500 so the upper bound is larger enough
b = 500

# Temperature
T = 273.15

# Define W(T) as function of Temperature
def W(T):
    W = 2*np.pi*k**4/h**3/c**2*simpson(a,b,N)*T**4
    return W

print("The calculated W(T) when T = 273.15K is",W(T),'Wm^-2')
print("The actual W(T) when T = 273.15 K is ",constants.sigma*T**4,'Wm^-2')

# Calculating accuracy of the output W(T)
accuracy = ((constants.sigma*T**4)/W(T))*100
print("The accuracy of the calculated value is", accuracy,"%")

# Q3 (c)

# Calculate the Stefan-Boltzmann Constant
s_b = W(T)/T**4

print('The calculated Stefan-Boltzmann constant is',s_b,'Wm^-2K-^4')
print('The Stefan-Boltmann constant is',constants.sigma,'Wm^-2K-^4')
