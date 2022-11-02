# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 18:41:31 2022

@author: Seunghyun Park
Lab02_Q2 
"""
# Q2 (b)
import numpy as np
import time

def f(x):
    return 4/(1+x**2)


def trapezoidal(a,b,N):
    # a = lower bound of integral
    # b = upper bound of integral
    # N = number of slice
    h = (b-a)/N # Width of slice
    c = 0   # Assign initial value of f(a+k*h)
    for k in range(1,N):
        c += f(a+k*h)
    c = h*(0.5*f(a)+0.5*f(b)+c)
    return c

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

print('The exact value of the integral is', np.pi)
print("The value of the integral using Trapezoidal's rule is",trapezoidal(0,1,4))
print("The value of the integral using Simpson's rule is",simpson(0,1,4))

# Q2 (c)
# Estimate the number of slice required to obtain an error of O(10^-9) for trapezoidal
n_trape =  2**14
if np.pi - trapezoidal(0,1,n_trape) < 10**(-9):
    print('For the Trapezoidal method, the slices to approximate the integral with an error of O(10^-9)', n_trape, 'and the error is ',np.pi-trapezoidal(0, 1, n_trape))

# Estimate the number of slice required to obtain an error of O(10^-9) for simpson
n_simp =  2**5
if np.pi - simpson(0,1,n_simp) < 10**(-9):
    print('For the Simpson method, the slices of approximate the integral with an error of O(10^-9)', n_simp, 'and the error is ',np.pi-simpson(0, 1, n_simp))

# measure the time taken for calculating trapezoidal
t1 = time.time()  #Start time measure for trapezoidal method
trape = np.pi - trapezoidal(0,1,n_trape) #Find an error of O(10-9) for trapezoidal
t2 = time.time() #End time measure for trapezoidal method
t_trape = t2-t1 #Find time it took for trapezoidal method


# measure the time taken for calculating simpson
t3 = time.time() #Start time measure for simpson method
simp = np.pi - simpson(0,1,n_simp) #Find an error of O(10-9) for simpson
t4 = time.time() #End time measure for simpson method
t_simp = t4-t3 #Find time it took for simpson method

print("The times takes to compute the integral with an error of O(-9) for Trapezoidal is",t_trape)
print("The times takes to compute the integral with an error of O(-9) for Simpson is",t_simp)

# Q2 (d)
# define error estimation of trapezoidal method
def error_trapezoid(i1,i2):
    error = 1/3*(i2-i1)
    return error

i1 = trapezoidal(0,1,16)
i2 = trapezoidal(0,1,32)

print("The error estimation of Trapezoidal Method is ", error_trapezoid(i1, i2))

# Q2 (e)
# Define function of error estimation of simpson method
def error_simpson(i1,i2):
    error = 1/15*(i2-i1)
    return error

i1 = simpson(0,1,10)
i2 = simpson(0,1,2**6)
print("The error estimation of simpson's method is ",error_simpson(i1,i2))