# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 22:17:06 2022

All questions are done by Seunghyun Park
Lab03_Q1
"""
import numpy as np
import matplotlib.pyplot as plt


C = 1e-16

# Q1 (a)

# define function f(x)
def f(x):
    return np.exp(-x**2)

# Analytic Value of the derivative
def true_value(x):
    return -2*x*np.exp(-x**2)

# Define forward difference
def forward(x,h):
    # c is the value of derivative
    # dx/dy = (f(x+h) - f(x))/h
    c = (f(x+h) - f(x))/h
    # return c
    return c

# Calculate the derivative using for loop
for i in range(0,17):
    # h becomes 10 times bigger in each step
    der = forward(0.5, 10**(-16+i))
    print("The derivative of function f(x) using forward difference when h =",10**(-16+i),"is",np.format_float_scientific(der, precision=2))
print()


print('The analytic value of the derivative is',true_value(0.5))
print()
# Calculate the numerical error of the derivative with each value of h using for loop
error = []
for i in range(0,17):
    c = abs(true_value(0.5)- forward(0.5, 10**(-16+i)))
    error.append(c)
    print("The absolute value of error using forward difference when h =",10**(-16+i),"is",np.format_float_scientific(c, precision=2))
print()

# Q1 (b)

h = [10**-16,10**-15,10**-14,10**-13,10**-12,10**-11,10**-10,10**-9,10**-8,10**-7,10**-6,10**-5,10**-4,10**-3,10**-2,10**-1,1]
plt.figure()
plt.loglog(h,error)
plt.title("log scale Error of Forward Difference", fontsize = '14')
plt.xlabel("log(h)", fontsize = '14')
plt.ylabel("log scale Error", fontsize = '14')
plt.grid()


# find the index of minimum point on the graph
log_error = np.log(error)
minimum = np.where(log_error == np.amin(log_error))

print("The roundoff error dominates until", h[8], 'and the truancated error dominates after',h[8] )
print()

# Q1 (c)

# define central diffrence
def central(x,h):
    # c is the value of derivative
    # dx/dy = (f(x+h) - f(x))/h
    c = (f(x+h/2) - f(x-h/2))/h
    # return c
    return c

# Calculate the derivative using for loop
for i in range(0,17):
    der = central(0.5, 10**(-16+i))
    print("The derivative of function f(x) using central difference when h =",10**(-16+i),"is",np.format_float_scientific(der, precision=2))
print()
error_central = []
# Calculate the numerical error of central difference using for loop
for i in range(0,17):
    c = abs(true_value(0.5)- central(0.5, 10**(-16+i)))
    error_central.append(c)
    print("The absolute value of error using central difference when h =",10**(-16+i),"is",np.format_float_scientific(c, precision=2))

h = [10**-16,10**-15,10**-14,10**-13,10**-12,10**-11,10**-10,10**-9,10**-8,10**-7,10**-6,10**-5,10**-4,10**-3,10**-2,10**-1,1]

plt.figure()
plt.loglog(h,error, label = 'Forward difference')
plt.loglog(h,error_central, label = 'Central difference')
plt.title("log scale Error", fontsize = '14')
plt.xlabel("log(h)", fontsize = '14')
plt.ylabel("log scale Error", fontsize = '14')
plt.legend(loc = 'lower left',fontsize = '12')
plt.grid()















