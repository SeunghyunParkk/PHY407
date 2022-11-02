# -*- coding: utf-8 -*-
'''
PHY407 Lab04

Done by Seunghyun Park

Q3

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
#(a)
#define function x = 1-np.exp(-c*x)
def f(x):
    return 1-np.exp(-c*x)

accuracy = 1e-6
c = 2
#calculate x using for loop
for k in range(50):
    x1 = 1.0
    error = 1
    while error > accuracy:
        x1,x2 = f(x1),x1
        error = abs((x1-x2)/(1-(c*np.exp(-c*x2))**-1))    
    print("The solution to an accuracy of at least 1e-6 is",np.format_float_scientific(x1, precision=6))
    break


accuracy = 1e-6
cs = np.arange(0.00001,3,0.01)
xs =[]
#calculate x as function of c
for c in cs:
    x1 = 1.0
    error = 1.0
    while error > accuracy:
        x1,x2 = f(x1),x1
        error = abs((x1-x2)/(1-(c*np.exp(-c*x2))**-1))    
    xs.append(x1)
    

plt.figure()
plt.plot(cs,xs)
plt.title('Solution of x = 1-exp(-cx)', fontsize = '14')
plt.xlabel("c", fontsize = '14')
plt.ylabel("x", fontsize = '14')
plt.ylim(-0.05,1)


#(b)

accuracy = 1e-6
x1 = 1
it = 0
c = 2
#measure the iterations of relaxation method
for k in range(50):
    x1 = 1.0
    error = 1
    while error > accuracy:
        x1,x2 = f(x1),x1
        error = abs((x1-x2)/(1-(c*np.exp(-c*x2))**-1))    
        it +=1
    break
print()
print("Number of iteration takes for relaxation method:",it)
print()

ws = np.arange(0.5,1,0.01)
accuracy = 1e-6
its = []
#measure the iteration of overrelaxation method
for w in ws:
    it = 0
    c = 2
    x1 = 1.0
    error = 1.0
    while error > accuracy:
        x1,x2 = f(x1),x1
        error = abs((x1-x2)/(1-1/((1+w)*c*np.exp(-c*x2)-w)))
        it +=1
    its.append(it)
print("The number of iterations of overrelaxation method",its)
print()
#find the minimum iteration
min_it = np.where(its ==np.amin(its))
min_w = ws[min_it]
print("The minimum number of iteration of overrelaxation method:",min(its),"with w",min_w)
print()

# (c)
h = constants.h
c = constants.c
kb = constants.Boltzmann
#define fuction 5 * np.exp(-x) + x -5 = 0
def f(x):
    return 5 * np.exp(-x) + x -5

#define binary search method
def binary(x1,x2):    
    accuracy = 1e-6
    error = 1
    s_x1 = np.sign(f(x1))
    s_x2 = np.sign(f(x2))
    # if the they have opposite sign
    if s_x1 != s_x2:
        while error > accuracy:
            x_mid = 0.5*(x1+x2)
            f_mid = f(x_mid)
            s_mid = np.sign(f_mid)
            if s_mid == s_x1:
                x1 = x_mid
            else:
                x2 = x_mid
            error = abs(x1-x2)
    # if they have same sign print error
    else:
        print("Sign of f(x1) and f(x2) are same")
    return x_mid

print("The solution of the equation is ",np.format_float_scientific(binary(2,10), precision=6))
print()
#calculate wien displacement constant
b = h*c/kb/binary(2,10)
print("Wien displacement constant is", np.format_float_scientific(b, precision=6),"meter Kelvin")
print()

wave = 502e-9
#calculate the temperature of the Sun
T = b/wave
print("The temperature of the Sun is",np.format_float_scientific(T, precision=6),"K")
