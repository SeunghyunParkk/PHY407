# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 20:39:43 2022

@author: Juann Jeon & Seunghyun Park
"""

from scipy import constants
import matplotlib.pyplot as plt
import numpy as np

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

m = 1 #Mass, in kg
k = 12 #Spring constant, in N/m
c = constants.speed_of_light #Speed of light, in m/s

# Define the integrand
# For convenience define g(x) as 4*(original g(x))^-1
def g(x, x_0):
    fx = (x_0**2) - (x**2)
    f1 = k * fx
    f2 = (2 * m * c**2) + k * fx/2
    f3 = 2 * ((m * c**2) + k * fx/2)**2
    return 4*(c * ((f1 * f2)/f3)**0.5)**-1

#Q2 a) Done by Juann Jeon
a = 0 #Lower boundary for integral, always 0 in equation (7)
b = 0.01 #b = x_0, in m
x_0 = b 
N_1 = 8
N_2 = 16

x_1, w_1 = gaussxw(N_1)
xp_1 = 0.5 * (b - a) * x_1 + 0.5 * (b + a)
wp_1 = 0.5 * (b - a) * w_1
s_1 = 0.0
for i in range(N_1):
    s_1 += wp_1[i] * g(xp_1[i],x_0)    
print("Period Using Gaussian quadrature with N = 8:", np.format_float_scientific(s_1, precision=2)) 

x_2, w_2 = gaussxw(N_2)
xp_2 = 0.5 * (b - a) * x_2 + 0.5 * (b + a)
wp_2 = 0.5 * (b - a) * w_2
s_2 = 0.0
for i in range(N_2):
    s_2 += wp_2[i] * g(xp_2[i],x_0)
print("Period Using Gaussian quadrature with N = 16:", np.format_float_scientific(s_2, precision=2)) 


classical = 2 * constants.pi * (m/k)**0.5

print("Classical value of Period:", np.format_float_scientific(classical, precision=2))
print("Fractional error for N = 8:", np.format_float_scientific(abs(s_1 - classical)/classical, precision=2))
print("Fractional error for N = 16:", np.format_float_scientific(abs(s_2 - classical)/classical, precision=2))
print()

#Q2 b) Done by Juann Jeon
plt.figure()
plt.plot(xp_1, g(xp_1,x_0), label = 'N = 8')
plt.plot(xp_2, g(xp_2,x_0), label = 'N = 16')
plt.title('Integrands using Gaussian quadrature', fontsize = '14')
plt.xlabel('x (m)', fontsize = '14')
plt.ylabel('Integrands (s/m)', fontsize = '14')
plt.legend()
plt.grid()
plt.figure()
plt.plot(xp_1, wp_1 * g(xp_1,x_0), label = 'N = 8')
plt.plot(xp_2, wp_2 * g(xp_2,x_0), label = 'N = 16')
plt.title('Weighted value', fontsize = '14')
plt.xlabel('x (m)', fontsize = '14')
plt.ylabel('Weighted value', fontsize = '14')
plt.legend()
plt.grid()



#Q2 c) Done by Juann Jeon
x_c = (( c**2)/k)**0.5
print("For classical, to have speed of light, x_0 needs to be (in m)", np.format_float_scientific(x_c, precision=2))
print()

#Q2 d) Done by Seunghyun Park
N_200 = 200
x_200, w_200 = gaussxw(N_200)
xp_200 = 0.5 * (b - a) * x_200 + 0.5 * (b + a)
wp_200 = 0.5 * (b - a) * w_200
s_200 = 0.0
for i in range(N_200):
    s_200 += wp_200[i] * g(xp_200[i],x_0)

N_400 = 400
x_400, w_400 = gaussxw(N_400)
xp_400 = 0.5 * (b - a) * x_400 + 0.5 * (b + a)
wp_400 = 0.5 * (b - a) * w_400
s_400 = 0.0
for i in range(N_400):
    s_400 += wp_400[i] * g(xp_400[i],x_0)
    
print("Percentage error for N = 200:", (s_400-s_200)*100,'%')
print()



#Q2 e) Done by Seunghyun Park
N_4 = 200
a = 0
up = 10*x_c
x_0 = np.linspace(1,10*x_c)
T = []

for i in range(0,len(x_0)):
    b = x_0[i]
    x_41, w_41 = gaussxw(N_4)
    xp_41 = 0.5 * (b - a) * x_41 + 0.5 * (b + a)
    wp_41 = 0.5 * (b - a) * w_41
    s_41 = 0
    for i in range(N_4):
        s_41 += (wp_41[i] * g(xp_41[i],b))
    T.append(s_41)

plt.figure()
plt.plot(x_0,T, label = 'T')


# T approaches 2*pi*np.sqrt(m/k) in small amplitude
low = 2 * np.pi*np.sqrt(m/k)
plt.scatter(10*x_c,10*x_c*4/c, label='4x_0/c', c='red')
plt.scatter(1,low, label = '2$\pi$(m/k)^0.5' , c ='black')
plt.xlabel('x_0 (m)', fontsize = '14')
plt.ylabel('T (s)', fontsize = '14')
plt.title('T as a function of x_0', fontsize = '14')
plt.grid()
plt.legend()

    
    