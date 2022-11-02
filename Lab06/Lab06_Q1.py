# -*- coding: utf-8 -*-
"""
Lab06_Q1
Created by Juann Jeon
Trajectories of particle
using Verlet method
"""

import numpy as np
import matplotlib.pyplot as plt

#(b)

#Define a function that calculates x,y components of acceleration
def f(r, t): 
    # x component
    x = r[0] - r[2]
    # y component
    y = r[1] - r[3]
    # calculate x component acceleration of particle 1
    fx1 = -((24 * x) / (x**2 + y**2)**4) + ((48 * x) / (x**2 + y**2)**7)
    # calculate y component acceleration of particle 1
    fy1 = -((24 * y) / (x**2 + y**2)**4) + ((48 * y) / (x**2 + y**2)**7)
    # calculate x component acceleration of particle 2
    fx2 = -fx1
    # calculate y component acceleration of particle 2
    fy2 = -fy1
    #return acceleration
    return np.array([fx1, fy1, fx2, fy2], float) 

# time start
a = 0.0 
# time end
b = 1.0 
# number of time step
N = 100 
# step size
h = (b-a)/N 
#assign time
t_p = np.arange(a, b, h)

# assign empty list for position for the particles
x_p11 = []
y_p11 = []
x_p12 = []
y_p12 = []
#assign initial velocity
v1 = [0, 0, 0, 0]

# r1=[4, 4], r2=[5.2, 4]
r1 = np.array([4.0, 4.0, 5.2, 4.0], float)

# assign position of particles
x_p11.append(r1[0]) 
y_p11.append(r1[1])
x_p12.append(r1[2]) 
y_p12.append(r1[3])

# calculate v(t+h/2)
v1 += h/2 * f(r1, t_p[0])
# apply verlet algorithm
for t in t_p[1:]:  
    r1 += h * v1
    k = h * f(r1, t)
    v1 += k
    
    x_p11.append(r1[0]) 
    y_p11.append(r1[1])
    x_p12.append(r1[2]) 
    y_p12.append(r1[3])   

#plot trajectories of particles
plt.figure()
# plt.plot(t_p, x_p11, '.',label="r1") 
# plt.plot(t_p, y_p11, '.') 
# plt.plot(t_p, x_p12, '.',label ="r2") 
# plt.plot(t_p, y_p12, '.') 
# plt.title("Tragectories of particles, r1=[4,4]&r2=[5.2,4]",fontsize= "13")
# plt.xlabel("time (s)", fontsize ="13")
# plt.ylabel("x (arbitrary unit)", fontsize ="13")
plt.plot(x_p11, y_p11, '.', label ='r1')
plt.plot(x_p12, y_p12, '.', label = 'r2')     
plt.title("Tragectories of particles, r1=[4,4]&r2=[5.2,4]",fontsize= "13")
plt.xlabel("x (arbitrary unit)",fontsize= "13")
plt.ylabel("y (arbitrary unit)",fontsize= "13")
plt.legend()
plt.figure()

# assign empty list for position for the particles
x_p21 = []
y_p21 = []
x_p22 = []
y_p22 = []
#assign initial velocity
v2 = [0, 0, 0, 0]

# r1=[4.5, 4], r2=[5.2, 4]
r2 = np.array([4.5, 4.0, 5.2, 4.0], float)

# assign position of particles
x_p21.append(r2[0]) 
y_p21.append(r2[1])
x_p22.append(r2[2]) 
y_p22.append(r2[3])

# calculate v(t+h/2)
v2 += h/2 * f(r2, t_p[0])
# apply verlet algorithm
for t in t_p[1:]:  
    r2 += h * v2
    k = h * f(r2, t)
    v2 += k
    
    x_p21.append(r2[0]) 
    y_p21.append(r2[1])
    x_p22.append(r2[2]) 
    y_p22.append(r2[3])   
    
#plot trajectories of particles
plt.figure()
#plt.plot(t_p, x_p21, '.') 
#plt.plot(t_p, y_p21, '.') 
#plt.plot(t_p, x_p22, '.') 
#plt.plot(t_p, y_p22, '.') 
plt.plot(x_p21, y_p21, '.', label='r1')
plt.plot(x_p22, y_p22, '.', label ='r2')
plt.title("Tragectories of particles, r1=[4.5,4]&r2=[5.2,4]",fontsize= "13")
plt.xlabel("x (arbitrary unit)",fontsize= "13")
plt.ylabel("y (arbitrary unit)",fontsize= "13")
plt.legend()

# assign empty list for position for the particles
x_p31 = []
y_p31 = []
x_p32 = []
y_p32 = []
#assign initial velocity
v3 = [0, 0, 0, 0]

# r1=[2, 43, r2=[3.5, 4.5]
r3 = np.array([2, 3, 3.5, 4.4], float)

# assign position of particles
x_p31.append(r3[0]) 
y_p31.append(r3[1])
x_p32.append(r3[2]) 
y_p32.append(r3[3])

# calculate v(t+h/2)
v3 += h/2 * f(r3, t_p[0])
# apply verlet algorithm
for t in t_p[1:]:  
    r3 += h * v3
    k = h * f(r3, t)
    v3 += k
    
    x_p31.append(r3[0]) 
    y_p31.append(r3[1])
    x_p32.append(r3[2]) 
    y_p32.append(r3[3])   
    
#plot trajectories of particles
plt.figure()
# plt.plot(t_p, x_p31, '.') 
# plt.plot(t_p, y_p31, '.') 
# plt.plot(t_p, x_p32, '.') 
# plt.plot(t_p, y_p32, '.') 
plt.plot(x_p31, y_p31, '.',label='r1')
plt.plot(x_p32, y_p32, '.', label = 'r2')
plt.title("Tragectories of particles, r1=[2,3]&r2=[3.5,4.4]",fontsize= "13")
plt.xlabel("x (arbitrary unit)",fontsize= "13")
plt.ylabel("y (arbitrary unit)",fontsize= "13")
plt.legend()