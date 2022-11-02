# -*- coding: utf-8 -*-
"""
Lab07_Q1
Space Garbage
Done by Seunghyun Park
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as pc
import time

# Q1 (a)
# Constant values
a = pc.physical_constants['Bohr radius'][0]
E0 = pc.physical_constants['Rydberg constant times hc in eV'][0]
em = pc.m_e # electron mass
hbar = pc.hbar
ec = pc.e # elementary charge
e0 = pc.epsilon_0

# This code is from Newman_8-8.py
def rhs(r):
    """ The right-hand-side of the equations
    INPUT:
    r = [x, vx, y, vy] are floats (not arrays)
    note: no explicit dependence on time
    OUTPUT:
    1x2 numpy array, rhs[0] is for x, rhs[1] is for vx, etc"""
    M = 10.
    L = 2.

    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]

    r2 = x**2 + y**2
    Fx, Fy = - M * np.array([x, y], float) / (r2 * np.sqrt(r2 + .25*L**2))
    return np.array([vx, Fx, vy, Fy], float)


ftsz = 16


# %% This next part adapted from Newman's odesim.py --------------------------|
# initial time
a = 0.0
# end time
b = 10.0
# number of steps
N = 1000 
# step size
h = (b-a)/N

# time
tpoints = np.arange(a, b, h)

# Assign empty list for x,y,vx,vy
x1points = []
vx1points = []  # the future dx/dt
y1points = []
vy1points = []  # the future dy/dt

# targer error per second
te = 10**(-6)


# define a function that calculates position using adaptive step size 
def varystep(r,T,h):
    # perform two steps of size h 
    rh1 = r
    k11 = h*rhs(rh1)  
    k21 = h*rhs(rh1 + 0.5*k11)  
    k31 = h*rhs(rh1 + 0.5*k21)
    k41 = h*rhs(rh1 + k31)
    rh1 = rh1 + (k11 + 2*k21 + 2*k31 + k41)/6
    k11 = h*rhs(rh1)  
    k21 = h*rhs(rh1 + 0.5*k11)  
    k31 = h*rhs(rh1 + 0.5*k21)
    k41 = h*rhs(rh1 + k31)
    rh1 = rh1 + (k11 + 2*k21 + 2*k31 + k41)/6
    
    # perform one step of size 2h
    rh2 = r
    k12 = 2*h*rhs(rh2)  
    k22 = 2*h*rhs(rh2 + 0.5*k12)  
    k32 = 2*h*rhs(rh2 + 0.5*k22)
    k42 = 2*h*rhs(rh2 + k32)
    rh2 = rh2 + (k12 + 2*k22 + 2*k32 + k42)/6
    
    # Calculate error 
    ex = 1 / 30 * (rh1[0] - rh2[0])
    ey = 1 / 30 * (rh1[2] - rh2[2])
    # calculate ratio of the target accuracy and the error terms
    p = h * te / np.sqrt(ex**2+ey**2)
    
    
    hnew = h*p**(0.25)
    tnew = t+ 2*h
    
    # if the ratio is larger than 1, return the value with one step of size 2h
    if p > 1 :
        return rh2, tnew , hnew
    # if the ratio is smaller than 1, repeat the current step again with new step size
    else:
        return r, t, hnew
    

# assign r and initial time
r1 = np.array([1., 0., 0., 1.], float)
t = 0

# calculate x,y,vx,vy using varystep function
while t < b:
    x1points.append(r1[0])
    vx1points.append(r1[1])
    y1points.append(r1[2])
    vy1points.append(r1[3])
    r1,t,h = varystep(r1,t,h)


# initial time
a = 0.0
# end time
b = 10.0
# number of steps
N = 10000 
# step size
h = (b-a)/N

# time
tpoints = np.arange(a, b, h)
# empty lists for x,y,vx,vy
xpoints = []
vxpoints = []  
ypoints = []
vypoints = []  

# below: ordering is x, dx/dt, y, dy/dt
r = np.array([1., 0., 0., 1.], float)
# calculate the solutions using fixed step size
for t in tpoints:
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    k1 = h*rhs(r)  
    k2 = h*rhs(r + 0.5*k1)  
    k3 = h*rhs(r + 0.5*k2)
    k4 = h*rhs(r + k3)
    r += (k1 + 2*k2 + 2*k3 + k4)/6


# plot the trajectories obtained by adaptive & fixed  step size
plt.figure()
plt.plot(xpoints, ypoints, ':', label = 'fixed')
plt.plot(x1points, y1points, ':', label = 'adaptive')
plt.xlabel("$x$ (m)", fontsize=ftsz)
plt.ylabel("$y$ (m)", fontsize=ftsz)
plt.title('Trajectory of Space Garbage', fontsize=ftsz)
plt.axis('equal')
plt.grid()
plt.legend(loc="best", prop={'size': 12})
plt.tight_layout()
plt.show()


# Q1 (b)

# initial time
a = 0.0
# end time
b = 10.0
# number of steps
N = 10000  
# step size
h = 0.001

# empty list for x,y,vx,vy
x1points = []
vx1points = []  # the future dx/dt
y1points = []
vy1points = []  # the future dy/dt

# targer error per second
te = 10**(-6)

# measure the time it takes for adaptive step size method
r = np.array([1., 0., 0., 1.], float)
t = 0
t1 = time.time()
while t < b:
    x1points.append(r[0])
    vx1points.append(r[1])
    y1points.append(r[2])
    vy1points.append(r[3])
    r,t,h = varystep(r,t,h)
t2 = time.time()
time_a = t2-t1

# initial time
a = 0.0
# end time
b = 10.0
# number of steps
N = 10000  
# step size
h = 0.001
#time
tpoints = np.arange(a, b, h)
# empty list for x,y,vx,vy
xpoints = []
vxpoints = []  # the future dx/dt
ypoints = []
vypoints = []  # the future dy/dt

# measure the time it takes for fixed step size method
r = np.array([1., 0., 0., 1.], float)
t3 = time.time()
for t in tpoints:
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    k1 = h*rhs(r)  # all the k's are vectors
    k2 = h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
    k3 = h*rhs(r + 0.5*k2)
    k4 = h*rhs(r + k3)
    r += (k1 + 2*k2 + 2*k3 + k4)/6
t4 = time.time()
time_o = t4-t3

#print the times
print()
print("Time it takes for fixed step size", time_o,"seconds")
print()
print("Time it takes for adaptive step size", time_a,"seconds")

# Q1 (c)
# initial time
a = 0.0
# end time
b = 10.0
# number of steps
N = 10000  
# step size
h = 0.001
# empty list for x,y,vx,vy
x1points = []
vx1points = []  # the future dx/dt
y1points = []
vy1points = []  # the future dy/dt

# targer error per second
te = 10**(-6)

# step size as a function of time
r = np.array([1., 0., 0., 1.], float)
t = 0
# empty list of time and step size
tpoint = []
hpoint = []
#empty list of velocity
v= []
while t < b:
    tpoint.append(t)
    hpoint.append(h)
    x1points.append(r[0])
    vx1points.append(r[1])
    y1points.append(r[2])
    vy1points.append(r[3])
    v.append((r[1]**2+r[3]**2)**0.5)
    r,t,h = varystep(r,t,h)

#plot velocity of the system
plt.figure()
plt.plot(tpoint,v)
plt.title("Velocity of space garbage", fontsize=ftsz)
plt.xlabel("time (s)", fontsize=ftsz)
plt.ylabel("velocity (m/s)", fontsize=ftsz)
#plot adaptive step size as a function of time
plt.figure()
plt.plot(tpoint,hpoint)
plt.title("Size of Time Step as a Function of Time", fontsize=ftsz)
plt.xlabel("time (s)", fontsize=ftsz)
plt.ylabel("h (s)", fontsize=ftsz)