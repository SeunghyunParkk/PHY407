# -*- coding: utf-8 -*-
"""
Lab08 Question 2
Simulating the shallow water system
Done by Seunghyun Park

This code simulates the 1D shallow water system using FTCS scheme
"""


import numpy as np
import matplotlib.pyplot as plt

L = 1. # width of the system
dx = 0.02 # grid spacing
x = np.arange(0,L+dx,dx) # x values

J = 50 # number of divisions in grid
g = 9.81 # gravitational acceleration

eta_b = 0. # flat bottom topography
H = 0.01 # water column height at rest
u_0 = 0. #boundary condition at x = 0
u_L = 0. #boundary condition at x = L

dt = 0.01 #timestep, in second
A = 0.002 # in metre
mu = 0.5 # in metre
sigma = 0.05 # in metre

epsilon = dt/1000

# given times
t1 = 0.0 # t1 = 0
t2 = 1.0 # t2 = 1.
t3 = 4.0 # t3 = 4.
tend = t3 + epsilon # end time

# define F(u,eta)
def F(u,eta):
    return u**2/2+ g*eta , (eta-eta_b)*u

#create arrays
u = np.zeros(J+1,float)
u[0] = u_0
u[J] = u_L
u[1:J] = 0
u_new = np.zeros(J+1,float)
u_new[0] = u_0
u_new[J] = u_L
eta = H + A*np.exp(-(x-mu)**2/sigma**2) - np.mean(A*np.exp(-(x-mu)**2/sigma**2))
eta_new = np.zeros(J+1,float)


# Main loop
t = 0.0 # start time
water = [] # empty list for eta
while t < tend:
    # flux of the system
    F_u, F_eta =  F(u,eta)
    # calculate the new value of u and eta
    for j in range(1,J):
        u_new[0] = 0. # boundary condition of u, x = 0
        u_new[J] = 0. # boundary condition of u, x = L
        
        eta_new[0] = eta[0] - dt/dx*(F_eta[1] - F_eta[0]) # boundary condition of eta, x = 0
        eta_new[J] = eta[J] - dt/dx*(F_eta[J] - F_eta[J-1]) # boundary condition of eta, x =L
        
        u_new[j] = u[j] - dt/2/dx*(F_u[j+1]-F_u[j-1])
        eta_new[j] = eta[j] - dt/2/dx*(F_eta[j+1]-F_eta[j-1])
        
    
    u = np.copy(u_new)
    eta = np.copy(eta_new)
    t += dt
    # append eta at the given times
    if abs(t-t2) < epsilon:
        water.append(eta)
    if abs(t-t3) < epsilon:
        water.append(eta)
        
# calculate eta at t = 0
eta0 = H + A*np.exp(-(x-mu)**2/sigma**2) - np.mean(A*np.exp(-(x-mu)**2/sigma**2))

# plot the 1D shallow water system at the given times
plt.figure(figsize=(7,4))
plt.plot(x,eta0,  label = 't = 0s')
plt.plot(x,water[0], label = 't = 1s')
plt.plot(x,water[1], label = 't = 4s')
plt.xlabel('x (m)', fontsize = '13')
plt.ylabel('$\eta$ (m)', fontsize = '13')
plt.title('Implementation of 1D shallow water system using the FTCS scheme', fontsize = '13')
plt.legend()





























