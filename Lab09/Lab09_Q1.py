# -*- coding: utf-8 -*-
"""
Lab09_Q1
Time-Dependent Schrodinger Equation
Code the Crank-Nicolson method for the wave equation in the time-dependent Schrodinger equation
in the case of an infinite square well

Done by Juann Jeon
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import scipy

# Q1 (a)
#Constants
L = 1e-8 #Length of well, in meter
m = 9.109e-31 #Mass of electron, in kg
sigma = L/25 
kappa = 500/L 
x_0 = L/5 #Initial starting position
hbar = constants.hbar #Reduced Planck constant
pi = constants.pi #Famous pi
tau = 1e-18 #Time split, in seconds
N = 3000 #Number of time steps
T = N * tau #Total duration
P = 1024 #Number of segments for spatial interval
a = L/P #Step size

#ψ0, the normalization constant
psi_0 = np.sqrt(1/(2 * pi * sigma**2)**0.5)

#spatial_interval = np.linspace(-L/2, L/2, P)
x = np.arange(-L/2, L/2, a)


#Ψ(t), equivalent to equation (9)
psi = np.exp(-(((x[1:] - x_0)**2) / (4 * sigma**2))) * np.exp(1j * kappa * x[1:])
psi *= psi_0

H_D = np.zeros((P - 1, P - 1)) #Discretized Hamiltonian
I = np.identity(P - 1) #Identity matrix
A = -(hbar**2) / (2 * m * a**2)
#Get discretized Hamiltonian
for i in range(len(H_D)):
    for j in range(len(H_D)):
        if i == j:
            H_D[i, j] = -2 * A #Since V = 0, B = -2A
        elif i - 1 == j or i + 1 == j:
            H_D[i, j] = A

L = I + ((1j * tau) / (2 * hbar)) * H_D # this is matrix L in the linear system, L*psi(n+1) = R*psi(n)
L_inv = np.linalg.inv(L) #Inverse of L
R = I - ((1j * tau) / (2 * hbar)) * H_D # this is matrix R in the linear system, L*psi(n+1) = R*psi(n)
psi_current = np.linalg.solve(R, psi) #v = RΨ^n

# Q2 (b) & (c)
# create empty arrays to save results
psi_record_1 = np.array([])
psi_record_2 = np.array([])
psi_record_3 = np.array([])
psi_record_4 = np.array([])
psi_record_5 = np.array([])
psi_record_trajectary = []
energy_record = []
norm = []

# calculate psi(x,t) using crank-nicholson method
for i in range(N):
    # t = 0
    if i == 0:
        psi_record_1 = psi_current        
    # t = T/4
    elif i == ((N/4) - 1):
        psi_record_2 = psi_current  
    # t = T/2
    elif i == ((N/2) - 1):
        psi_record_3 = psi_current   
    # t= 3T/4
    elif i == ((3*N/4) - 1):
        psi_record_4 = psi_current  
    # t = T
    elif i == (N - 1):
        psi_record_5 = psi_current

        
    psi_next = np.linalg.solve(L_inv, psi_current) #Ψ^{n + 1}
    psi_current = np.linalg.solve(R,psi_next) #v = RΨ^n
    
    # find the trajectory of the particle
    integrand = np.conjugate(psi_next) * x[1:] *psi_next
    psi_record_trajectary.append(np.real(scipy.integrate.trapezoid(integrand,x[1:])))
    
    # find the energy of the system
    energy = (np.conjugate(psi_next) *np.linalg.solve(H_D, psi_next))
    energy_record.append(abs(scipy.integrate.trapezoid(energy,x[1:])))
    
    # calculate the probability of finding the particle in the well at a time
    i2 = np.conjugate(psi_next)*psi_next
    norm.append(abs(scipy.integrate.trapezoid(i2,x[1:])))

# plot the real part wavefunction at given time
plt.figure(figsize=(10,7))
plt.plot(x[1:], np.real(np.linalg.solve(L_inv, psi_record_1)), label = "t = 0")
plt.plot(x[1:], np.real(np.linalg.solve(L_inv, psi_record_2)), label = "t = T/4")
plt.plot(x[1:], np.real(np.linalg.solve(L_inv, psi_record_3)), label = "t = T/2")
plt.plot(x[1:], np.real(np.linalg.solve(L_inv, psi_record_4)), label = "t = 3T/4")
plt.plot(x[1:], np.real(np.linalg.solve(L_inv, psi_record_5)), label = "t = T")
plt.title("Real part of ψ(x)", fontsize = "16")
plt.xlabel("L (m)", fontsize = "16")
plt.ylabel("Re{ψ(x)}", fontsize = "16")
plt.legend()

t = np.arange(0,N*tau,tau)

# plot the trajectory of the particle
plt.figure()
plt.plot(t,psi_record_trajectary)
plt.title("The expectation value of particle position", fontsize = "13")
plt.xlabel("time (s)", fontsize = "13")
plt.ylabel("<X>(t) (m)", fontsize = "13")
plt.tight_layout()

# plot the energy of the system
plt.figure()
plt.plot(t,energy_record)
plt.title("Energy of the system", fontsize = "13")
plt.xlabel("time (s)", fontsize = "13")
plt.ylabel("energy (arbitrary unit)", fontsize = "13")
plt.ylim(6.4e16,7.0e16)
plt.tight_layout()

# plot the probability of finding the particle in the well
plt.figure()
plt.plot(t,norm)
plt.title("Probability to find the particle at a time in the well", fontsize = "13")
plt.xlabel("time (s)", fontsize = "13")
plt.ylabel("Probability", fontsize = "13")
plt.ylim(0.990,1.010)
plt.tight_layout()