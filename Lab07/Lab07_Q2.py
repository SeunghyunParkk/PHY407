# -*- coding: utf-8 -*-
"""
Lab07_Q2
Hydrogen Atom
Done by Juann Jeon
"""

import numpy as np
from scipy import constants
import scipy
import matplotlib.pyplot as plt

#Q2 (a)
#Constants 
m = constants.m_e #Mass of electron 
hbar = constants.hbar #Planck's constant over 2*pi 
e = constants.eV #Electron charge 
a = 5.2918e-11 #Bohr radius 
epsilon_0 = constants.epsilon_0 #Vacuum permitivity
pi = constants.pi #Famous 3.14...
h = 0.002 * a #Time step


n = 1
l = 0
r_inf = np.arange(h, 20 * a, h)

#Potential energy for hydrogen atom
def V(x): 
    return (-e**2)/(4 * pi * epsilon_0 * x)

# define a function that return dR/dr, dS/dr
def f(r, x, E): 
    R = r[0] 
    S = r[1] 
    fR = S 
    fS = ((2 * m / hbar**2) * (V(x) - E) * R) + (l * (l + 1) * R / x**2) - (2 * S / x)
    return np.array([fR, fS], float) 

#Calculate the wavefunction for a particular energy 
def solve(E): 
    R = 0.0 
    S = 1.0 
    r = np.array([R, S], float) 
    for x in r_inf: 
        k1 = h * f(r, x, E) 
        k2 = h * f(r + 0.5 * k1, x + 0.5 * h, E) 
        k3 = h * f(r + 0.5 * k2, x + 0.5 * h, E) 
        k4 = h * f(r + k3, x + h, E) 
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6 
    return r[0] 


# Q2 (b)
#Main program to find the energy using the secant method 
# ground state with l = 0
E1 = -15 * e / n**2
E2 = -13 * e / n**2
R2 = solve(E1) 

# target convergence
target = e / 1000 

while abs(E1 - E2) > target: 
    R1, R2 = R2, solve(E2) 
    E1, E2 = E2, E2 - R2 * (E2 - E1) / (R2 - R1) 
    
print("Ground state with l = 0, E = ", E2 / e, "eV", sep = '')

# first excited state with n = 2 , l = 0
n = 2
l = 0
#Main program to find the energy using the secant method 
E1 = -15 * e / n**2
E2 = -13 * e / n**2
R2 = solve(E1) 



while abs(E1 - E2) > target: 
    R1, R2 = R2, solve(E2) 
    E1, E2 = E2, E2 - R2 * (E2 - E1) / (R2 - R1) 

# print energy of First excited state with l = 0
print("First excited state with l = 0, E = ", E2 / e, "eV", sep = '')

# first excited state with n = 2 , l = 1
n = 2
l = 1
#Main program to find the energy using the secant method 
E1 = -15 * e / n**2
E2 = -13 * e / n**2
R2 = solve(E1) 


while abs(E1 - E2) > target: 
    R1, R2 = R2, solve(E2) 
    E1, E2 = E2, E2 - R2 * (E2 - E1) / (R2 - R1) 
    
# print energy of First excited state with l = 1
print("First excited state with l = 1, E = ", E2 / e, "eV", sep = '')


# Q2(c)
r_inf = np.arange(h, 20 * a, h)

# define a function tha returns R as a function of r
def rlist(E): 
    R_list = []
    R = 0.0 
    S = 1.0 
    r = np.array([R, S], float) 
    for x in r_inf: 
        k1 = h * f(r, x, E) 
        k2 = h * f(r + 0.5 * k1, x + 0.5 * h, E) 
        k3 = h * f(r + 0.5 * k2, x + 0.5 * h, E) 
        k4 = h * f(r + k3, x + h, E) 

        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6   
        R_list.append(r[0])
    return R_list

# define a function that returns normalized absolute(R)^2
def normalize_R(R, x):
    y = np.abs(R)**2
    norm_factor = scipy.integrate.simps(y,x)
    return y/(norm_factor)

# obtain R values for n =1, l=0
n = 1
l = 0
E1 = -15 * e / n**2
r10 = rlist(E1)


# obtain R values for n =2, l=0
n = 2
l = 0
E1 = -15 * e / n**2
r20 = rlist(E1)

# obtain R values for n =2, l=1
n = 2
l = 1
E1 = -15 * e / n**2
r21 = rlist(E1)

# plot normalized R
plt.figure()
plt.plot(r_inf[:2000], normalize_R(r10, r_inf)[:2000], label = "n=1,l=0")
plt.title("Normalized Eigenfunction R, n=1,l=0", fontsize = "16")
plt.xlabel("r (m)", fontsize = "16")
plt.ylabel("normalized $|R(r)|^2$", fontsize = "16")
plt.tight_layout()
plt.figure()
plt.plot(r_inf[:5000], normalize_R(r20, r_inf)[:5000], label = "n=2,l=0")
plt.title("Normalized Eigenfunction R, n=2,l=0", fontsize = "16")
plt.xlabel("r (m)", fontsize = "16")
plt.ylabel("normalized $|R(r)|^2$", fontsize = "16")
plt.tight_layout()
plt.figure()
plt.plot(r_inf[:5000], normalize_R(r21, r_inf)[:5000], label = "n=2,l=1")
plt.title("Normalized Eigenfunction R, n=2,l=1", fontsize = "16")
plt.xlabel("r (m)", fontsize = "16")
plt.ylabel("normalized $|R(r)|^2$", fontsize = "16")
plt.tight_layout()