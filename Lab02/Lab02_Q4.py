# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 23:05:21 2022

@author: Juann Jeon
Lab02_Q4
"""

import numpy as np
import matplotlib.pyplot as plt

#Define function p(u) of order 8
def p(u):
    return (1 - u)**8

#Expended form of p(u), mathematically q(u) is same as p(u)
def q(u):
    return 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8

C = 1e-16 #Error constant
N = 500 #Size of the array for u
u = np.linspace(0.98, 1.02, num = N) #Generate number close to 1 with length N

#Q4a)
plt.figure()
plt.plot(u, p(u), label = 'p(u)')
plt.plot(u, q(u), label = 'q(u)')
plt.title('Computation of p(u) and q(u)')
plt.xlabel('u')
plt.ylabel('output')
plt.legend()


#Q4b)

# The standard deviation when x = p(u) - q(u)
std_pq = np.std(p(u) - q(u)) #Calculate much accurate standard deviation
sigma_pq = C * (N**0.5) * np.sqrt(np.mean(np.square(p(u) - q(u)))) #Calculate standard deviation using Eq (3)
sigma_pq_normal = np.sqrt(np.mean(np.square(p(u) - q(u)))) #Calculate standard deviation without using np.std
digit_pq = C * (N**0.5)
print("standard deviation of p(u) - q(u) using np.std():", std_pq) 
print("Using Equation (3):", sigma_pq)
print("standard deviation of p(u) - q(u) without C * sqrt(N):", sigma_pq_normal) 
print("Digits of C * sqrt(N):", digit_pq)


# The standard deviation when x is the sum of a series of terms, where each term is power
# of u or of (1-u).
# then x = 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8 - (1-u)**8

# since we are doing when u is very close to 1
n = int(N/2) 
x = [1,- 8*u[n], 28*u[n]**2, - 56*u[n]**3, 70*u[n]**4, -56*u[n]**5, 28*u[n]**6, -8*u[n]**7, 1*u[n]**8, -(1 - u[n])**8]
x_sq = np.square(x)
x_mean_sq = np.sum(x_sq)/len(x) # Calculate the mean of squared value x 
sigma_pq_e = C*np.sqrt(len(x))*np.sqrt(x_mean_sq) # Calculate the standard deviation
print("Standard deviation if we treat p(u)-q(u) as the sum of series:", sigma_pq_e)

plt.figure()
plt.plot(u, p(u) - q(u), label = 'p(u) - q(u)')
plt.title('p(u) - q(u) vs u')
plt.xlabel('u')
plt.ylabel('p(u) - q(u)')

plt.figure()
plt.hist(p(u) - q(u))
plt.title('Histogram of p(u) - q(u)')
plt.xlabel('p(u) - q(u)')
plt.ylabel('Frequency')


#Q4c)
# Assign u starting from 0.980 to 0.984 with 50 points
# the number of points were 500 when u = [0.980,1.020]
# However the interval becomes 10 times smaller than original
# N becomes 10 times smaller too
u = np.linspace(0.98, 0.984, 50)
error = abs(p(u) - q(u))/abs(p(u))

plt.figure()
plt.plot(u, error, label = 'Error')
plt.title('Fractional Error')
plt.xlabel('u')
plt.ylabel('Error in decimal (1.0 = 100%)')

# Find a point when the error approaches 100% using for loop
for i in range(0,len(error)):
    if error[i] > 1:
        print("The error approaches 100% at",u[i])
        break


#Q4d)
# Assign u same as part (a)
u = np.linspace(0.98, 1.02, 500)
f = u**8/((u**4)*(u**4))

plt.figure()
plt.plot(u, f - 1, label = 'Error')
plt.title('f - 1 vs u')
plt.xlabel('u')
plt.ylabel('f - 1')

print("The standard deviation of f using numpy.std is ",np.std(f))

# error estimation 
sigma_round = np.sqrt(2)*C*u
print("The mean of estimated error using equation (4.5) on p.131 of the textbook is", sigma_round[250])