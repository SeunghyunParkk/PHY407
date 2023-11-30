# -*- coding: utf-8 -*-
"""
Lab10 Question 3
Importance Sampling

Calculating given integral using two different method: 
Mean value method & Importance sampling method.

Done by Seunghyun park
"""

import numpy as np
import random
import matplotlib.pyplot as plt

# Q3 (a)
# set seed
random.seed(815)

# define the integrand
def f(x):
    return x**(-1/2)/(1 + np.exp(x))

# number of sample points
N = 10000

# define the mean value method
def mean_value(f,a,b,N):
    sum_f = 0 # initial value of sum f(x)
    for i in range(N):
        x = a + random.random() * (b-a) # generate random number between upper bound and lower bound of integral
        sum_f += f(x) 
    value = (b-a)/N* sum_f
    return value

# create a empty list for integral value calculated using mean value method 
mv_int = []
# repeat calculation 100 times
for i in range(0,100):
    mv_int.append(mean_value(f,0,1,1000))

# print the value of integral calculated using mean value method
print("The value of integral calculated using mean value method with 10,000 sample points is",mv_int)
print()


# Q3 (b)

# define weighting function
def weight(x):
    return x**(-1/2)

# define the importance sampling method
def ism(f,a,b,N):
    sumfw = 0 
    for i in range(N):
        # integral of probability distribution = x^(1/2) = z
        # random number z^(2) = x
        x = (random.random() * (b-a))**2 
        sumfw += f(x) / weight(x)
        
    # integral of weighting function = 2
    value = 1 / N * sumfw * 2
    return value

# create a empty list for integral value calculated using importance sampling method
is_int = []
# repeat calculation 100 times
for i in range(0,100):
    is_int.append(ism(f,0,1,1000))
print("The value of integral calculated using importance sampling method with 10,000 sample points is",is_int)
print()

# Q3 (c)

# make histograms of the integral values for each method
plt.figure(figsize = (8,4))
plt.hist(mv_int, 10, range=[0.8, 0.88],edgecolor="k")
plt.title("Histogram of the integral value using Mean Value Method", fontsize = "15")
plt.xlabel("Value of Integral", fontsize = "13")
plt.ylabel("Frequency", fontsize = "13")
plt.show()
plt.figure(figsize = (8,4))
plt.hist(is_int, 10, range=[0.8, 0.88],edgecolor="k")
plt.title("Histogram of the integral value using Importance Sampling Method", fontsize = "15")
plt.xlabel("Value of Integral", fontsize = "13")
plt.ylabel("Frequency", fontsize = "13")
plt.show()



