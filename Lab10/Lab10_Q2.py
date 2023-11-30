# -*- coding: utf-8 -*-
"""
Lab10 Q2
Volume of a 10-dimensional Hypersphere

Done by Juann Jeon
"""
import numpy as np

np.random.seed(815)

N = int(1e6) #Size of the points
dimension = 10 #Dimension we are interested in

def f(x):
    r_square = np.zeros(N) #r^2 = x_1^2 + x_2^2 + ... + x_n^2
    for i in range(dimension):
        r_square += np.square(x[i]) #Adds up x_i
    
    #Condision of f(x_1, x_2, ..., x_n)
    for i in range(N):
        if r_square[i] <= 1:
            r_square[i] = 1
        else:
            r_square[i] = 0
        
    return r_square

sample_x = np.random.uniform(-1, 1, (dimension, N)) #Uniformly distributed random numbers from range -1 to 1

#Calculate the volume of sphere using Monte Carlo Method
print()
print("The Volume of a sphere of unit radius in", dimension, "dimensions is:", ((2**dimension) / N) * np.sum(f(sample_x)))