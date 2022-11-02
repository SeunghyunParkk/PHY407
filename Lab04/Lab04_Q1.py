# -*- coding: utf-8 -*-
'''
PHY407 Lab04

Done by Seunghyun Park

Q1

'''
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve
from numpy.random import rand
from SolveLinear import GaussElim
from SolveLinear import PartialPivot
import time

#(a)
# Equation 6.2, Assign A and v
A = np.array([[2,1,4,1],[3,4,-1,-1],[1,-4,1,5],[2,-2,1,3]],float)
v = np.array([-4,3,9,7],float)
#Solve the equation with partialpivot
x  = PartialPivot(A,v)
print()
print("Answer to Equation 6.2 using Partial Pivot function:", x)



#(b)
#Assign N value
NS = np.arange(5,300)
#Assign empty lists for time
t_G = []
t_P = []
t_L = []
#Assign empty lists for error
error_G = []
error_P = []
error_L = []

#Measure the time it takes and error using for loop
for N in NS:
    #generate random matrix
    v = rand(N)
    A = rand(N,N)
    
    #Measure time and error of gaussian elimination
    t1 = time.time()
    x_G = GaussElim(A,v)
    t2 = time.time()
    t_G.append(t2-t1)
    vsol_G = np.dot(A,x_G)
    error_G.append( np.mean(abs(v-vsol_G)))
    
    #Measure time and error of partialpivot
    t3 = time.time()
    x_P = PartialPivot(A,v)
    t4 = time.time()
    t_P.append(t4-t3)
    vsol_P = np.dot(A,x_P)
    error_P.append( np.mean(abs(v-vsol_P)))
    
    #Measure time and error of LU decomposition
    t5 = time.time()
    x_L = solve(A,v)
    t6 = time.time()
    t_L.append(t6-t5)
    vsol_L = np.dot(A,x_L)
    error_L.append( np.mean(abs(v-vsol_L)))

plt.figure()
plt.loglog(NS,t_G, label = 'GaussElim')
plt.loglog(NS,t_P, label = 'PartialPivot')
plt.loglog(NS,t_L, label = 'Lu')
plt.title("The time takes to solve for x in logarithmic scale", fontsize = '14')
plt.xlabel("N in logarithmic scale", fontsize = '14')
plt.ylabel("time in logarithmic scale (s)", fontsize = '14')
plt.legend()
plt.figure()
plt.loglog(NS,error_G,label = 'GaussElim' )
plt.loglog(NS,error_P, label = 'PartialPivot')
plt.loglog(NS,error_L, label = 'Lu')
plt.title("Error in logarithmic scale", fontsize = '14')
plt.xlabel("N in logarithmic scale", fontsize = '14')
plt.ylabel("error in logarithmic scale", fontsize = '14')
plt.legend()