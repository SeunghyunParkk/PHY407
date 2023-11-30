# -*- coding: utf-8 -*-
"""
Lab11 Question 1
Simulated Annealing Optimization
Done by Juann Jeon
"""
from math import sqrt, exp, cos, pi
import numpy as np
from random import random, randrange, seed, gauss
import matplotlib.pyplot as plt
#from visual import sphere, curve, display, rate 

seed(10)
N = 25 
R = 0.02 

Tmax = 10.0 
Tmin = 1e-3
tau = 1e3

#Q1 a)
# Function to calculate the magnitude of a vector 
def mag(x): 
    return sqrt(x[0]**2 + x[1]**2) 

# Function to calculate the total length of the tour 
def distance(): 
    s = 0.0 
    for i in range(N): 
        s += mag(r[i + 1] - r[i]) 
    return s 

# Choose N city locations and calculate the initial distance 
r = np.empty([N + 1, 2], float) 
for i in range(N): 
    r[i, 0] = random() 
    r[i, 1] = random() 
r[N] = r[0] 
D = distance()
 
# Set up the graphics 
#display(center = [0.5, 0.5]) 
#for i in range(N): 
    #sphere(pos = r[i], radius=R) 
#l = curve(pos = r, radius = R/2)
 
# Main loop 
t = 0 
T = Tmax 
seed()
while T > Tmin: 
    t += 1 
    T = Tmax * exp(-t / tau) # Cooling
    #if t % 100 == 0: 
        #l.pos = r # Update the visualization every 100 moves 
        #rate(25) 

    # Choose two cities to swap and make sure they are distinct 
    i, j = randrange(1, N), randrange(1, N) 
    while i == j: 
        i, j = randrange(1, N), randrange(1, N) 
        
    # Swap them and calculate the change in distance 
    oldD = D 
    r[i, 0], r[j, 0] = r[j, 0], r[i, 0] 
    r[i, 1], r[j, 1] = r[j, 1], r[i, 1] 
    D = distance() 
    deltaD = D - oldD 

    # If the move is rejected, swap them back again 
    if random() >= exp(-deltaD / T): 
        r[i, 0], r[j, 0] = r[j, 0], r[i, 0] 
        r[i, 1], r[j, 1] = r[j, 1], r[i, 1] 
        D = oldD 

print("τ:", tau)
print("D:", D)

x, y = r.T
plt.figure()
plt.plot(x, y, ".", markersize = 12)
plt.plot(x, y)
plt.title("Solution to TSM with D = {}".format(D), fontsize = "14")
plt.xlabel("x", fontsize = "14")
plt.ylabel("y", fontsize = "14")

#Q1 b)
seed(815)
print()
tau = 1e4

# define function f(x,y)
def f(x, y):
    return x**2 - cos(4 * pi * x) + (y - 1)**2

x = 2
y = 2
newf = f(x, y)
t = 0 
T = Tmax 
record_x = [x]
record_y = [y]
#  program to confirm (x,y) = (0,1) using simulated annealing
while T > Tmin:
    t += 1
    T = Tmax * exp(-t / tau) # cooling rate
    
    oldx = x
    oldy = y
    oldf = newf
    
    x += gauss(0, 1) #gaussian distribution with mean 0, standard deviation 1
    y += gauss(0, 1) #gaussian distribution with mean 0, standard deviation 1
    newf = f(x, y)
    
    deltaf = newf - oldf
    
    if random() >= exp(-deltaf / T): 
        x = oldx
        y = oldy
        newf = oldf
        
        record_x.append(x)
        record_y.append(y)
        
print("The global min happens at x =", x)
print("The global min happens at y =", y)
print("The global min of f(x, y) is", newf)

# plot x, y
plt.figure()
plt.plot(record_x, ".", label = "x")
plt.plot(record_y, ".", label = "y")
plt.title("Plot of (x, y) → (x + δx, y + δy)", fontsize = "14")
plt.xlabel("t", fontsize = "14")
plt.ylabel("Values of x and y", fontsize = "14")
plt.legend()
# plt.savefig('b.png')



#Q1 c)
seed(1611)
print()
tau = 1e4
# define function f(x,y), equation (5)
def f_hard(x, y):
    if x > 0 and x < 50 and y > -20 and y < 20:
        return cos(x) + cos(sqrt(2) * x) + cos(sqrt(3) * x) + (y - 1)**2
    else:
        return 100000000

x_hard = 2
y_hard = 2
newf_hard = f_hard(x_hard, y_hard)
t = 0 
T = Tmax 
record_x_hard = [x_hard]
record_y_hard = [y_hard]

#  program to confirm (x,y) = (16,1) using simulated annealing
while T > Tmin:
    t += 1
    T = Tmax * exp(-t / tau)
    
    oldx_hard = x_hard
    oldy_hard = y_hard
    oldf_hard = newf_hard
    
    x_hard += gauss(0, 1) #gaussian distribution with mean 0, standard deviation 1
    y_hard += gauss(0, 1) #gaussian distribution with mean 0, standard deviation 1
    newf_hard = f_hard(x_hard, y_hard)
    
    deltaf_hard = newf_hard - oldf_hard
    
    if random() >= exp(-deltaf_hard / T): 
        x_hard = oldx_hard
        y_hard = oldy_hard
        newf_hard = oldf_hard
        
        record_x_hard.append(x_hard)
        record_y_hard.append(y_hard)
        
print("The global min happens at x =", x_hard)
print("The global min happens at y =", y_hard)
print("The global min of f(x, y) is", newf_hard)

# plot x,y
plt.figure()
plt.plot(record_x_hard, ".", label = "x")
plt.plot(record_y_hard, ".", label = "y")
plt.title("Scatter plot of (x, y) → (x + δx, y + δy)", fontsize = "14")
plt.xlabel("y", fontsize = "14")
plt.ylabel("Values of x and y", fontsize = "14")
plt.legend()
# plt.savefig('c.png')