"""
Lab10_Q1
This program simulates Brownian motion in the presence of walls
Note that the physical behaviour would be to stick to walls,
which is the purpose of Q1a.
Author: Nico Grisouard, University of Toronto
Modified by: Juann Jeon
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import random

random.seed(95)

#For Q1a)
def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = right, 3 = left"""
    
    blocked = [] #Blocks for going certain direction if there is wall at there
    if y == Lp:
        blocked.append(0) #Blocks going up
    if y == 0:
        blocked.append(1) #Blocks going down
    if x == Lp:
        blocked.append(2) #Blocks going right
    if x == 0:
        blocked.append(3) #Blocks going left
     
    direction = random.choice([i for i in range(0, 4) if i not in blocked]) #Decide which direction to go
    
    if direction == 0: #Move up
        y += 1
    elif direction == 1: #Move down
        y -= 1 
    elif direction == 2: #Move right
        x += 1 
    elif direction == 3: #Move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y

#For Q1b)
def nextmove_DLA(x, y, check, anchored):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = right, 3 = left"""
    direction = random.randint(0, 3) #Decide which direction to go
    
    if y == Lp or y == 0 or x == Lp or x == 0 or check[x, y + 1] == True or check[x, y - 1] == True or check[x + 1, y] == True or check[x - 1, y] == True:
        check[x, y] = True
        anchored = True
    else:
        if direction == 0:
            y += 1 #Move up
        elif direction == 1:
            y -= 1 #Move down
        elif direction == 2:
            x += 1 #Move right
        elif direction == 3:
            x -= 1 #Move left
        else:
            print("error: direction isn't 0-3")
        
    return x, y, check, anchored

font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)

# %% main program starts here ------------------------------------------------|
#Q1a)

plt.ion()

Lp = 101  #Size of domain
Nt = 5000  #Number of time steps

centre_point = (Lp-1)//2  #Middle point of domain
xp = centre_point
yp = centre_point

#Arrays to record the trajectory of the particle
xp_record = [xp]
yp_record = [yp]

for i in range(1, Nt):
    xp, yp = nextmove(xp, yp) #Move to next point
    xp_record.append(xp) #Record x position
    yp_record.append(yp) #Record y position
    
plt.figure()
plt.plot(xp_record)
plt.title("Particle's i position vs t")
plt.xlabel("time step")
plt.ylabel("i position")

plt.figure()
plt.plot(yp_record)
plt.title("Particle's j position vs t")
plt.xlabel("time step")
plt.ylabel("j position")

plt.figure(figsize=(7,7))
plt.plot(xp_record, yp_record)
plt.title("Trajectory of a Particle in Brownian Motion")
plt.xlabel("i")
plt.ylabel("j")
plt.xlim(-1,101)
plt.ylim(-1,101)
# plt.axis('scaled')

#Q1b)
anchored = False #Boolean expression to check if particle is anchored
check = np.zeros((Lp + 1, Lp + 1), dtype = bool) #Boolean grid check board to record if particle is anchored at given coordinate
count = 0 #Count the number of particles

plt.figure(figsize=(7,7))
plt.xlim([0, Lp]) #Limits the x size of the plot
plt.ylim([0, Lp]) #Limits the y size of the plot

while check[centre_point, centre_point] == False:
    xp, yp, check, anchored = nextmove_DLA(xp, yp, check, anchored)
    if anchored == True:
        plt.plot(xp, yp, ".")
        anchored = False
        xp = centre_point
        yp = centre_point
        count += 1
    
plt.title("DLA run for {} particles".format(count))
plt.xlabel("i")
plt.ylabel("j")
plt.axis('scaled')
