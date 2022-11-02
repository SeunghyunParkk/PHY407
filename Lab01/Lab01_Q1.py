'''
PHY407 Lab 1 Question 1

This .py allows user to see the displacement and velocity of the classical 
Newtonian gravitational force and the modern version of with general relativity.
Here, we assumed Mercury is orbiting around the Sun and they are the only two
masses in the system.

By Seunghyun Park & Juann Jeon
All parts of codes are done together
'''

import numpy as np
import matplotlib.pyplot as plt

# Q1 (c)

# Constants
x_i = 0.47 #Initial x position of Mercury, in AU
y_i = 0.0 #Initial y position of Mercury, in AU
v_xi = 0.0 #Initial x velocity of Mercury, in AU/yr
v_yi = 8.17 #Initial y velocity of Mercury, in AU/yr
dt = 0.0001 #Time step, in yr
t = np.arange(0.0, 1.0 , dt) #Time interval of 1 year
M_s = 2e30 #Mass of the sun, in kg
G = 39.5 / M_s #Gravitational constant, in AU^3 M_s^-1 yr^-2
a = 0.01 #Some constant, in AU^2


#np array to store values of x and y position and x and y velocity for Q1 c)
#velocity of mercury in x-direction. Create np array with same length as t, with all indexes 0
v_x = np.zeros(len(t)) 
v_x[0] = v_xi #Assign initial value
#position of mercury in x-direction. Create np array with same length as t, with all indexes 0 for x position
x = np.zeros(len(t)) 
x[0] = x_i #Assign initial value
#velocity of mercury in y-direction. Create np array with same length as t, with all indexes 0
v_y = np.zeros(len(t)) 
v_y[0] = v_yi #Assign initial value
#position of mercury in y-direction. Create np array with same length as t, with all indexes 0
y = np.zeros(len(t)) 
y[0] = y_i #Assign initial value
#orbital radius of mercury. Create np array with same length as t, with all indexes 0
r = np.zeros(len(t)) 
r[0] = (x[0]**2+y[0]**2)**0.5 #Assign initial value of orbital radius of Mercury

for i in range(len(t)-1):
    v_x[i+1] = v_x[i] - G * M_s * x[i] / r[i]**(3) * dt 
    x[i+1] = x[i] + v_x[i+1] * dt
    v_y[i+1] = v_y[i] - G * M_s * y[i] / r[i]**(3) * dt 
    y[i+1] = y[i] + v_y[i+1]*dt
    r[i+1] = (x[i+1]**2 + y[i+1]**2)**0.5

# Angular Momentum Conservation
L = np.zeros(len(t)) # Create np array with same length as t, with all indexes 0
# Mass of Mercury is not given in the lab instruction
# Calculate Angular Momentum of Mercury divided by its mass
L[0] = x[0]*v_y[0]-y[0]*v_x[0] # Assign initial value of angluar momentum
for i in range(len(t)):
    L[i] = x[i]*v_y[i]-y[i]*v_x[i]

plt.figure()
plt.plot(t, v_x, label ='v_x')
plt.plot(t,v_y, label ='v_y')
plt.title('Velocity of Mercury orbiting around the Sun')
plt.xlabel('Time (year)')
plt.ylabel('Velocity (AU/year)')
plt.grid()
plt.legend()


plt.figure()
plt.plot(x,y)
plt.scatter(0,0, c ='red', label = 'The Sun')
plt.title('Position of Mercury orbiting around the Sun')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.axis('scaled')
plt.grid()

plt.figure()
plt.plot(t,L)
plt.ylim(-10,10)
plt.title('Angular Momentum of Mercury divided by its Mass')
plt.xlabel('Time (year)')
plt.ylabel('Angular Momentum divided by Mass ($AU^{2}/year$)')
plt.grid()



# Q1 (d)

#np array to store values of x and y position and x and y velocity for Q1 d)
#velocity of mercury in x-direction. Create np array with same length as t, with all indexes 0
v_x2 = np.zeros(len(t)) 
v_x2[0] = v_xi #Assign initial value
#position of mercury in x-direction. Create np array with same length as t, with all indexes 0
x2 = np.zeros(len(t)) #Create np array with same length as t, with all indexes 0
x2[0] = x_i #Assign initial value
#velocity of mercury in y-direction. Create np array with same length as t, with all indexes 0
v_y2 = np.zeros(len(t)) 
v_y2[0] = v_yi #Assign initial value
#position of mercury in y-direction. Create np array with same length as t, with all indexes 0
y2 = np.zeros(len(t)) 
y2[0] = y_i #Assign initial value
#orbital radius of mercury. Create np array with same length as t, with all indexes 0
r2 = np.zeros(len(t)) 
r2[0] = (x2[0]**2+y2[0]**2)**0.5 #Assign initial value of orbital radius of Mercury

#Calculate x and y position and x and y velocity for Q1 d) using Euler-Cromer method
for i in range(len(t)-1):
    v_x2[i+1] = v_x2[i] - G * M_s * x2[i] / r2[i]**3 * dt * (1 + a/(r2[i] ** 2)) #x velocity
    x2[i+1] = x2[i] + v_x2[i+1] * dt #x position
    v_y2[i+1] = v_y2[i] - G * M_s * y2[i] / r2[i]**3 * dt * (1 + a/(r2[i] ** 2)) #y velocity
    y2[i+1] = y2[i] + v_y2[i+1] * dt #y position
    r2[i+1] = (x2[i+1]**2 + y2[i+1]**2)**0.5

plt.figure()
plt.plot(t, v_x2, label ='v_x')
plt.plot(t,v_y2, label ='v_y')
plt.title('Velocity of Mercury orbiting around the Sun with general relativity')
plt.xlabel('Time (year)')
plt.ylabel('Velocity (AU/year)')
plt.grid()
plt.legend()

plt.figure()
plt.plot(x2,y2)
plt.scatter(0,0, c ='red', label = 'The Sun')
plt.title('Position of Mercury orbiting around the Sun with general relativity')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.axis('scaled')
plt.grid()