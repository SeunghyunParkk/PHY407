# '''
# Question 2
# By Seunghyun Park & Juann Jeon
# All parts of codes are done together
# '''

import numpy as np
import matplotlib.pyplot as plt


# Question 2 (a)

dt = 0.0001 #Time step, in yr
t = np.arange(0.0, 10.0 , dt) #Time interval of 10 year
M_s = 2e30 #Mass of the sun, in kg
G = 39.5 / M_s #Gravitational constant, in AU^3 M_s^-1 yr^-2
M_j = 10**(-3) * M_s #Mass of Jupiter


#np array to store values of x and y position and x and y velocity for Q1 c)

v_jx = np.zeros(len(t)) #Velocity of Jupiter in x-direction. Create np array with same length as t, with all indexes 0
x_j = np.zeros(len(t)) #Position of Jupiter in x-direction. Create np array with same length as t, with all indexes 0
v_jy = np.zeros(len(t)) #Velocity of Jupiter in y-direction. Create np array with same length as t, with all indexes 0
y_j = np.zeros(len(t)) #Position of Jupiter in y-direction. Create np array with same length as t, with all indexes 0
r_j = np.zeros(len(t)) #Radius of Jupiter's orbit. Create np array with same length as t, with all indexes 0

x_j[0] = 5.2 #Initial x position of Jupiter, in AU
y_j[0] = 0.0 #Initial y position of Jupiter, in AU
v_jx[0] = 0.0 #Initial x velocity of Jupiter, in AU/yr
v_jy[0] = 2.63 #Initial y velocity of Jupiter, in AU/yr
r_j[0] = (x_j[0]**2+y_j[0]**2)**0.5 #initial value of orbital radius of Jupiter


# Jupiter position & velocity due to Sun
for i in range(len(t)-1):
    v_jx[i+1] = v_jx[i] - G * M_s * x_j[i] / r_j[i]**(3) * dt 
    x_j[i+1] = x_j[i] + v_jx[i+1] * dt
    v_jy[i+1] = v_jy[i] - G * M_s * y_j[i] / r_j[i]**(3) * dt 
    y_j[i+1] = y_j[i] + v_jy[i+1] * dt
    r_j[i+1] = (x_j[i+1]**2 + y_j[i+1]**2)**0.5

v_exf = np.zeros(len(t)) # velocity of Earth in x-direction
v_exf[0] = 0.0 #Initial x velocity of Earth, in AU/yr
v_eyf = np.zeros(len(t)) # velocity of Earth in y-direction
v_eyf[0] = 6.18 #Initial y velocity of Earth, in AU/yr
x_ef = np.zeros(len(t)) # position of Earth in x-direction, in AU
x_ef[0] = 1.0 #Initial x position of Earth, in AU
y_ef = np.zeros(len(t)) # position of Earth in y-direction, in AU
y_ef[0] = 0.0 #Initial y position of Earth, in AU
r_ef = np.zeros(len(t)) # distance between Earth & Sun
r_ef[0] = (x_ef[0]**2+y_ef[0]**2)**0.5 #Initial distance between Earth & Sun
r_ej = np.zeros(len(t)) # distance between Earth & Jupiter
r_ej[0] = ((x_ef[0]-x_j[0])**2 + (y_ef[0]-y_j[0])**2)**0.5 #Initial distance between Earth & Jupiter

# Earth position & velocity due to Sun & Jupiter
for i in range(len(t)-1):
    # Velocity of Earth in x-direction
    v_exf[i+1] = v_exf[i] - G * M_s * x_ef[i] / r_ef[i]**3 * dt - G * M_j * (x_ef[i] - x_j[i]) / (r_ej[i])**3 * dt
    # Position of Earth in x-direction
    x_ef[i+1] = x_ef[i] + v_exf[i+1] * dt
    # Velocity of Earth in y-direction
    v_eyf[i+1] = v_eyf[i] - G * M_s * y_ef[i] / r_ef[i]**3 * dt - G * M_j * (y_ef[i] - y_j[i]) / (r_ej[i])**3 * dt
    # Position of Earth in y-direction
    y_ef[i+1] = y_ef[i] + v_eyf[i+1]*dt
    # distance between Earth & Sun
    r_ef[i+1] = (x_ef[i+1]**2 + y_ef[i+1]**2)**0.5
    # distance between Earth & Jupiter
    r_ej[i+1] = ((x_ef[i+1]-x_j[i+1])**2 + (y_ef[i+1]-y_j[i+1])**2)**0.5

plt.figure()
plt.scatter(0,0, c ='red', label = 'The Sun')
plt.plot(x_j,y_j, label = 'Jupiter')
plt.plot(x_ef,y_ef, label = ' Earth')
plt.title('Orbit of Jupiter and Earth')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.legend(loc='upper right')
plt.axis('scaled')
plt.grid()


# # Q2 (b)

M_j = 2e30 # Mass of Jupiter, in kg
dt = 0.0001 #Time step, in yr
t = np.arange(0, 3 , dt) #Time interval of 3 year
M_s = 2e30 #Mass of the sun, in kg
G = 39.5 / M_s #Gravitational constant, in AU^3 M_s^-1 yr^-2

#np array to store values of x and y position and x and y velocity for Q1 c)
v_jx = np.zeros(len(t)) #Velocity of Jupiter in x-direction. Create np array with same length as t, with all indexes 0
x_j = np.zeros(len(t)) #Position of Jupiter in x-direction. Create np array with same length as t, with all indexes 0
v_jy = np.zeros(len(t)) #Velocity of Jupiter in y-direction. Create np array with same length as t, with all indexes 0
y_j = np.zeros(len(t)) #Position of Jupiter in y-direction. Create np array with same length as t, with all indexes 0
r_j = np.zeros(len(t)) #Radius of Jupiter's orbit. Create np array with same length as t, with all indexes 0

x_j[0] = 5.2 #Initial x position of Jupiter, in AU
y_j[0] = 0.0 #Initial y position of Jupiter, in AU
v_jx[0] = 0.0 #Initial x velocity of Jupiter, in AU/yr
v_jy[0] = 2.63 #Initial y velocity of Jupiter, in AU/yr
r_j[0] = (x_j[0]**2+y_j[0]**2)**0.5 #initial value of orbital radius of Jupiter


# Jupiter position & velocity due to Sun
for i in range(len(t)-1):
    v_jx[i+1] = v_jx[i] - G * M_s * x_j[i] / r_j[i]**(3) * dt 
    x_j[i+1] = x_j[i] + v_jx[i+1] * dt
    v_jy[i+1] = v_jy[i] - G * M_s * y_j[i] / r_j[i]**(3) * dt 
    y_j[i+1] = y_j[i] + v_jy[i+1] * dt
    r_j[i+1] = (x_j[i+1]**2 + y_j[i+1]**2)**0.5


v_exf = np.zeros(len(t)) # velocity of Earth in x-direction
v_exf[0] = 0.0 #Initial x velocity of Earth, in AU/yr
v_eyf = np.zeros(len(t)) # velocity of Earth in y-direction
v_eyf[0] = 6.18 #Initial y velocity of Earth, in AU/yr
x_ef = np.zeros(len(t)) # position of Earth in x-direction, in AU
x_ef[0] = 1.0 #Initial x position of Earth, in AU
y_ef = np.zeros(len(t)) # position of Earth in y-direction, in AU
y_ef[0] = 0.0 #Initial y position of Earth, in AU
r_ef = np.zeros(len(t)) # distance between Earth & Sun
r_ef[0] = (x_ef[0]**2+y_ef[0]**2)**0.5 #Initial distance between Earth & Sun
r_ej = np.zeros(len(t)) # distance between Earth & Jupiter
r_ej[0] = ((x_ef[0]-x_j[0])**2 + (y_ef[0]-y_j[0])**2)**0.5 #Initial distance between Earth & Jupiter

# Earth position & velocity due to Sun & Jupiter
for i in range(len(t)-1):
    # Velocity of Earth in x-direction
    v_exf[i+1] = v_exf[i] - G * M_s * x_ef[i] / r_ef[i]**3 * dt - G * M_j * (x_ef[i] - x_j[i]) / (r_ej[i])**3 * dt
    # Position of Earth in x-direction
    x_ef[i+1] = x_ef[i] + v_exf[i+1] * dt
    # Velocity of Earth in y-direction
    v_eyf[i+1] = v_eyf[i] - G * M_s * y_ef[i] / r_ef[i]**3 * dt - G * M_j * (y_ef[i] - y_j[i]) / (r_ej[i])**3 * dt
    # Position of Earth in y-direction
    y_ef[i+1] = y_ef[i] + v_eyf[i+1]*dt
    # distance between Earth & Sun
    r_ef[i+1] = (x_ef[i+1]**2 + y_ef[i+1]**2)**0.5
    # distance between Earth & Jupiter
    r_ej[i+1] = ((x_ef[i+1]-x_j[i+1])**2 + (y_ef[i+1]-y_j[i+1])**2)**0.5

plt.figure()
plt.scatter(0,0, c ='red', label = 'The Sun')
plt.plot(x_j,y_j, label = 'Jupiter')
plt.plot(x_ef,y_ef, label = ' Earth')
plt.title('Orbit of Jupiter and Earth when Mass of Jupiter gets 1000 times bigger')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.legend(loc='upper right')
plt.axis('scaled')
plt.grid()


# # Q2 (c)

M_j = 10**(-3) * M_s
t = np.arange(0.0, 20.0 , dt) #Time interval of 20 year

#np array to store values of x and y position and x and y velocity for Q1 c)
v_jx = np.zeros(len(t)) #Velocity of Jupiter in x-direction. Create np array with same length as t, with all indexes 0
x_j = np.zeros(len(t)) #Position of Jupiter in x-direction. Create np array with same length as t, with all indexes 0
v_jy = np.zeros(len(t)) #Velocity of Jupiter in y-direction. Create np array with same length as t, with all indexes 0
y_j = np.zeros(len(t)) #Position of Jupiter in y-direction. Create np array with same length as t, with all indexes 0
r_j = np.zeros(len(t)) #Radius of Jupiter's orbit. Create np array with same length as t, with all indexes 0

x_j[0] = 5.2 #Initial x position of Jupiter, in AU
y_j[0] = 0.0 #Initial y position of Jupiter, in AU
v_jx[0] = 0.0 #Initial x velocity of Jupiter, in AU/yr
v_jy[0] = 2.63 #Initial y velocity of Jupiter, in AU/yr
r_j[0] = (x_j[0]**2+y_j[0]**2)**0.5 #initial value of orbital radius of Jupiter


# Jupiter position & velocity due to Sun
for i in range(len(t)-1):
    v_jx[i+1] = v_jx[i] - G * M_s * x_j[i] / r_j[i]**(3) * dt 
    x_j[i+1] = x_j[i] + v_jx[i+1] * dt
    v_jy[i+1] = v_jy[i] - G * M_s * y_j[i] / r_j[i]**(3) * dt 
    y_j[i+1] = y_j[i] + v_jy[i+1]*dt
    r_j[i+1] = (x_j[i+1]**2 + y_j[i+1]**2)**0.5


    
v_exf = np.zeros(len(t)) # velocity of asteriod in x-direction
v_exf[0] = 0.0 #Initial x velocity of asteriod, in AU/yr
v_eyf = np.zeros(len(t)) # velocity of asteriod in y-direction
v_eyf[0] = 3.46 #Initial y velocity of asteriod, in AU/yr
x_ef = np.zeros(len(t)) # position of asteriod in x-direction, in AU
x_ef[0] = 3.3 #Initial x position of asteriod, in AU
y_ef = np.zeros(len(t)) # position of asteriod in y-direction, in AU
y_ef[0] = 0.0 #Initial y position of asteriod, in AU
r_ef = np.zeros(len(t)) # distance between asteriod & Sun
r_ef[0] = (x_ef[0]**2+y_ef[0]**2)**0.5 #Initial distance between asteriod & Sun
r_ej = np.zeros(len(t)) # distance between asteriod & Jupiter
r_ej[0] = ((x_ef[0]-x_j[0])**2 + (y_ef[0]-y_j[0])**2)**0.5 #Initial distance between asteriod & Jupiter


# asteriod position & velocity due to Sun & Jupiter
for i in range(len(t)-1):
    # Velocity of asteriod in x-direction
    v_exf[i+1] = v_exf[i] - G * M_s * x_ef[i] / r_ef[i]**3 * dt - G * M_j * (x_ef[i] - x_j[i]) / (r_ej[i])**3 * dt
    # Position of asteriod in x-direction
    x_ef[i+1] = x_ef[i] + v_exf[i+1] * dt
    # Velocity of asteriod in y-direction
    v_eyf[i+1] = v_eyf[i] - G * M_s * y_ef[i] / r_ef[i]**3 * dt - G * M_j * (y_ef[i] - y_j[i]) / (r_ej[i])**3 * dt
    # Position of asteriod in y-direction
    y_ef[i+1] = y_ef[i] + v_eyf[i+1]*dt
    # distance between asteriod & Sun
    r_ef[i+1] = (x_ef[i+1]**2 + y_ef[i+1]**2)**0.5
    # distance between asteriod & Jupiter
    r_ej[i+1] = ((x_ef[i+1]-x_j[i+1])**2 + (y_ef[i+1]-y_j[i+1])**2)**0.5

plt.figure(figsize =(5,5))
plt.scatter(0,0, c ='red', label = 'The Sun')
plt.plot(x_ef,y_ef, label = ' asteriod')
plt.title('Orbit of asteriod')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.legend(loc='upper right')
plt.axis('scaled')
plt.grid()


