# -*- coding: utf-8 -*-
"""
Lab05_Q1
Revisiting the relativistic spring

Done by Seunghyun Park
"""

from scipy import constants
import matplotlib.pyplot as plt
import numpy as np

m = 1 #Mass, in kg
k = 12 #Spring constant, in N/m
c = constants.speed_of_light #Speed of light, in m/s
#xc from Lab03 2c), approximately 8.65107m
xc = c*(m/k)**0.5

#define a function that calculates x and dx/dt using Euler-Cromer method
def relativistic(x0,t,dt):
    #Assign zero list with length of t for u1, u2, du1, du2
    u1 = np.zeros(len(t))
    u2 = np.zeros(len(t))
    du1 = np.zeros(len(t))
    du2 = np.zeros(len(t))
    #Assign initial value of u1
    u1[0] = x0
    #Calculate the position and velocity of relativistic spring using Euler-Cromer method
    for i in range(len(t)-1):
        du2[i+1] = -k/m*u1[i]*(1-u2[i]**2/c**2)**(3/2)
        u2[i+1] = u2[i]+du2[i+1]*dt
        du1[i+1] = u2[i+1]
        u1[i+1] = u1[i] + du1[i+1]*dt 
    return u1,u2

# Assign time
t = np.linspace(0,130,2**20)
dt = t[1]-t[0]

#Calculate the position and velocity of relavistic spring with differenct x0
#and plot them
u11 ,u12 = relativistic(1,t,dt)
plt.figure(figsize = (10,2))
plt.plot(t,u11)
plt.title("Position of relativistic spring, x0 = 1",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("Position (m)",fontsize = "13")
plt.figure(figsize = (10,2))
plt.plot(t,u12)
plt.title("Velocity of relativistic spring, x0 = 1",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("velocity (m/s)",fontsize = "13")
    

u21 ,u22 = relativistic(xc,t,dt)
plt.figure(figsize = (10,2))
plt.plot(t,u21)
plt.title("Position of relativistic spring, x0 = xc",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("Position (m)",fontsize = "13")
plt.figure(figsize = (10,2))
plt.plot(t,u22)
plt.title("Velocity of relativistic spring, x0 = xc",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("velocity (m/s)",fontsize = "13")



u31 ,u32 = relativistic(10*xc,t,dt)
plt.figure(figsize = (10,2))
plt.plot(t,u31)
plt.title("Position of relativistic spring, x0 = 10xc",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("Position (m)",fontsize = "13")
plt.figure(figsize = (10,2))
plt.plot(t,u32)
plt.title("Velocity of relativistic spring, x0 = 10xc",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("velocity (m/s)",fontsize = "13")

#(b)
#Apply fourier transform on position
fftx1 = np.fft.rfft(u11)
fftx2 = np.fft.rfft(u21)
fftx3 = np.fft.rfft(u31)
#find the corresponding angular frequency
fre = np.fft.rfftfreq(len(u11),dt)

#plot fft
plt.figure()
plt.plot(fre,abs(fftx1/max(fftx1)),label = 'x0=1')
plt.plot(fre,abs(fftx2/max(fftx2)), label = 'x0=xc')
plt.plot(fre,abs(fftx3/max(fftx3)),label = 'x0=10xc')
plt.title("Scaled fourier component",fontsize = "13")
plt.legend()
plt.xlim(0,1)
plt.xlabel("Freqeuncy (Hz)",fontsize = "13")
plt.ylabel("$|\hat{x}(\omega)|/|\hat{x}(\omega)|_{max}$",fontsize = "13")


#(c)
#define guassian quadrature integration
def gaussxw(N): 
    # Initial approximation to roots of the Legendre polynomial 
    a = np.linspace(3, 4*N-1, N)/(4*N+2) 
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))  
    # Find roots using Newton's method 
    epsilon = 1e-15 
    delta = 1.0 
    while delta>epsilon: 
        p0 = np.ones(N, float) 
        p1 = np.copy(x) 
        for k in range(1, N): 
            p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1) 
        dp = (N+1)*(p0-x*p1)/(1-x*x) 
        dx = p1/dp 
        x -= dx 
        delta= max(abs(dx)) 
    # Calculate the weights 
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp) 
    return x,w 

#define integrand from lab03
def g(x, x_0):
    fx = (x_0**2) - (x**2)
    f1 = k * fx
    f2 = (2 * m * c**2) + k * fx/2
    f3 = 2 * ((m * c**2) + k * fx/2)**2
    return 4*(c * ((f1 * f2)/f3)**0.5)**-1

m = 1 #Mass, in kg
k = 12 #Spring constant, in N/m
c = constants.speed_of_light #Speed of light, in m/s

xc = c*(m/k)**0.5


# assign lower bound, upper bound of integral and number of sample, for x0=1m
a1 = 0 #Lower boundary for integral, always 0 in equation (7)
b1 = 1 #b = x_0, in m
N1 = 200

#calculate sample point and weight using guassxw
x_1, w_1 = gaussxw(N1)
#map the sample point and weight to required integration domain
xp_1 = 0.5 * (b1 - a1) * x_1 + 0.5 * (b1 + a1)
wp_1 = 0.5 * (b1 - a1) * w_1
s_1 = 0.0

# calculate the integral
for i in range(N1):
    s_1 += wp_1[i] * g(xp_1[i],b1)    
w1 = 1/s_1

# assign lower bound, upper bound of integral and number of sample, for x0=xc
a2 = 0 #Lower boundary for integral, always 0 in equation (7)
b2 = xc #b = x_0, in m
N2 = 200

#calculate sample point and weight using guassxw
x_2, w_2 = gaussxw(N2)

#map the sample point and weight to required integration domain
xp_2 = 0.5 * (b2 - a2) * x_2 + 0.5 * (b2 + a2)
wp_2 = 0.5 * (b2 - a2) * w_2

s_2 = 0.0
# calculate the integral
for i in range(N2):
    s_2 += wp_2[i] * g(xp_2[i],b2)    
w2 =1/s_2


# assign lower bound, upper bound of integral and number of sample, for x0=10xc
a3 = 0 #Lower boundary for integral, always 0 in equation (7)
b3 = 10*xc #b = x_0, in m
N3 = 200
x_3, w_3 = gaussxw(N3)
xp_3 = 0.5 * (b3 - a3) * x_3 + 0.5 * (b3 + a3)
wp_3 = 0.5 * (b3 - a3) * w_3
s_3 = 0.0

# calculate the integral
for i in range(N3):
    s_3 += wp_3[i] * g(xp_3[i],b3)    
w3 = 1/s_3


#plot fft and predicted period obtained by gaussian quadrature integration
plt.figure()
plt.plot(fre,abs(fftx1/max(fftx1)),label = 'x0=1')
plt.plot(fre,abs(fftx2/max(fftx2)), label = 'x0=xc')
plt.plot(fre,abs(fftx3/max(fftx3)),label = 'x0=10xc')
plt.title("Scaled fourier component and frequencies by gaussian quadrature",fontsize = "13")
plt.axvline(w1,color ='b',label ="x0=1")
plt.axvline(w2,color ='darkorange',label = "x0=xc")
plt.axvline(w3,color = 'g',label = "x0=10xc")
plt.xlabel("Freqeuncy (Hz)",fontsize = "13")
plt.ylabel("$|\hat{x}(\omega)|/|\hat{x}(\omega)|_{max}$",fontsize = "13")
plt.xlim(0,1)
plt.legend()
