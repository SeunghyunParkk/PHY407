# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 05:11:19 2022

@author: Park
"""

import numpy as np
import matplotlib.pyplot as plt


h = 0.01
t = np.linspace(0,100*h,h)

def f(r):
    return 4*d*((sigma/r)**12-(sigma/r)**6)
v = np.zeros(len(t))
x = np.zeros(len(t))
y = np.zeros(len(t))
v[0] = 0

def verlet(r,t,dt,f):
    v[i+1] = v[0] + h/2*f(r[0])
    
    
    
def verlet(v0,r0,h,N,fn):
	"""
	v0 - initial particle velocities
	r0 - initial particle positions
	h - timestep size
	N - number of timesteps
	fn - a function that returns the derivative of each velocity component
	returns positions and velocities for each particle at each timestep
	"""
	# Calculate half step velocity with Euler
	vhalf = v0 + 0.5*h*fn(r0)
	# Create empty arrays to hold resulting positions and velocities
	vs = []
	rs = []
	# Start iteration counter at zero
	i = 0
	while i < N:
		# Add position to the array
		rs.append(np.copy(r0))
		# Update position estimate
		r0 += h*vhalf
		k = h*fn(r0)
		# Update velocity estimate and add it to the array
		v0 = vhalf + 0.5*k
		vs.append(np.copy(v0))
		# Update half step velocity
		vhalf += k
		i+=1
	return np.array(rs),np.array(vs)