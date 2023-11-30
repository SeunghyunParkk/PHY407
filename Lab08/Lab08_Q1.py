# -*- coding: utf-8 -*-
"""
Lab08 Question 1
Electrostatics and Laplace's equations
Done by Juann Jeon
"""
import numpy as np
import matplotlib.pyplot as plt
import time

grid = 100 #Size of the n * n matrix
V = 1 #Voltage
precision = 1e-6 #Desired precision

#Matrix for electronic capacitor in n * n matrix
phi_1 = np.zeros([grid, grid], float) 
phi_2 = np.zeros([grid, grid], float) 
phi_3 = np.zeros([grid, grid], float) 
z_axis = np.linspace(-V, V, num = 41) #Z axis for voltage

#Overrelaxation constant
omega_2 = 0.1 
omega_3 = 0.5

#Q 1a)
t1 = time.time() #Start recording time it takes
error = precision * 2
while error > precision:
    error = 0
    #Row of the matrix
    for i in range(1, grid - 1):
        #Column of the matrix
        for j in range(1, grid - 1):
            #Positive plate with negligible thickness
            if i > 19 and i < 79 and j == 19:
                phi_1[i, j] = V
            #Negative plate with negligible thickness
            elif i > 19 and i < 79 and j == 79:
                phi_1[i, j] = -V
            #Gauss-Seidel method without over-relaxation
            else: 
                gauss_seidel = (1/4) * (phi_1[i + 1, j] + phi_1[i - 1, j] + phi_1[i, j + 1] + phi_1[i, j - 1])
                
                #Check the error
                if gauss_seidel - phi_1[i, j] > error:
                    error = gauss_seidel - phi_1[i, j]
                
                phi_1[i, j] = gauss_seidel

t2 = time.time() #End recording time it takes
print("Total time it took for Q 1a) (in s):", t2 - t1)

plt.figure()
plt.contourf(phi_1, z_axis)
plt.xlabel('x (mm)', fontsize = '13')
plt.ylabel('y (mm)', fontsize = '13')
plt.title('Contour plot of the potential', fontsize = '13')
plt.colorbar(label = "Potential (V)")
plt.tight_layout()


#Q 1b)
t3 = time.time() #Start recording time it takes
error = precision * 2
while error > precision:
    error = 0
    #Row of the matrix
    for i in range(1, grid - 1):
        #Column of the matrix
        for j in range(1, grid - 1):
            #Positive plate with negligible thickness
            if i > 19 and i < 79 and j == 19:
                phi_2[i, j] = V
            #Negative plate with negligible thickness
            elif i > 19 and i < 79 and j == 79:
                phi_2[i, j] = -V
            #Gauss-Seidel method overrelaxation ω = 0.1
            else: 
                gauss_seidel_over = ((1 + omega_2)/4) * (phi_2[i + 1, j] + phi_2[i - 1, j] + phi_2[i, j + 1] + phi_2[i, j - 1]) - (omega_2 * phi_2[i, j])
                
                #Check the error
                if gauss_seidel_over - phi_2[i, j] > error:
                    error = gauss_seidel_over - phi_2[i, j]
                
                phi_2[i, j] = gauss_seidel_over

t4 = time.time() #End recording time it takes
print("Total time it took for omega = 0.1 (in s):", t4 - t3)

plt.figure()
plt.tight_layout()
plt.contourf(phi_2, z_axis)
plt.xlabel('x (mm)', fontsize = '13')
plt.ylabel('y (mm)', fontsize = '13')
plt.title('Contour plot of the potential using over-relaxation, ω = 0.1', fontsize = '13')
plt.colorbar(label = "Potential (V)")
plt.tight_layout()


error = precision * 2
t5 = time.time() #Start recording time it takes
while error > precision:
    error = 0
    #Row of the matrix
    for i in range(1, grid - 1):
        #Column of the matrix
        for j in range(1, grid - 1):
            #Positive plate with negligible thickness
            if i > 19 and i < 79 and j == 19:
                phi_3[i, j] = V
            #Negative plate with negligible thickness
            elif i > 19 and i < 79 and j == 79:
                phi_3[i, j] = -V
            #Gauss-Seidel method overrelaxation ω = 0.5
            else: 
                gauss_seidel_over = ((1 + omega_3)/4) * (phi_3[i + 1, j] + phi_3[i - 1, j] + phi_3[i, j + 1] + phi_3[i, j - 1]) - (omega_3 * phi_3[i, j])
                
                if gauss_seidel_over - phi_3[i, j] > error:
                    error = gauss_seidel_over - phi_3[i, j]
                
                phi_3[i, j] = gauss_seidel_over

t6 = time.time() #End recording time it takes
print("Total time it took for omega = 0.5 (in s):", t6 - t5)

plt.figure()
plt.contourf(phi_3, z_axis)
plt.xlabel('x (mm)', fontsize = '13')
plt.ylabel('y (mm)', fontsize = '13')
plt.title('Contour plot of the potential using over-relaxation, ω = 0.5', fontsize = '13')
plt.colorbar(label = "Potential (V)")
plt.tight_layout()

