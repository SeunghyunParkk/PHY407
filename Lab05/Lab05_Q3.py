# -*- coding: utf-8 -*-
"""
Lab05_Q3
Analysis of sea level pressure

Done by Juann Jeon
"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.pyplot import contourf, xlabel, ylabel, title, colorbar

#load text 
SLP = np.loadtxt('SLP.txt')
Longitude = np.loadtxt('lon.txt')
Times = np.loadtxt('times.txt')

#apply 2D fourier transform on SLP data
fft2D = np.fft.rfft2(SLP)

#assign zero matrix 
wave3 = np.zeros((len(fft2D), len(fft2D[0])), dtype = 'complex')
wave4 = np.zeros((len(fft2D), len(fft2D[0])), dtype = 'complex')
wave5 = np.zeros((len(fft2D), len(fft2D[0])), dtype = 'complex')


for i in range(len(fft2D)):
    wave3[i][3] = fft2D[i][3]
    wave4[i][4] = fft2D[i][4]
    wave5[i][5] = fft2D[i][5]
    
# apply inverse 2D fourier transform
ifft2D3 = np.fft.irfft2(wave3)
ifft2D4 = np.fft.irfft2(wave4)
ifft2D5 = np.fft.irfft2(wave5)

#plot the inversed 2D fourier transform
plt.figure()
plt.contourf(Longitude, Times, SLP)
plt.xlabel('longitude(degrees)',fontsize = "13")
plt.ylabel('days since Jan. 1 2015',fontsize = "13")
plt.title('SLP anomaly (hPa)',fontsize = "13")
plt.colorbar()

plt.figure()
plt.contourf(Longitude, Times, ifft2D3)
plt.xlabel('longitude(degrees)',fontsize = "13")
plt.ylabel('days since Jan. 1 2015',fontsize = "13")
plt.title('SLP anomaly (hPa) of wavenumber 3',fontsize = "13")
plt.colorbar()

plt.figure()
plt.contourf(Longitude, Times, ifft2D4)
plt.xlabel('longitude(degrees)',fontsize = "13")
plt.ylabel('days since Jan. 1 2015',fontsize = "13")
plt.title('SLP anomaly (hPa) of wavenumber 4',fontsize = "13")
plt.colorbar()

plt.figure()
plt.contourf(Longitude, Times, ifft2D5)
plt.xlabel('longitude(degrees)',fontsize = "13")
plt.ylabel('days since Jan. 1 2015',fontsize = "13")
plt.title('SLP anomaly (hPa) of wavenumber 5',fontsize = "13")
plt.colorbar()
