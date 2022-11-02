# -*- coding: utf-8 -*-
"""
Lab05_Q2
Audio filtering

Done by Seunghyun Park
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read, write

# import the wav file
sample, data = read('GraviteaTime.wav')

channel_0 = data[:, 0]
channel_1 = data[:, 1]

#(b)
# N is the length of data
N = len(channel_0)

# assign dt
dt= 1/float(sample)
# assign time
t = np.arange(0,N/sample,dt)




# plot original time series of channel 0 & channel 1
plt.figure()
plt.subplot(2,1,1)
plt.plot(t,channel_0)
plt.title("Channel 0 vs time",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")
plt.subplot(2,1,2)
plt.plot(t,channel_1)
plt.title("Channel 1 vs time",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")
plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.5, 
                    hspace=0.7)


#(c)
# find index where t > 0.02 & t < 0.05
index1 = np.where(np.logical_and(0.02< t, t < 0.05))
print("The index of t in the range of 20-50ms is", index1)

# plot the original time series where t is 0.02 < t < 0.05 seconds
plt.figure()
plt.subplot(2,1,1)
plt.plot(t[882:2204],channel_0[882:2204])
plt.title("Channel 0 vs time",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")
plt.subplot(2,1,2)
plt.plot(t[882:2204],channel_1[882:2204])
plt.title("Channel 1 vs time",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")
plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.5, 
                    hspace=0.7)

#(d)
# define filter that coefficient with frequencies greater than 880 Hz make zero.
def slp(x,y):
    
    fc = 880

    for i in range(0,len(y)):
        if abs(y[i]) > fc:
            x[i] = 0
    return x

#apply fourier transform on data
ft0 = (np.fft.rfft(channel_0))
ft1 = (np.fft.rfft(channel_1))
fre = np.fft.rfftfreq(N,dt)

#plot fourier transform 
plt.figure()
plt.subplot(2,1,1)
plt.plot(fre, abs(ft0))
plt.title(" Amplitude of original coefficient of channel 0",fontsize = "13")
plt.xlabel("Frequency (Hz)",fontsize = "13")
plt.ylabel("Amplitude",fontsize = "13")
plt.subplot(2,1,2)
plt.plot(fre, abs(ft1))
plt.title(" Amplitude of original coefficient of channel 1",fontsize = "13")
plt.xlabel("Frequency (Hz)",fontsize = "13")
plt.ylabel("Amplitude",fontsize = "13")
plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.5, 
                    hspace=0.7)

#filter the fourier transform
fil_ft0 = slp(ft0,fre)
fil_c0 = np.fft.irfft(fil_ft0)
fil_ft1 = slp(ft1,fre)
fil_c1 = np.fft.irfft(fil_ft1)

#plot the filtered fourier transform
plt.figure()
plt.subplot(2,1,1)
plt.plot(fre, abs(fil_ft0) )
plt.title(" Amplitude of filtered coefficient channel 0",fontsize = "13")
plt.xlabel("Frequency (Hz)",fontsize = "13")
plt.ylabel("Amplitude",fontsize = "13")
plt.subplot(2,1,2)
plt.plot(fre,abs(fil_ft1) )
plt.title(" Amplitude of filtered coefficient of channel 1",fontsize = "13")
plt.xlabel("Frequency (Hz)",fontsize = "13")
plt.ylabel("Amplitude",fontsize = "13")
plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.5, 
                    hspace=0.7)


#plot the original data
plt.figure()
plt.subplot(2,1,1)
plt.plot(t[882:2204],channel_0[882:2204])
plt.title("Original time series of channel 0",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")

plt.subplot(2,1,2)
plt.plot(t[882:2204],channel_1[882:2204])
plt.title("Original time series of channel 1",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")
plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.5, 
                    hspace=0.7)

#plot the filtered time series
plt.figure()
plt.subplot(2,1,1)
plt.plot(t[882:2204],fil_c0[882:2204])
plt.title("filtered time series of channel 0",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")
plt.subplot(2,1,2)
plt.plot(t[882:2204],fil_c1[882:2204])
plt.title("filtered time series of channel 1",fontsize = "13")
plt.xlabel("time (s)",fontsize = "13")
plt.ylabel("data (arbitrary units)",fontsize = "13")
plt.subplots_adjust(left=0.125,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.5, 
                    hspace=0.7)

# save the filtered data as wav file
data_out = np.empty(data.shape, dtype = data.dtype)
data_out[:, 0] = np.real(fil_c0)
data_out[:, 1] = np.real(fil_c1)
write('GraviteaTime_lpf.wav', sample, data_out)