# -*- coding: utf-8 -*-
"""
@author: Seunghyun Park
Lab02_Q1 

"""
import numpy as np

# Q1 (b)
# import the data
data = np.loadtxt("cdata.txt")
# Calculating correct answer by using numpy.std
data_std = np.std(data, ddof =1)

# define standard devation function (1)
def std_1(x):
    n = len(x) # length of data
    x_sum = 0 # Assign initial value of sum of data
    c = 0. # Assign initial value of (x[i]-x_mean)**2
    std = 0. # Assign initial value of standard deviation
    for i in range(0,n):
        x_sum += x[i] # Add up all components of data
    x_mean = x_sum/n # Calculate the mean value of the data
    for i in range(0,n):
        c += (x[i]-x_mean)**2 # Calculate (x[i]-x_mean)**2
    std = np.sqrt(1./(n-1.)*c)
    return std

def std_2(x):
    n = len(x)# length of data
    x_sum =0. # Assign initial value of sum of data
    c = 0. # Assign initial value of (x[i])**2
    std = 0.# Assign initial value of standard deviation
    for i in range(0,n):
        c += x[i]**2
        x_sum += x[i] # Add up all components of data
    x_mean = x_sum / n # Calculate the mean value of the data
    if c < n*x_mean**2: # print warning if it takes the square root of a negative number
        print("Warning, this calculation contains the square root of a negative number")
    else:
        std = np.sqrt(1./(n-1.)*(c - n*x_mean**2))
    return std


print('The standard deviation of speed of light using numpy.std is', data_std)
print('The standard deviation of speed of light using equation 1 is', std_1(data))
print('The standard deviation of speed of light using equation 2 is', std_2(data))

#The relative error of the calculated standard deviation using equation 1
re_error1 = np.abs((std_1(data) - data_std)/ data_std)
#The relative error of the calculated standard deviation using equation 2
re_error2 = np.abs((std_2(data) - data_std)/ data_std)

print('The relative error of the calculated standard deviation using equation 1: ', re_error1)
print('The relative error of the calculated standard deviation using equation 2: ', re_error2)
# The relative error of second equation is larger in magnitude


# Q1 (c)

# Set seed
np.random.seed(815)
# Normal distribution with mean = 0, standard deviation = 1, and length of 2000
norm_dist1 = np.random.normal(0., 1.0, 2000)
# Normal distribution with mean = 1.e7, standard deviation = 1, and length of 2000
norm_dist2 = np.random.normal(1.e7, 1.0, 2000)

# Calculate the correct standard deviation using numpy.std
norm_dist1_std = np.std(norm_dist1,ddof=1)
norm_dist2_std = np.std(norm_dist2,ddof=1)

print('The standard deviation of normal distribution with mean = 0 using numpy.std is', norm_dist1_std)
print('The standard deviation of normal distribution with mean = 0 using equation 1 is', std_1(norm_dist1))
print('The standard deviation of normal distribution with mean = 0 using equation 2 is', std_2(norm_dist1))

print('The standard deviation of normal distribution with mean = 1.e7 using numpy.std is', norm_dist2_std)
print('The standard deviation of normal distribution mean = 1.e7 using equation 1 is', std_1(norm_dist2))
print('The standard deviation of normal distribution mean = 1.e7 using equation 2 is', std_2(norm_dist2))


#The relative error of the calculated standard deviation of normal distribution with mean = 0 using equation 1
re_error3 = np.abs((std_1(norm_dist1) - norm_dist1_std)/ norm_dist1_std)
#The relative error of the calculated standard deviation of normal distribution with mean = 0 using equation 2
re_error4 = np.abs((std_2(norm_dist1) - norm_dist1_std)/ norm_dist1_std)

#The relative error of the calculated standard deviation of normal distribution with mean = 1.e7 using equation 1
re_error5 = np.abs((std_1(norm_dist2) - norm_dist2_std)/ norm_dist2_std)
#The relative error of the calculated standard deviation of normal distribution with mean = 1.e7 using equation 2
re_error6 = np.abs((std_2(norm_dist2) - norm_dist2_std)/ norm_dist2_std)

print('The relative error of std of normal with mean 0 using equation 1: ', re_error3)
print('The relative error of std of normal with mean 0 using equation 2: ', re_error4)
print('The relative error of std of normal with mean 1e7 using equation 1: ', re_error5)
print('The relative error of std of normal with mean 1e7 using equation 2: ', re_error6)


# Q1 (d)

# def fix_std_2(x):
    
def fix_std_2(x):
    n = len(x)# length of data
    x_sum =0.# Assign initial value of sum of data
    c = 0.# Assign initial value of sum of x[i]**2
    std = 0.# Assign initial value of standard deviation
    x = x-np.mean(x) # Make x close to 0 so it becomes more accurate
    for i in range(0,n):
        c += x[i]**2
        x_sum += x[i] # Add up all components of data
    x_mean = x_sum / n # Calculate the mean value of the data
    if c < n*x_mean**2: # print warning if it takes the square root of a negative number
        print("Warning, this calculation contains the square root of a negative number")
    else:
        std = np.sqrt(1./(n-1.)*(c - n*x_mean**2))
    return std

re_error7 = np.abs((fix_std_2(norm_dist2) - norm_dist2_std)/ norm_dist2_std)
print("The standard deviation of normal distribution with mean 1.e7 using fixed equation 2 is",fix_std_2(norm_dist2))
print('The relative error of std of normal with mean 1e7 using fixed equation 2: ', re_error7)







