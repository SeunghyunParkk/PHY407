# Question 3
# Done by Seunghyun Park 

import numpy as np
import time
import matplotlib.pyplot as plt


# Define matrix function to generate Matrix
def matrix(s):
    M = np.ones([s,s], float)*3
    return M

# Assign N, 2 to a few hundred
N = [2,10,50,100,200,300]

t_f = np.zeros(len(N)) # time takes using forloop in the textbook 
t_d = np.zeros(len(N)) # time takes using numpy dot

for i in range(len(N)):
    # Generate matrices A & B
    A = matrix(N[i])*3
    B = matrix(N[i])*3
    
    #Create a N x N matrix with all indexes 0
    C = np.zeros([N[i],N[i]],float)

    t1 = time.time() # initial time (for loop)
    for a in range(N[i]): 
        for b in range(N[i]):
            for c in range(N[i]):
                C[a,b] += A[a,c]*B[c,b]
    t2 = time.time() # final time (for loop)
    
    t3 = time.time() # initial time (numpy dot) 
    C = np.dot(A,B)
    t4 = time.time() # final time (numpy dot)
    t_f[i] = (t2-t1) 
    t_d[i] = (t4-t3)
    
for i in range(len(N)):
    print('The time taken by for loop matrix multiplication when N = '+str(N[i])+ ' is ' + str(t_f[i]) +" seconds")
    print('The time taken by numpy dot when N = '+ str(N[i]) + ' is ' + str(t_d[i]) +" seconds")
      
plt.figure()
plt.plot(N, t_f, 'r', label = 'For loop')
plt.plot(N, t_d, 'b', label = 'Numpy Dot')
plt.legend(loc='upper left')   
plt.xlabel('N')
plt.ylabel('time (s)')
plt.title("Speed test ( t vs N)")

N3 = np.power(N,3) # N^3 using numpy power function

plt.figure()
plt.plot(N3, t_f, 'r', label = 'For loop')
plt.plot(N3, t_d, 'b', label = 'Numpy Dot')
plt.legend(loc='upper left')   
plt.xlabel('$N^{3}$')
plt.ylabel('time (s)')
plt.title("Speed test ( t vs $N^{3}$)")

