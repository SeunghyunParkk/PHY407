
"""
Lab09_Q2
Resonant EM cavity
reconstruct electric and magnetic fields using Crank-Nicolson Method

Done by Seunghyun Park
"""
import dcst
import numpy as np
import matplotlib.pyplot as plt


# Q2 (a)

def dsct2(f):
    M = f.shape[0] # Number of rows
    N = f.shape[1] # Number of columns
    a = np.zeros((M, N)) # Intermediate array
    b = np.zeros((M, N)) # Final array
    # Take transform along x
    for j in range(N):
        a[:,j] = dcst.dct(f[:,j])
    for i in range(M):
        b[i,:] = dcst.dst(a[i,:])
    return b

def idsct2(b):
    M = b.shape[0] # Number of rows
    N = b.shape[1] # Number of columns
    a = np.zeros((M, N)) # Intermediate array
    f = np.zeros((M, N)) # Final array
    for i in range(M):
        a[i,:] = dcst.idst(b[i,:])

    for j in range(N):
        f[:,j] = dcst.idct(a[:,j])

    return f

def dcst2(f):
    M = f.shape[0] # Number of rows
    N = f.shape[1] # Number of columns
    a = np.zeros((M, N)) # Intermediate array
    b = np.zeros((M, N)) # Final arra
    # Take transform along x
    for j in range(N):
        a[:,j] = dcst.dst(f[:,j])
    for i in range(M):
        b[i,:] = dcst.dct(a[i,:])
    return b

def idcst2(b):
    M = b.shape[0] # Number of rows
    N = b.shape[1] # Number of column
    a = np.zeros((M, N)) # Intermediate array
    f = np.zeros((M, N)) # Final array
    for i in range(M):
        a[i,:] = dcst.idct(b[i,:])
    for j in range(N):
        f[:,j] = dcst.idst(a[:,j])
    return f

# create an array for testing the functions
test = np.array([[0.,0.,0.],[0.,1.,1.],[0.,1.,1.]])
result = idsct2(dsct2(test))
print("The original 2D array is")
print(test)
print()
print("The inverse fourier transform of fourier transform of f is")
print(result)

test = np.array([[0.,0.,0.],[0.,1.,1.],[0.,1.,1.]])
result = idcst2(dcst2(test))
print()
print("The inverse fourier transform of fourier transform of f is")
print(result)

# Q2 (b)

tau = 0.01 # time step
N = 2000 # number of time step
time = np.arange(0, N * tau, tau)
Lx = 1 # Lx
Ly = 1 # Ly
J0 = 1 # value of J0
m = 1 # value of m
n = 1 # value of n
c = 1 # value of c
P = 32 # number of grid
w = 3.75 # angular frequency
ax = Lx/P # step size of x component
ay = Ly/P # step size of y component

Dx = np.pi * c * tau / (2*Lx)
Dy = np.pi * c * tau / (2*Ly)


x = np.arange(0, Lx, ax) # x
y = np.arange(0, Ly, ay) # y



# define function J_z
def J_z(x, y, t):
    jz = np.empty((len(x), len(y)))
    for i in range(len(x)):
        for j in range((len(y))):
            jz[i, j] = J0 * np.sin(m * np.pi * x[i] / Lx) * np.sin(n * np.pi * y[j] / Ly) * np.sin(w * t)
    return jz

# initial values
Ez = np.zeros((P+1, P+1))
Hx = np.zeros((P+1, P+1))
Hy = np.zeros((P+1, P+1))
Jz = J_z(x, y, time[0])

# fourier transfrom Ez,Hx,Hy,Jz
f_Ez = dcst.dst2(Ez)
f_Hx = dcst2(Hx)
f_Hy = dsct2(Hy)
f_Jz = dcst.dst2(Jz)


#create arrays
E_new = np.zeros((P+1, P+1), float)
X_new = np.zeros((P+1, P+1), float)
Y_new = np.zeros((P+1, P+1), float)

# create arrays for results
Ezfinal = []
Hxfinal = []
Hyfinal = []

# calculate Ez,Hx,Hy using for loop
# k is the index of time
for k in range(len(time)) :
    Jz = J_z(x, y, time[k])
    f_Jz = dcst.dst2(Jz)
    # calculathe the Ez,Hx,Hy using crank-nicolson method
    for i in range(1, P):
        for j in range(1, P):
            p = i
            q = j
            E_new[i, j] = ((1 - p**2 * Dx**2 - q**2 * Dy**2) * f_Ez[i, j] + 2 * q * Dy * f_Hx[i, j] + 2 * p * Dx * f_Hy[i, j] + tau * f_Jz[i, j]) / (1 + p**2 * Dx**2 + q**2 * Dy**2)
            X_new[i, j] = f_Hx[i, j] - q * Dy * (E_new[i, j] + f_Ez[i, j])
            Y_new[i, j] = f_Hy[i ,j] - p * Dx * (E_new[i, j] + f_Ez[i, j])
            f_Ez[i, j] = np.copy(E_new[i, j])
            f_Hx[i, j] = np.copy(X_new[i, j])
            f_Hy[i, j] = np.copy(Y_new[i, j])
    # apply inverse fourier transform to reconstruct Ez, Hx, Hy
    Ezfinal.append(dcst.idst2(f_Ez)[16, 16])
    Hxfinal.append(idcst2(f_Hx)[16,0])
    Hyfinal.append(idsct2(f_Hy)[0,16])

    
t = np.arange(0, N * tau, tau)

# plot Ez
plt.figure()
plt.plot(t,Ezfinal)
plt.xlabel('time (s)', fontsize = '13')
plt.ylabel('E_z (arbitrary unit)', fontsize = '13')
plt.title('Trace of E_z(x=0.5,y=0.5)', fontsize = '13')
# plot Hx
plt.figure()
plt.plot(t,Hxfinal)
plt.xlabel('time (s)', fontsize = '13')
plt.ylabel('H_x (arbitrary unit)', fontsize = '13')
plt.title('Trace of H_x(x=0.5,y=0.0)', fontsize = '13')
# plot Hy
plt.figure()
plt.plot(t,Hyfinal)
plt.xlabel('time (s)', fontsize = '13')
plt.ylabel('H_y (arbitrary unit)', fontsize = '13')
plt.title('Trace of H_y(x=0.0,y=0.5)', fontsize = '13')



