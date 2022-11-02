'''
Lab06_Q2
Done by Seunghyun Park
Molecular Dynamics simulation
'''
import numpy as np
import matplotlib.pyplot as plt

#part (a)
# define verlet algorithm
def verlet(r,h,t):
    #create empty list for position and velocity
    r_value = []
    v_value = []
    #assign initial velocity
    v0 = 0.
    #calculate v(t+h/2)
    v_h2 = v0+ h/2*acceleration(r)
    #assign initial position
    r0 = r
    #calculate position and velocity using for loop
    for i in range(len(t)):
        rval = r0 + h*v_h2
        k = h*acceleration(rval)
        vval = v_h2 + 0.5*k
        v_h2 = v_h2 + k
        r0 = rval
        r_value.append(rval)
        v_value.append(vval)
    #return position and velocity
    return np.array([r_value, v_value],float)

# define a function that calculates acceleration
def acceleration(r):
    #create empty list for acceleration
    acc = []
    i = 0
    #particle i, total 16 particles
    while i < int(len(r)):
        a_x = 0
        a_y = 0
        j = 0
        # j is the particles that affect particle i
        while j < int(len(r)):
            # if i == j then it is same particle, therefore pass
            if i == j:
                pass
            # calculate acceleration of particle i due to particle j
            else:
                x_d = r[:,0][i] - r[:,0][j]
                y_d = r[:,1][i] - r[:,1][j]
                #sum up all the accelerations on particle i
                a_x += 48*(x_d)/(x_d**2+y_d**2)**7 - 24*(x_d)/(x_d**2+y_d**2)**4
                a_y += 48*(y_d)/(x_d**2+y_d**2)**7 - 24*(y_d)/(x_d**2+y_d**2)**4
            j += 1
        i += 1
        acc.append([a_x,a_y])
    #return acceleration
    return np.array(acc,float)

#this is the code snippet from lab instruction
N = 16
Lx = 4.0
Ly = 4.0
dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)
x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)
xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()

#time  step
dt = 0.01
#time
T = np.arange(0,1000*dt,dt)  

#make a array of positions of particles
r = []
for i in range(len(x_initial)):
    r.append([x_initial[i],y_initial[i]])
r = np.array(r)

#apply verlet algorithm 
r11 ,v11 = verlet(r,dt,T)

#make a list of color
c = ["black","gray","rosybrown","brown","red","salmon","peru","darkorange","yellow","olive","limegreen","tab:olive","violet","mediumaquamarine","steelblue","royalblue"]

#plot the trajectories of 16 particles
plt.figure()
for i in range(len(x_initial)):
    plt.plot(r11[:,i][:,0],r11[:,i][:,1],'.', color = c[i])
plt.title("Tragectories of particles",fontsize= "13")
plt.xlabel("x (arbitrary unit)",fontsize= "13")
plt.ylabel("y (arbitrary unit)",fontsize= "13")

#part (b)
#define a function that calculates kinetic energy of the system
def kinetic(v,T):
    ke = []
    for i in range(len(T)):
        k_e = 0
        for j in range(16):
            # KE = 1/2*m*v^2
            k_e += 0.5*(v11[:,j][:,0][i] + v11[:,j][:,1][i])**2
        ke.append(k_e)
    return ke

#define a function that calculates potential energy of the system
def potential(r,T):
    pe = []
    for i in range(len(T)):
        p_e = 0
        for j in range(16):
            for k in range(16):
                if j !=k:
                    #lennard-jones potential 
                    x_d = (r11[:,j][:,0][i] - r11[:,k][:,0][i])
                    y_d = (r11[:,j][:,1][i] - r11[:,k][:,1][i])
                    r_tot = (x_d**2+y_d**2)**0.5
                    p_e += 4/(r_tot)**12 - 4/r_tot**6
        pe.append(p_e)
    return pe

# calculate kinetic and potential energy of the system
ke = kinetic(v11,T)
pe = potential(r11,T)
# calculate the energy
energy = []
for i in range(len(T)):
    # potential energy of a paricles is V(r)/2 therefore divide by 2
    energy.append(ke[i]+pe[i]/2)

#plot the energy plot
plt.figure()
plt.plot(T,energy)
plt.title("Energy of N=16 system",fontsize= "13")
plt.xlabel("Time (s)",fontsize= "13")
plt.ylabel("Energy (arbitrary unit)",fontsize= "13")


