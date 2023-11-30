"""
Starter code for protein folding
Author: Nicolas Grisuard, based on a script by Paul Kushner


Lab11_Q1
Protein folding
This code is designed to find the energy and the structure of the protein.
Modified by Seunghyun Park
"""
from random import random, randrange, seed
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def calc_energy(monomer_coords, monomer_array):
    """ Compute energy of tertiary structure of protein """
    energy = 0.0

    # compute energy due to all adjacencies (incl. directly bonded monomers)
    for i in range(N):
        for nghbr in [[-1, 0], [1, 0], [0, -1], [0, 1]]:  # 4 neighbours
            nghbr_monomer = monomer_array[monomer_coords[i, 0] + nghbr[0],
                                          monomer_coords[i, 1]+nghbr[1]]

            if nghbr_monomer == 1:  # check neighbour is not empty
                energy += eps

    # divide by 2 to correct for double-counting
    energy = .5*energy

    # correct energy to not count directly bonded monomer neighbours
    energy -= (N-1)*eps

    return energy


def dist(position1, position2):
    """ Compute distance """
    return ((position1[0]-position2[0])**2+(position1[1]-position2[1])**2)**.5


font = {'family': 'DejaVu Sans', 'size': 14}  # adjust fonts
rc('font', **font)
dpi = 150

#Question 2 (d)

seed(8573145) # set seed
eps = -5.0  # interaction energy
N = 30  # length of protein
T_f = 0.5  # temperature for Monte Carlo
n = int(2e6)  # number of Monte Carlo steps
T_steps = 4 
T_i = T_f + T_steps - 1 # Initial temperature
T_array = np.zeros(n) # initialize a array for temperature 
for step in range(T_steps): # temperature steps
    T_array[step*n//T_steps:(step+1)*n//T_steps] = \
        (T_i-T_f)*(1-step/(T_steps-1)) + T_f


energy_array = np.zeros(n)  # initialize array to hold energy

# initialize arrays to store protein information
# 1st column is x coordinates, 2nd column is y coordinates, of all N monomers
monomer_coords = np.zeros((N, 2), dtype='int')

# initialize position of polymer as horizontal line in middle of domain
monomer_coords[:, 0] = range(N//2, 3*N//2)
monomer_coords[:, 1] = N

# 2D array representing lattice,
# equal to 0 when a lattice point is empty,
# and equal to 1 when there is a monomer at the lattice point
monomer_array = np.zeros((2*N+1, 2*N+1), dtype='int')

# fill lattice array
for i in range(N):
    monomer_array[monomer_coords[i, 0], monomer_coords[i, 1]] = 1

# calculate energy of initial protein structure
energy = calc_energy(monomer_coords, monomer_array)

# do Monte Carlo procedure to find optimal protein structure
for j in range(n):
    energy_array[j] = energy

    # move protein back to centre of array
    shift_x = int(np.mean(monomer_coords[:, 0])-N)
    shift_y = int(np.mean(monomer_coords[:, 1])-N)
    monomer_coords[:, 0] -= shift_x
    monomer_coords[:, 1] -= shift_y
    monomer_array = np.roll(monomer_array, -shift_x, axis=0)
    monomer_array = np.roll(monomer_array, -shift_y, axis=1)

    # pick random monomer
    i = randrange(N)
    cur_monomer_pos = monomer_coords[i, :]

    # pick random diagonal neighbour for monomer
    direction = randrange(4)

    if direction == 0:
        neighbour = np.array([-1, -1])  # left/down
    elif direction == 1:
        neighbour = np.array([-1, 1])  # left/up
    elif direction == 2:
        neighbour = np.array([1, 1])  # right/up
    elif direction == 3:
        neighbour = np.array([1, -1])  # right/down

    new_monomer_pos = cur_monomer_pos + neighbour

    # check if neighbour lattice point is empty
    if monomer_array[new_monomer_pos[0], new_monomer_pos[1]] == 0:
        # check if it is possible to move monomer to new position without
        # stretching chain
        distance_okay = False
        if i == 0:
            if dist(new_monomer_pos, monomer_coords[i+1, :]) < 1.1:
                distance_okay = True
        elif i == N-1:
            if dist(new_monomer_pos, monomer_coords[i-1, :]) < 1.1:
                distance_okay = True
        else:
            if dist(new_monomer_pos, monomer_coords[i-1, :]) < 1.1 \
                and dist(new_monomer_pos, monomer_coords[i+1, :]) < 1.1:
                distance_okay = True

        if distance_okay:
            # calculate new energy
            new_monomer_coords = np.copy(monomer_coords)
            new_monomer_coords[i, :] = new_monomer_pos

            new_monomer_array = np.copy(monomer_array)
            new_monomer_array[cur_monomer_pos[0], cur_monomer_pos[1]] = 0
            new_monomer_array[new_monomer_pos[0], new_monomer_pos[1]] = 1

            new_energy = calc_energy(new_monomer_coords, new_monomer_array)

            if random() < np.exp(-(new_energy-energy)/T_array[j]):
                # make switch
                energy = new_energy
                monomer_coords = np.copy(new_monomer_coords)
                monomer_array = np.copy(new_monomer_array)


print()
print('Energy averaged over last quarter of simulations is: {0:.2f}'
      .format(np.mean(energy_array[3*n//4:])))

# %% Question 2 (e) -----------------------------------------------------------|


seed(8573145)
eps = -5.0  # interaction energy
N = 30  # length of protein
T_f = 0.5  # temperature for Monte Carlo
n = 500000  # number of Monte Carlo steps
T_steps = 20
T_i = T_f + T_steps - 10.5 # initial temperature T = 10
T_array = np.zeros(n) # initialize an array for temperature

for step in range(T_steps): # temperature steps
    T_array[step*n//T_steps:(step+1)*n//T_steps] = \
        (T_i-T_f)*(1-step/(T_steps-1)) + T_f

energy_array = np.zeros(n)  # initialize array to hold energy

# initialize arrays to store protein information
# 1st column is x coordinates, 2nd column is y coordinates, of all N monomers
monomer_coords = np.zeros((N, 2), dtype='int')

# initialize position of polymer as horizontal line in middle of domain
monomer_coords[:, 0] = range(N//2, 3*N//2)
monomer_coords[:, 1] = N

# 2D array representing lattice,
# equal to 0 when a lattice point is empty,
# and equal to 1 when there is a monomer at the lattice point
monomer_array = np.zeros((2*N+1, 2*N+1), dtype='int')

# fill lattice array
for i in range(N):
    monomer_array[monomer_coords[i, 0], monomer_coords[i, 1]] = 1

# calculate energy of initial protein structure
energy = calc_energy(monomer_coords, monomer_array)

# do Monte Carlo procedure to find optimal protein structure
for j in range(n):
    energy_array[j] = energy

    # move protein back to centre of array
    shift_x = int(np.mean(monomer_coords[:, 0])-N)
    shift_y = int(np.mean(monomer_coords[:, 1])-N)
    monomer_coords[:, 0] -= shift_x
    monomer_coords[:, 1] -= shift_y
    monomer_array = np.roll(monomer_array, -shift_x, axis=0)
    monomer_array = np.roll(monomer_array, -shift_y, axis=1)

    # pick random monomer
    i = randrange(N)
    cur_monomer_pos = monomer_coords[i, :]

    # pick random diagonal neighbour for monomer
    direction = randrange(4)

    if direction == 0:
        neighbour = np.array([-1, -1])  # left/down
    elif direction == 1:
        neighbour = np.array([-1, 1])  # left/up
    elif direction == 2:
        neighbour = np.array([1, 1])  # right/up
    elif direction == 3:
        neighbour = np.array([1, -1])  # right/down

    new_monomer_pos = cur_monomer_pos + neighbour

    # check if neighbour lattice point is empty
    if monomer_array[new_monomer_pos[0], new_monomer_pos[1]] == 0:
        # check if it is possible to move monomer to new position without
        # stretching chain
        distance_okay = False
        if i == 0:
            if dist(new_monomer_pos, monomer_coords[i+1, :]) < 1.1:
                distance_okay = True
        elif i == N-1:
            if dist(new_monomer_pos, monomer_coords[i-1, :]) < 1.1:
                distance_okay = True
        else:
            if dist(new_monomer_pos, monomer_coords[i-1, :]) < 1.1 \
                and dist(new_monomer_pos, monomer_coords[i+1, :]) < 1.1:
                distance_okay = True

        if distance_okay:
            # calculate new energy
            new_monomer_coords = np.copy(monomer_coords)
            new_monomer_coords[i, :] = new_monomer_pos

            new_monomer_array = np.copy(monomer_array)
            new_monomer_array[cur_monomer_pos[0], cur_monomer_pos[1]] = 0
            new_monomer_array[new_monomer_pos[0], new_monomer_pos[1]] = 1

            new_energy = calc_energy(new_monomer_coords, new_monomer_array)

            if random() < np.exp(-(new_energy-energy)/T_array[j]):
                # make switch
                energy = new_energy
                monomer_coords = np.copy(new_monomer_coords)
                monomer_array = np.copy(new_monomer_array)




mean_energy = [] # create empty list for mean energy
std_energy = [] # create empty list for standard deviation of the energy
T_mean = [] # create empty list for temperature

# calculate the mean energy, standard deviation, and temperature 
for i in range(0,20):
    start = 25000*i
    end = start+24999
    mean_energy.append(np.mean(energy_array[start:end]))
    std_energy.append(np.std(energy_array[start:end]))
    T_mean.append(T_array[start])

# plot mean energy vs temperature with errorbar
plt.figure()
plt.title('Mean energy vs Temperature with erorrbar')
plt.plot(T_mean,mean_energy,label="Mean Energy")
plt.errorbar(T_mean,mean_energy,yerr = std_energy, fmt= 'None', ecolor='red', label="standard deviation")
plt.ylabel('Energy')
plt.xlabel('Temperature')
plt.legend()
plt.grid()
plt.tight_layout()
# plt.savefig('2e.png',dpi=dpi)









