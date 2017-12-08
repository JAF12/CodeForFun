"""
    Created:      Sun/10/22/17 (Home)
    Last update:  Sun/10/22/17 (Home)
    Author:       Jose Flores
    
"""

#######################################################################

"""This program shows how to use the Velocity Verlet method for an orbit"""

#######################################################################

#...Import Modules

import numpy as np
import matplotlib.pyplot as plt

######################################################################

#---- Function that returns the acceleration of the system.

def acceleration(coordinates,params,t):
    
    #---- Positions
    x  = coordinates[0]
    y  = coordinates[1]
    
    #---- Parameter
    MG = params[0]
    
    #---- Radial Distance
    r  = np.sqrt(x**2+y**2)
    
    #---- Acceleration (Newton's Law)
    a  = MG/(r**2)
    
    #---- Acceleration Components
    a_x = -a*x/r
    a_y = -a*y/r
    
    return np.array([a_x,a_y])

######################################################################

#---- Function that calculates the energy of the system

def energy(r,velocities,params,t):
    
    #---- Velocities
    v_x = velocities[0]
    v_y = velocities[1]
    
    #---- Parameter
    MG = params[0]
    
    #---- Total Velcoties
    v  = np.sqrt(v_x**2 + v_y**2)
    
    #---- Kinetic and Potential Energy
    KE = .5 * v**2
    U  = -1.0 * (MG/r)
    
    return KE,U

######################################################################

#---- Function that calculates the angluar momentum of the system

def angular_momentum(coordinates,velocities):
    
    #---- Positions
    x  = coordinates[0]
    y  = coordinates[1]
    
    #---- Velocities
    v_x = velocities[0]
    v_y = velocities[1]
    
    #---- Angular Momentum (Cross products)
    L = x*v_y - y*v_x
    
    return L

######################################################################

#---- Variables

#...Sun (At origin and stationary)
x_position_sun = 0.0
y_position_sun = 0.0
v_x_sun = 0.0
v_y_sun = 0.0

#...Planet 1
x_position_1 = 0.0
y_position_1 = 0.09
v_x_1 = 0.0
r_0 = np.sqrt(x_position_1**2 + y_position_1**2)

#...Circular Oribit
v_y_1 = 2.09 * np.pi / np.sqrt(r_0)

#...Non Circular Orbit
#v_y_1 = 3.0 * np.pi / np.sqrt(r_0)

#---- Parameters
MG = 4*(np.pi**2)
params = np.array([MG],dtype='float')

#---- Time/Step Information
numtimes = 4*10**3
dt         = 0.0001
half_dt    = 0.5 * dt
half_dt_sq = 0.5 * (dt*dt)
times      = np.linspace(0, (numtimes-1)*dt, numtimes)

#---- Degree of Freedom
dof = 2

#---- Empty of Data Array
coordinates    = np.zeros([numtimes,dof])
velocities     = np.zeros([numtimes,dof])

#---- Initialization of Data Arrays
coordinates[0] = np.array([x_position_1,y_position_1])
velocities[0]  = np.array([v_y_1,v_x_1])

#---- Initial Acceleration
a_old = acceleration(coordinates[0],params,0)

######################################################################

#---- Velocity Verlet Calculations
for t in range(1,numtimes):
    
    coordinates[t]  = coordinates[t-1]+velocities[t-1] * dt + half_dt_sq * a_old
    a_new           = acceleration(coordinates[t],params,times[t-1])
    velocities[t]   = velocities[t-1] + (a_new + a_old) * half_dt
    a_old = a_new

#---- Empty + Initialization of distance and angle list
r_dist = [r_0]
theta  = [0]

#---- Updating radial distance and angle list
for i in range(1,numtimes):
    r_dist.append(np.sqrt(coordinates[i,0]**2 + coordinates[i,1]**2))
    theta.append(np.arctan(coordinates[i,1]/coordinates[i,0]))

#---- Empty Lists
Kinetic_Energy   = []
Potential_Energy = []
Angular_Momentum = []

#---- Updating Empty List Above
for i in range(numtimes):
    
    #---- Calculation of Energies and Anglular Momentum
    KE, U = energy(r_dist[i],velocities[i],params,t)
    L     = angular_momentum(coordinates[i],velocities[i])
    
    #---- Updating List
    Kinetic_Energy.append(KE)
    Potential_Energy.append(U)
    Angular_Momentum.append(L)

#Total Energy
total_energy = np.array(Kinetic_Energy) + np.array(Potential_Energy)
print "KE:", np.mean(Kinetic_Energy)
print "U/2", np.mean(Potential_Energy)/2

######################################################################

#---- Plotting

#...Energy vs Time
plt.plot(times,total_energy, label='$<Total-Energy>:$' + str(np.mean(total_energy)))
plt.xlabel('$Time$ $[yr]$', fontsize=18)
plt.ylabel('$Energy$ $[J]$', fontsize=15)
plt.title('$Energy$ $vs$ $Time$', fontsize=18)
plt.legend(loc='upper right')
plt.show()

#...Angular Momentum vs Time
plt.plot(times,Angular_Momentum, label='$<L>=$' + str(np.mean(Angular_Momentum)))
plt.xlabel('$Time$ $[yr]$', fontsize=18)
plt.ylabel(r'$Angluar$ $Momentum$ $[AU^{2}kgs^{-1}]$', fontsize=15)
plt.title('$L$ $vs$ $Time$', fontsize=18)
plt.legend(loc='upper right')
plt.show()

#...Y-pos vs X-pos
plt.plot(coordinates[:,0],coordinates[:,1])
plt.xlabel('$X-Pos$ $[AU]$', fontsize=18)
plt.ylabel('$Y-Pos$ $[AU]$', fontsize=18)
plt.title('$Y-pos$ $vs$ $X-pos$', fontsize=18)
plt.legend(loc='upper right')
plt.show()

#...KE vs Time
plt.plot(times, Kinetic_Energy, label='$<KE>=$' + str(np.mean(Kinetic_Energy)))
plt.xlabel('$Time$ $[Yr]$', fontsize=18)
plt.ylabel('$Kinetic$ $Energy$ $[J]$', fontsize=18)
plt.title('$KE$ $vs$ $Time$', fontsize=18)
plt.legend(loc='upper right')
plt.show()

#...U vs Time
plt.plot(times, Potential_Energy, label='$<U(r)>=$' + str(np.mean(Potential_Energy)))
plt.xlabel('$Time$ $[Yr]$', fontsize=18)
plt.ylabel('$Potential$ $Energy$ $[J]$', fontsize=18)
plt.title('$U(r)$ $vs$ $Time$', fontsize=18)
plt.legend(loc='upper right')
plt.show()








