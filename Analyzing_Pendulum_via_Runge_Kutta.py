"""
    Created:      Thu/10/12/17 (Class)
    Last update:  Thu/10/12/17 (Class)
    Author:       Jose Flores
    
"""

#######################################################################


import numpy as np
import matplotlib.pyplot as plt

#######################################################################


def derivatives(variables, params, t):
	
    #---- Unpacking Variables
    g               = params[0]
    length_pendulum = params[2]
    theta           = variables[0]
    omega           = variables[1]
    
    #---- d/dt (theta) = w
    theta_dot = omega
    
    #---- d^2/dt^2 (theta) = -g/l * sin(theta)
    omega_dot  = -(g*np.sin(theta))/length_pendulum
    
    return np.array([theta_dot,omega_dot])

#######################################################################


def Energy(variables, params, t):
    
    
    #---- Unpacking Variables
    g                   = params[0]
    mass                = params[1]
    length_pendulum     = params[2]
    Mom_Inertia_pt_part = params[3]
    theta               = variables[0]
    omega               = variables[1]
    
    #---- Kinetic and Potential Energy
    KE = 0.5 * Mom_Inertia_pt_part * omega**2
    U  = (mass*g*length_pendulum)*(1-np.cos(theta))
    
    return KE+U

#######################################################################

#---- Variables/Initial Conditions

#...Gravitional Acce (The 4*pi**2 is because we work in units where omega = 2pi, so the period of small oscillations is 1)
g = 4 * np.pi**2

#....Pendulum Variables
mass = 1.0
length_pendulum = 1.0
Mom_Inertia_pt_part = (mass)*(length_pendulum**2)

#...Angular Velocity
omega_0 = 0.0

#....Intial Angle
theta_0 = 175 * (np.pi/180)

#...Degrees of Freedom
dof = 2

#---- Time Steps
t_0 = 0
dt = 0.01
numtimes = 4*10**3
times = np.linspace(t_0, (numtimes-1)*dt, numtimes)

#---- Parameters
params    = np.array([g, mass, length_pendulum, Mom_Inertia_pt_part])

#---- Establishing Data Arrays
variables = np.zeros([numtimes,dof],dtype=float)
Energies = np.zeros(numtimes,dtype=float)

#---- Initialize the arrays
variables[0] = np.array([theta_0,omega_0])
Energies[0]  = Energy(variables[0],params,times[0])

#######################################################################

#---- Runge-Kutta Method

for step in range(1,numtimes):
    
	#----Calculate derivatives at different times, multiply the derivatives by the time step, to get the changes in the variables.
    k1 = dt*derivatives(variables[step-1],params,times[step-1])
    k2 = dt*derivatives(variables[step-1]+k1/2,params,times[step-1]+dt/2)
    k3 = dt*derivatives(variables[step-1]+k2/2,params,times[step-1]+dt/2)
    k4 = dt*derivatives(variables[step-1]+k3,params,times[step-1]+dt)

    #---- Runge-Kutta estimate of the new values of theta and d(theta)/dt
    variables[step] = variables[step-1] + k1/6 + k2/3 + k3/3 + k4/6
    Energies[step]  = Energy(variables[step],params,times[step])

#######################################################################

#---- First Plot
plt.figure(1)
plt.subplot(211)
plt.plot(times,variables[:,0],label='Theta (RK4)')

#---- Extra Plotting Things
plt.xlabel("Time")
plt.ylabel("Theta")
plt.legend(loc = 'upper right')

#---- Second Plot
plt.subplot(212)
plt.plot(times,Energies,label='Energy (RK4)')

#---- Extra Plotting Things
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend(loc = 'upper right')
plt.show()

#######################################################################

plt.plot(times, Energies)
plt.title('$Runge-Kutta$', fontsize=20)
plt.ylabel('$Energy$ $[J]$', fontsize=16)
plt.xlabel('$Time$ $[sec]$', fontsize=20)
plt.show()


