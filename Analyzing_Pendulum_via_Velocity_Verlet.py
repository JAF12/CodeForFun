"""
    Created:      Tue/10/10/17 (Class)
    Last update:  Tue/10/10/17 (Class)
    Author:       Jose Flores
    
"""

#######################################################################

import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fourier

#######################################################################

"""This program shows how to use the Velocity Verlet method for a pendulum"""

#######################################################################

def acceleration(coordinates,params,t):
    
    #---- Unpacking Variables
    g               = params[0]
    length_pendulum = params[2]
    theta           = coordinates[0]
    
    #---- d^2/dt^2 (theta) = -g/l * sin(theta)
    vdot = -(g*np.sin(theta))/length_pendulum
    
    return vdot

#######################################################################

def energy(coordinates,velocities,params,t):

    #---- Unpacking Variables
    g                   = params[0]
    mass                = params[1]
    length_pendulum     = params[2]
    Mom_Inertia_pt_part = params[3]
    theta               = coordinates[0]
    omega               = velocities[0]

    #---- Kinetic and Potential Energy
    KE = 0.5 * Mom_Inertia_pt_part * omega**2
    U  = (mass*g*length_pendulum)*(1-np.cos(theta))
    
    return KE+U

#######################################################################


#---- Variables/Initial Conditions

#...Gravitional Acce (The 4*pi**2 is because we work in units where omega = 2pi, so the period of small oscillations is 1)
g = 4*np.pi**2

#....Pendulum Variables
mass = 1.0
length_pendulum = 1.0
Mom_Inertia_pt_part = (mass)*(length_pendulum**2)

#...Angular Velocity
omega_0 = 0.0

#....Intial Angle
theta_0 = 95 * (np.pi/180)

#...Degrees of Freedom
dof = 1

#---- Time Steps
dt = 0.01
numtimes = 4*10**3
half_dt = 0.5 * dt
half_dt_sq = 0.5 * (dt**2)
times = np.linspace(0, (numtimes-1)*dt, numtimes)

#---- Parameters
params = np.array([g, mass, length_pendulum, Mom_Inertia_pt_part])


#---- Establishing Data Arrays
coordinates   = np.zeros([numtimes,dof])
velocities    = np.zeros([numtimes,dof])
accelerations = np.zeros([numtimes,dof])
energies      = np.zeros(numtimes)

#---- Initialize the arrays
i=0
coordinates[i] = theta_0
velocities[i]  = omega_0
accelerations[i] = acceleration(coordinates[i], params, 0)
energies[i]    = energy(coordinates[i],velocities[i], params,0)

#######################################################################

#---- Velocity Verlet Method

for t in range(1,numtimes):
    
    #---- Updating Coordinates
    coordinates[t] = coordinates[t-1] + velocities[t-1]*dt + half_dt_sq * accelerations[t-1]

    #---- Current Acceleration
    acc_new = acceleration(coordinates[t], params, times[t])
    
    #---- Updating Velocities
    velocities[t] = velocities[t-1] + half_dt * (acc_new + accelerations[t-1])

    #---- Calculating Energy
    energies[t] = energy(coordinates[t],velocities[t], params,times[t])

    accelerations[t] = acc_new


#######################################################################

#---- Computing Freq Spectrum

#...Normalizing the spectrum
spectrum = fourier.fft(coordinates[:,0]/numtimes)

#...Shift things so that the center frequency is zero
spectrum = fourier.fftshift(spectrum)

#...Create an array for the frequency
freq = np.linspace(-0.5/dt,0.5/dt,numtimes)

#######################################################################

#---- Plotting Info

plt.figure()
plt.subplot(311)
plt.plot(times,coordinates[:,0],label=r'$\theta$')  #Make a plot
plt.xlabel("Time")
plt.ylabel("Angle")
plt.legend(loc = 'upper right')
plt.subplot(312)
plt.plot(times,energies,label='Energy (Euler-Cromer)')
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend(loc = 'upper right')
plt.subplot(313)
plt.plot(freq, abs(spectrum)**2,label='Amplitude Spectrum')
plt.xlim(-3.1,3.1)
plt.xlabel("Frequency")
plt.ylabel("Power")
plt.legend(loc = 'upper right')
plt.show()

#######################################################################
plt.plot(times, energies)
plt.title('$Velocity-Verlet$', fontsize=20)
plt.ylabel('$Energy$ $[J]$', fontsize=16)
plt.xlabel('$Time$ $[sec]$', fontsize=20)
plt.show()
