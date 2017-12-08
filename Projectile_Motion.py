"""
    Created:      Tue/09/26/17 (Class)
    Last update:  Tue/09/26/17 (Class)
    Author:       Jose Flores

"""

#######################################################################

#...Import key external routines:

import numpy as np
import matplotlib.pyplot as plt

#######################################################################


"""This code plots the trajectory an object that is lauched at ground level and lands at ground level experiences as time passes. Multiple functions do different things, hope the comments above them explain well."""

#######################################################################

#---- Parameters

#---- Gravitational Acceleration in m/s^2
g = 9.8
    
#---- Angle + velocity in which the projectile is launched
theta_0 = 45.0  # in degrees
v_0     = 20.0  # in m/s
    
#---- Time intervals
dt   = 0.005

#---- Initial Position
x_0 = 0.0
y_0 = 0.0

#######################################################################

#---- Functions: Get either total displacement of the projectile (Range), path (trajectory) or x and y displacements in time.

def Range(v_0, theta_0, g=9.80):
    return (v_0**2 * np.sin(2 * theta_0 * np.pi/180.0)) / g


def x_dist(x_0, v_0, theta_0, time):
    return x_0 + (v_0 * np.cos(theta_0 * np.pi/180))*time


def y_dist(y_0, v_0, theta_0, g, time):
    return y_0 + (v_0 * np.sin(theta_0 * np.pi/180.0))*time - 0.5*g*(time**2)


def Path(theta_0, x, g, v_0):
    return ((np.tan(theta_0 * np.pi/180.0)) * x) - ((g*(x**2))/2*(v_0*np.cos(theta_0 * np.pi/180.0)))


#######################################################################

#---- Function that calculates the path a projectile takes as it's launced from ground level given the parameters above.

def trajectory(theta_0, v_0, dt):

    #---- Initial Time
    time = 0.0

    #---- Empty list for x and y positions
    x          = []
    y          = []
    total_time = []

    #---- Range Equation to calculate the total distance the projectile travels
    distance = Range(v_0, theta_0, g)
    print "Range:", distance

    #---- Looping Equations of motion through multiple time intervals
    while time >=0.0:
    
        #---- Calculating the y position as time increases by dt
        y_distance = y_dist(y_0, v_0, theta_0, g, time)
        
        #---- Storing times into a list
        total_time.append(time)
    
        #---- Appending non-negative values for y postion
        if y_distance >=0.0:
            y.append(y_distance)

        #---- Break loop at negative y positions
        if y_distance < 0.0:
            break

        x_distance = x_dist(x_0, v_0, theta_0, time)
        x. append(x_distance)

        #---- Break loop once the x position is greater than distances traveled
        if x_distance >= Range(v_0, theta_0, g):
            break

        #---- Increase time every dt intervals
        time+=dt

#    print "x positions:", x
#    print "y positions:", y
#    print "time:", total_time
    return  x,y, total_time

#######################################################################

#---- This function plots a single trajectory of an object being thrown from ground level, landing again on ground level.

def plot_single_trajectory():
    
    #---- Plot
    plt.plot(x,y, label = r'$\theta_{o} =$' + str(theta_0) + '$,v_{o} = $' + str(v_0) + '$[ms^{-1}]$')
    plt.title('$Projectile$ $Motion$', fontsize=20)
    
    #---- Axis
    plt.legend(loc='upper right', frameon = False)
    plt.ylabel(r'$Height$ $[m]$', fontsize=20)
    plt.xlabel(r'$Displacement$ $[m]$', fontsize=20)
    plt.show()


#######################################################################

#---- This function plots the projectile motion by having a fixed v_0 and a varying theta_0.

def varrying_angle():
    
    #---- Different angles
    angles = np.array([20,40,60, 70,80])

    for theta_0 in angles:

        if theta_0 == 20:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'k-', linewidth=3.0,  label=r'$\theta$ $=$' + str(theta_0)+'$^{o}$', )

        if theta_0 == 40:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'b-', linewidth=3.0, label=r'$\theta$ $=$' + str(theta_0)+'$^{o}$')

        if theta_0 == 60:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'g-', linewidth=3.0, label=r'$\theta$ $=$' + str(theta_0)+'$^{o}$')
        
        if theta_0 == 70:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'c-', linewidth=3.0, label=r'$\theta$ $=$' + str(theta_0)+'$^{o}$')
        
        if theta_0 == 80:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'r-', linewidth=3.0, label=r'$\theta$ $=$' + str(theta_0)+'$^{o}$')
    
    plt.legend(loc = 'upper right', frameon=False)
    plt.title('$Projectile$ $Motion,$ $v_{o} = 20[m s^{-1}]$', fontsize=20)
    plt.ylabel(r'$Height$ $[m]$', fontsize=20)
    plt.xlabel(r'$Displacement$ $[m]$', fontsize=20)
    plt.show()

#######################################################################

#---- This function plots the projectile motion by having a fixed theta_0 and a varying v_0.

def varrying_velocity():
    
    #---- Different velocities
    velocities = np.array([20,40,60,70,80])

    for v_0 in velocities:
    
        if v_0 == 20:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'k-', linewidth=2.0,  label=r'$v_{o}$ $=$' + str(v_0)+'$[ms^{-1}]$')
    
        if v_0 == 40:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'b-', linewidth=2.0, label=r'$v_{o}$ $=$' + str(v_0)+'$[ms^{-1}]$')
    
        if v_0 == 60:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'g-', linewidth=2.0, label=r'$v_{o}$ $=$' + str(v_0)+'$[ms^{-1}]$')
        
        if v_0 == 70:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'c-', linewidth=2.0, label=r'$v_{o}$ $=$' + str(v_0)+'$[ms^{-1}]$')
    
        if v_0 == 80:
            x,y, total_time = trajectory(theta_0, v_0, dt)
            plt.plot(x, y, 'r-', linewidth=2.0, label=r'$v_{o}$ $=$' + str(v_0)+'$[ms^{-1}]$')
    

    plt.legend(loc = 'upper right', frameon=False, prop={'size': 13})
    plt.title(r'$Projectile$ $Motion,$ $\theta=$' + '$45^{o}$' , fontsize=20)
    plt.ylabel(r'$Height$ $[m]$', fontsize=20)
    plt.xlabel(r'$Displacement$ $[m]$', fontsize=20)
    plt.show()

#######################################################################

#---- Calling functions

#---- Calculates trajectory
x,y, total_time = trajectory(theta_0, v_0, dt)
print "x-pos:", x
print "y-pos:", y
print "time:",  total_time

#---- Produces a single trajectory with the parameters above
plot_single_trajectory()

#---- Varrying velocity and angles plots
varrying_velocity()
varrying_angle()

#######################################################################

