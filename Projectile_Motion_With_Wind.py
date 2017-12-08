"""
    Created:      Sat/10/07/17 (Home)
    Last update:  Sat/10/07/17 (Home)
    Author:       Jose Flores
    
"""

#######################################################################
import numpy as np
import matplotlib.pyplot as plt

#######################################################################

"""This code plots the trajectory a projectile experiences when there is a drag force. Code is separated into two sections. First half does the calculation without an external wind on the projectile. The second half of the code calculates the trajectory with an external wind facing the direction of the projectile."""


#######################################################################


def derivatives(variables,params,t): #We are defining a function that returns derivatives

    #---- Unpacking variables + parameters
    x  = variables[0]
    y  = variables[1]
    v_x = variables[2] * np.cos(theta_0)
    v_y = variables[3] * np.sin(theta_0)
    g          = parameters[0]
    drag_coeff = parameters[1]


    #---- Derivative Calculations
    x_dot = v_x
    y_dot = v_y
    
    vx_dot = -1*drag_coeff*v_x
    vy_dot = -1*(drag_coeff*v_y+g)

    
    return np.array([x_dot, y_dot, vx_dot, vy_dot])

#######################################################################

def projectil(v_0, theta_0):

    #---- Time steps
    dt = 0.001
    numtimes = 10000
    time = np.linspace(0,(numtimes-1)*dt,numtimes)

    #---- Variables + Parameters
    variables = np.array([0.0,0.0,v_0,v_0])

    #---- Number of variables
    num_variables = len(variables)

    #---- Numbtimes x Variables (10000x4) Columns: (X-pos, Y-pos, Vx, Vy)
    data = np.zeros((numtimes, num_variables))

    #---- Inserting initial values to our data matrix
    data[0] = variables

    #---- Loop that calculates the new values in the matrix
    for t in range(1,numtimes):

        #---- Using Euler's Method to solve the differential equation for a particle experinencing a drag force
        data[t] = data[t-1] + (dt*derivatives(data[t-1], parameters, dt*t))

    #---- Conainers for x an y positions
    x_pos = []
    y_pos = []

    #---- Appending containers for the trajectory of the particle to have non-negative values
    for i in range(len(data[::,1])):

        if data[::,1][i] > 0.0:
            x_pos.append(data[::,0][i])
            y_pos.append(data[::,1][i])

    return x_pos, y_pos

#######################################################################

"""Calculations without wind"""

#---- For testing mutliple inital speeds
v_0 = np.array([0.05,0.2,0.5, 0.8])
v_0_new = v_0 + v_0*0.01
theta_0 = 75 * (np.pi/180)

#---- Parameters
parameters = np.array([1.0,1.0])

#---- Colors for plotting
colors = np.array(['r', 'k', 'b', 'g'])
colors2 = np.array(['r--', 'k--', 'b--', 'g--'])

#---- Loop that generates the plots
for i in range(len(v_0)):
#for i in range(len(theta_0s)):

    x_pos, y_pos = projectil(v_0[i], theta_0)
    x_pos2, y_pos2 = projectil(v_0_new[i], theta_0)
    plt.plot(x_pos, y_pos, colors[i], label = r'$v_{o} / v_{T} =$' + str(v_0[i]))
    plt.plot(x_pos2, y_pos2, colors2[i], label = r'$v_{o} / v_{T} =$' + '0.01' + r'$v_{o}/v_{T}$')
    print"1%", i, x_pos2
#---- Plotting Extras
plt.ylabel('$y-Position$ $[m]$', fontsize = 17)
plt.xlabel('$x-Position$ $[m]$', fontsize = 17)
plt.title('$Projectile$ $Motion$ $w/$ $Drag$ $@$' + ' ' +  r'$\theta_{0} = 75^{o}$', fontsize=17)
plt.legend(loc='upper right', frameon=False)
plt.show()
#######################################################################



def derivatives_wind(variables,params,t): #We are defining a function that returns derivatives

    #---- Unpacking variables + parameters
    x  = variables[0]
    y  = variables[1]
    
    #---- Second term corresponds to the wind
    v_x = variables[2] * np.cos(theta_0) + 0.01 * variables[2] * np.cos(theta_0)
    v_y = variables[3] * np.sin(theta_0) + 0.01 * variables[2] * np.cos(theta_0)
    
    
    g          = parameters[0]
    drag_coeff = parameters[1]


    #---- Calculations
    x_dot = v_x
    y_dot = v_y
    
    vx_dot = -1*drag_coeff*v_x
    vy_dot = -1*(drag_coeff*v_y+g)

    
    return np.array([x_dot, y_dot, vx_dot, vy_dot])



#######################################################################

def projectile_wind(v_0, theta_0):

    #---- Time steps
    dt = 0.001
    numtimes = 10000
    time = np.linspace(0,(numtimes-1)*dt,numtimes)

    #---- Variables + Parameters
    variables = np.array([0.0,0.0,v_0,v_0])

    #---- Number of variables
    num_variables = len(variables)

    #---- Numbtimes x Variables (10000x4) Columns: (X-pos, Y-pos, Vx, Vy)
    data = np.zeros((numtimes, num_variables))

    #---- Inserting initial values to our data matrix
    data[0] = variables

    #---- Loop that calculates the new values in the matrix
    for t in range(1,numtimes):

        #---- Using Euler's Method to solve the differential equation for a particle experinencing a drag force
        data[t] = data[t-1] + (dt*derivatives_wind(data[t-1], parameters, dt*t))

    #---- Conainers for x an y positions
    x_pos = []
    y_pos = []

    #---- Appending containers for the trajectory of the particle to have non-negative values
    for i in range(len(data[::,1])):

        if data[::,1][i] > 0.0:
            x_pos.append(data[::,0][i])
            y_pos.append(data[::,1][i])

    return x_pos, y_pos

#######################################################################

"""Calculations with wind added"""

#---- For testing mutliple inital speeds
v_0 = np.array([0.05,0.2,0.5, 0.8])
#v_0_new = v_0 + v_0*0.01
theta_0 = 75 * (np.pi/180)

#---- Parameters
parameters = np.array([1.0,1.0])

#---- Colors for plotting
colors = np.array(['r', 'k', 'b', 'g'])
colors2 = np.array(['r--', 'k--', 'b--', 'g--'])

#---- Loop that generates the plots
for i in range(len(v_0)):

    x_pos, y_pos = projectil(v_0[i], theta_0)
    x_pos2, y_pos2 = projectile_wind(v_0[i], theta_0)
    plt.plot(x_pos, y_pos, colors[i], label = r'$v_{o} / v_{T} =$' + str(v_0[i]))
    plt.plot(x_pos2, y_pos2, colors2[i], label = r'$v_{o} / v_{T} = $' +str(v_0[i])+ r'+$v_{wind}$')
#    print"no", i, x_pos
#    print"yes", i, x_pos2

#---- Plotting Extras
plt.ylabel('$y-Position$ $[m]$', fontsize = 17)
plt.xlabel('$x-Position$ $[m]$', fontsize = 17)
plt.title('$Projectile$ $Motion$ $w/$ $Drag$ $+Wind$ $@$' + ' ' +  r'$\theta_{0} = 75^{o}$', fontsize=17)
plt.legend(loc='upper right', frameon=False)
plt.show()

#######################################################################











