"""
    Created:      Thu/11/02/17 (Class)
    Last update:  Thu/11/02/17 (Class)
    Author:       Jose Flores
    
"""

#######################################################################

"""This program simulates a random walks of tmax steps"""

#######################################################################

import numpy as np
import matplotlib.pyplot as plt

#######################################################################


def random_walk(a=2):


    #---- Number of walks
    numwalks = 10**4

    #---- Number of steps in walk
    tmax = 200

    x2   = np.zeros(tmax)
    y2   = np.zeros(tmax)
    xave = np.zeros(tmax)
    yave = np.zeros(tmax)

    #---- Seed random generator
    np.random.seed(1)

    #---- Random Walk Code
    for j in range(numwalks):
    
        #---- Stores most recent walk
        x = np.zeros([tmax])
        y = np.zeros([tmax])
    
        #---- At time 0
        t=0
        for t in range(tmax):
        
            #---- Random angle from 0 to 2pi and x and y displacements
            p  = np.random.rand() * 2*np.pi
            dx = (a)*np.cos(p)
            dy = (a)*np.sin(p)

            #---- Update Position
            x[t] = x[t-1]+dx
            y[t] = y[t-1]+dy
        
            #---- Update running tally of x^2(t)
            x2[t]   = x2[t]+(x[t])**2
            y2[t]   = y2[t]+(y[t])**2
            xave[t] = xave[t]+x[t]
            yave[t] = yave[t]+y[t]


    x2=x2/numwalks
    y2=y2/numwalks
    xave=xave/numwalks
    yave=yave/numwalks
    r2 = x2+y2

    return x2,y2,xave,yave,r2,x,y, tmax

#######################################################################

#---- Different a parameters
a_s = np.array([2,4])

#---- Empty data container for slopes calculated below
slopes = []

for i in np.arange(len(a_s)):
    print i
    #---- Perform random walk
    x2,y2,xave,yave,r2,x,y, tmax = random_walk(a=a_s[i])

    #---- Step Array used for slope calculation
    step_array = np.linspace(0,tmax,len(r2))

    #---- Slope and Intercept of MSD
    slope, intercept = np.polyfit(step_array,r2,1)
    slopes.append(slope)

    """If you want to see trajectories, simply run the code. If you want to see the MSD plots, comment lines 107-112 and uncomment lines 96-104"""
    
#---- Plotting Linear Plots

#    plt.plot(r2,label=r'$\langle r^{2}_{a= 10}\rangle$')
#    plt.plot(yave,label=r'$\langle y_{a=10}\rangle$')
#    plt.xlabel(r'$Steps$', fontsize = 20)
#    plt.ylabel(r'$\langle y\rangle$ $,$ $\langle r^2\rangle$', fontsize=20)
#    plt.legend(loc='upper right')
#plt.show()

#print "a:", a_s
#print "slope:", slopes

#---- Plotting Trajectory
    plt.plot(x,y,label='a=' + str(a_s[i]))
    plt.xlabel(r'$X-Pos$', fontsize = 20)
    plt.ylabel(r'$Y-Pos$', fontsize=20)
    plt.legend(loc='upper right')
    plt.axis('equal')
plt.show()
