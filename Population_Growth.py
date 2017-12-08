"""
    Created:      Thu/10/05/17 (Class)
    Last update:  Fri/10/06/17 (Class)
    Author:       Jose Flores
    
"""

#######################################################################

#...Import key external routines:

import numpy as np
import matplotlib.pyplot as plt

###################################################

#---- This function returns the derivative of how a population changes with time.

def derivatives(variables, a,b,t):

    #---- Unpacking variables + parameters
    N  = variables[0]

#    a  = parameters[0]
#    b  = parameters[1]


    #---- Derivative that describes how the number of individuals in a population N vary with time. dN/dt = aN -bN^2
    N_dot = (a * N) - (b * (N**2))
    
    return np.array([N_dot])

###################################################

def part_a():

    #---- Time steps
    dt = 0.001
    numtimes = 50000

    #---- Variables + Parameters
    variables = np.array([10**(-1*2)])

    #---- Establishing time array
    time = np.linspace(0,(numtimes-1)*dt,numtimes)

    #---- Number of variables
    num_variables = len(variables)

    #---- Numbtimes x Variables NxM Matrix with zeros
    Population = np.zeros((numtimes, num_variables))

    #---- Initializing first row of the data matrix with our initial variables
    Population[0] = variables

    #---- Loop that updates the data matrix using Euler's Method

    for t in range(1,numtimes):

        Population[t] = Population[t-1] + (dt*derivatives(Population[t-1], a, b, dt*t))

    #---- Plots Results
    plt.plot(time,Population, label = r'$N(0) = $' + str(variables[0]))

    #---- Extra Plotting Techniques
    plt.ylim(0,2)
    plt.legend( loc='lower right')
    plt.ylabel('$Population,$ $N$', fontsize = 21)
    plt.xlabel('$Time,$ $t$', fontsize = 21)
    plt.title(r'$ dN/dt = aN - bN^{2}$', fontsize = 21)
    plt.text(38, 4, '$a=1,$ $b=1$', fontsize = 18)
    plt.show()

###################################################

#---- This part of the code executes the answers to part a of the HW. Which is a curve that demostrates the solution to the differential equation is exponential at early times.
#---- I decided to get rid of a parameters array and instead insert them individually. This makes the code much easier to use when plotting as seen in part c.
#parameters = np.array([1,1])
a=1
b=1
part_a()

###################################################

def part_b():
    
    #---- Time steps
    dt = 0.001
    numtimes = 50000
    
    #---- Variables + Parameters
    N1 = np.array([300])
    N2 = np.array([600])
    N3 = np.array([900])

    parameters = np.array([1,1])
    
    #---- Establishing time array
    time = np.linspace(0,(numtimes-1)*dt,numtimes)
    
    #---- Number of variables
    num_variables = len(N1)
    
    #---- Numbtimes x Variables NxM Matrix with zeros
    Population1 = np.zeros((numtimes, num_variables))
    Population2 = np.zeros((numtimes, num_variables))
    Population3 = np.zeros((numtimes, num_variables))
    
    #---- Initializing first row of the data matrix with our initial variables
    Population1[0] = N1
    Population2[0] = N2
    Population3[0] = N3

    #---- Loop that updates the data matrix using Euler's Method
    
    for t in range(1,numtimes):
    
    
        Population1[t] = Population1[t-1] + (dt*derivatives(Population1[t-1], a,b, dt*t))
        Population2[t] = Population2[t-1] + (dt*derivatives(Population2[t-1], a,b, dt*t))
        Population3[t] = Population3[t-1] + (dt*derivatives(Population3[t-1], a,b, dt*t))

    #---- Plots Results
    plt.semilogy(time,Population1, label = r'$N=$' + str(N1))
    plt.semilogy(time,Population2, label = r'$N=$' + str(N2))
    plt.semilogy(time,Population3, label = r'$N=$' + str(N3))

    #---- Extra Plotting Techniques
    #plt.ylim(0,10)
    #plt.xlim(0,10)
    plt.legend(loc='lower right')
    plt.ylabel('$Population,$ $N$', fontsize = 21)
    plt.xlabel('$Time,$ $t$', fontsize = 21)
    plt.title(r'$ dN/dt = aN - bN^{2}$', fontsize = 21)
    plt.text(38, 10, '$a=1,$ $b=1$', fontsize = 18)
    plt.show()

###################################################

#---- This part of the code executes the answers to part b of the HW. Creates the Fig 4 in my report. Can plot indivually, just comment out some lines in the function.

#parameters = np.array([1,1])
a=1
b=1
part_b()

###################################################

def part_c(a,b):
    
    #---- Time steps
    dt = 0.001
    numtimes = 50000
    
    #---- Variables + Parameters
    variables = np.array([0.5])
    
    #---- Establishing time array
    time = np.linspace(0,(numtimes-1)*dt,numtimes)
    
    #---- Number of variables
    num_variables = len(variables)
    
    #---- Numbtimes x Variables NxM Matrix with zeros
    Population = np.zeros((numtimes, num_variables))
    
    #---- Initializing first row of the data matrix with our initial variables
    Population[0] = variables
    
    #---- Loop that updates the data matrix using Euler's Method
    
    for t in range(1,numtimes):
        
        Population[t] = Population[t-1] + (dt*derivatives(Population[t-1], a, b, dt*t))
    
        #---- Plots Results
#        plt.semilogy(time,Population, label = r'$b = $' + str(b[1]))
    return time, Population

###################################################

#---- This part of the code executes the answers to part c of the HW. Creates  Fig 5 in my report.

#---- Parameters
a = 1.0
b = np.arange(0,1.1,0.1)

#---- Array with colors used for plotting
colors = np.array(['b-', 'g-', 'r-', 'c-', 'm-', 'y-', 'k-', 'k--', 'r--', 'g--', 'm--'])

#---- Calculationi + plotting
for i in range(len(b)):
    time, Population = part_c(a,b[i])
    plt.semilogy(time,Population, colors[i], label = r'$N(0)=0.5,$ $b = $' + str(b[i]))

#---- Plotting Parameters
plt.ylim(0,20)
plt.xlim(0,30)
plt.legend(loc='upper right')
plt.ylabel('$Population,$ $N$', fontsize = 21)
plt.xlabel('$Time,$ $t$', fontsize = 21)
plt.title(r'$ dN/dt = aN - bN^{2}$', fontsize = 21)
plt.show()

###################################################


