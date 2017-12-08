"""
    Created:      Sat/09/30/17 (Home)
    Last update:  Sat/09/30/17 (Home)
    Author:       Jose Flores
    
    """

#######################################################################

#...Import key external routines:

import numpy as np
import matplotlib.pyplot as plt

#######################################################################

"""This code calculates chaotic behavior of a system and plots logistic maps for the given parameters."""

#######################################################################

#---- This function calculates values for logistic maps.

def calculate(N, mu):

    #---- Construct Nx1 array
    x = np.zeros((N,1))
    
    #---- Insert first element of Nx1 array and create an n array
    x[0] = 0.01
    nth = [0]

    #---- Calculation of logistic map.
    for n in np.arange(N-1):
        x[n+1] = mu * x[n] * (1-x[n])
        nth.append(n+1)
    
    return x, nth, mu

#######################################################################

#---- This function plots a single logistic map curve. Try mu = 3.2, 3.5, 3.6, 3.7, 3.8

def single_plot():
    
    #---- Number of elements in array
    N = 100
    
    #---- Parameter related to "food" that is available according to the book.
    mu = input('Insert mu: ')
    
    #---- Calculates chatic behavior
    x, nth, mu = calculate(N, mu)

    #---- Plots logistic map
    plt.plot(nth,x, 'r-', label=r'$\mu$ $=$' + str(mu) )
    
    plt.legend(loc='lower right', frameon=False, prop={'size': 13})
    plt.title(r'$x_{n+1}$ $=$ $\mu$' + '$(x_{n})$' + '$(1-x_{n})$', fontsize=20)
    plt.ylabel(r'$x_{n}$', fontsize=25)
    plt.xlabel(r'$n$', fontsize=25)
    plt.ylim(0,1)
    plt.show()

#######################################################################

#---- This function plots logistic map curves. mu = 0.5,1.0,1.5,2.0,2.5,3.0

def multiple_mu():
    
    #---- Number of elements in array
    N=100
    
    #---- Parameters related to "food" that is available according to the book.
    total_mu = np.array([0.5,1.0,1.5,2.0,2.5,3.0])
    
    #---- Calculates chaotic behavior + Plots logistic map
    for mu in total_mu:

        if mu == 0.5:
            x, nth, mu = calculate(N, mu)
            plt.plot(nth, x, 'k-', linewidth=1.3,  label=r'$\mu$ $=$' + str(mu))

        if mu == 1.0:
            x, nth, mu = calculate(N,mu)
            plt.plot(nth, x, 'b-', linewidth=1.3,  label=r'$\mu$ $=$' + str(mu))

        if mu == 1.5:
            x, nth, mu = calculate(N, mu)
            plt.plot(nth, x, 'g-', linewidth=1.3,  label=r'$\mu$ $=$' + str(mu))

        if mu == 2.0:
            x, nth, mu = calculate(N,mu)
            plt.plot(nth, x, 'r-', linewidth=1.3,  label=r'$\mu$ $=$' + str(mu))

        if mu == 2.5:
            x, nth, mu = calculate(N, mu)
            plt.plot(nth, x, 'c-', linewidth=1.3,  label=r'$\mu$ $=$' + str(mu))

        if mu == 3.0:
            x, nth, mu = calculate(N,mu)
            plt.plot(nth, x, 'y-', linewidth=1.3,  label=r'$\mu$ $=$' + str(mu))

    #---- Plotting Parameters
    plt.legend(loc='lower right', frameon=False, prop={'size': 13})
    plt.title(r'$x_{n+1}$ $=$ $\mu$' + '$(x_{n})$' + '$(1-x_{n})$', fontsize=20)
    plt.ylabel(r'$x_{n}$', fontsize=25)
    plt.xlabel(r'$n$', fontsize=25)
    plt.show()


#######################################################################

#---- This function plots the logistic maps when mu = 2 and x[0] = 0-1 in intervals of 0.2

def plot2():
    
    #---- Number of elements in array
    N = 100

    #---- Construct Nx1 array
    x = np.zeros((N,1))

    #---- List of the different x[0] I want to use
    x_0 = np.arange(0,1.2,0.2)
    
    #---- Array with the given n's I use
    nth = np.arange(0,100,1)

    #---- Parameter related to "food" that is available according to the book.
    mu = 2

    #---- Loop that runs for each x[0]
    for i in np.arange(len(x_0)):

        x[0] = x_0[i]

        #---- Calculation of logistic map.
        for n in np.arange(N-1):
            x[n+1] = mu * x[n] * (1-x[n])

        #---- Plots logistic maps
        plt.plot(nth,x, '-', label=r'$x_{0}$ $=$' + str(x_0[i]) )
    
    #---- Plotting parameters
    plt.ylim(-1.0,1.5)
    plt.legend(loc='lower right', frameon=False, prop={'size': 13})
    plt.title(r'$x_{n+1}$ $=$ $\mu$' + '$(x_{n})$' + '$(1-x_{n})$', fontsize=20)
    plt.ylabel(r'$x_{n}$', fontsize=25)
    plt.xlabel(r'$n$', fontsize=25)
    plt.ylim(0,1)
    plt.show()

#######################################################################

#---- Calling functions above.
multiple_mu()
single_plot()
plot2()

#######################################################################



