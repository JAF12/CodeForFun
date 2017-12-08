"""
    Created:      Sat/11/04/17 (Home)
    Last update:  Sat/11/04/17 (Home)
    Author:       Jose Flores
    
"""

#######################################################################

"""This program simulates NumWalks random walk on a grid of tmax steps"""

#######################################################################

import numpy as np
import matplotlib.pyplot as plt

#######################################################################

#---- Number of walks
numwalks = 5*10**3

#---- Number of steps in walk
tmax     = 10**2

#---- Empty Data Containers
x2   = np.zeros(tmax)
y2   = np.zeros(tmax)
xave = np.zeros(tmax)
yave = np.zeros(tmax)

#---- Seed random generator
np.random.seed(1)

#---- Random Walk Code
for j in range(numwalks):
    
        #---- Empty Data Containers
        x = np.zeros([tmax])
        y = np.zeros([tmax])
    
    #---- At time 0
	t=0
	for t in range(tmax):
        
        #---- Random Number generator of probablities (0-1)
		p = np.random.rand()
		
        #---- Right Probability
		if p<=0.3:
			dx = 1
			dy = 0
        
        #---- Left Probability
		elif p<=0.5:
			dx = -1
			dy = 0
        
        #---- Up Probability
		elif p<=0.8:
			dx = 0
			dy = 1
        
        #---- Down Probability
		else:
			dx = 0
			dy = -1

        #---- Update Position
		x[t] = x[t-1]+dx
		y[t] = y[t-1]+dy
        
        #---- Update running tally of x^2(t)
		x2[t] = x2[t]+(x[t])**2
		y2[t] = y2[t]+(y[t])**2
		xave[t] = xave[t]+x[t]
		yave[t] = yave[t]+y[t]


x2=x2/numwalks
y2=y2/numwalks
xave=xave/numwalks
yave=yave/numwalks
r2 = x2+y2

#---- Step Array used for plotting and slope calculation
step_array = np.linspace(0,tmax,len(r2))

#---- Used for <x^2> - <x>^2 vs steps
x_ave_sq = xave**2
y_ave_sq = yave**2

x_sub = x2 - x_ave_sq
y_sub = y2 - y_ave_sq

#######################################################################

#---- Plotting

#---- Plot: <x>, <y>, <x^2>, <y^2> vs steps
plt.plot(step_array,xave, color='#4B0082',           label=r'$\langle x\rangle$')
plt.plot(step_array,yave, color='#FFFF00', ls ='--', label=r'$\langle y\rangle$')
plt.plot(step_array,x2,   color='#7CFC00', label=r'$\langle x^2\rangle$')
plt.plot(step_array,y2,   color='c',       label=r'$\langle y^2\rangle$')
plt.plot(step_array, x_sub, '--',          label=r'$\langle x^2\rangle$ $-$ $\langle x\rangle^{2}$')
plt.plot(step_array, y_sub, 'k--',         label=r'$\langle y^2\rangle$ $-$ $\langle y\rangle^{2}$')
plt.legend(loc='upper left')
plt.xlabel(r'$Steps$', fontsize = 20)
plt.ylabel(r'$\langle x\rangle$ $,$ $\langle y\rangle$ $,$ $\langle x^2\rangle$ $,$ $\langle y^2\rangle$ $,$ $\langle x\rangle^{2}$ $,$ $\langle y\rangle^{2}$', fontsize=20)
plt.show()


#---- Random walk on grid
plt.plot(x,y,label='Random Walk')
plt.xlabel(r'$X-Pos$', fontsize = 20)
plt.ylabel(r'$Y-Pos$', fontsize=20)
plt.legend(loc='upper left')
plt.show()

