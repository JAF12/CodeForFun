"""
    Created:      Thu/10/26/17 (Class)
    Last update:  Thu/10/26/17 (Class)
    Author:       Jose Flores
    
"""

##################################################################

"""This code uses Successive Over-Relaxation to solve Laplace's equation on 2D square grids with arbitrary boundaries.
"""

##################################################################

import numpy as np
import matplotlib.pyplot as plt

##################################################################

#---- Height of Window
H = 100

#---- Width of Window
W = H

#---- Half Length of boundary
R = H/5

#---- Position of circle center
xc = 0.5*W

#---- Position of circle center
yc = 0.5*H

#---- Potential on square
V_square = 1.0

#---- Want to calculate this so that the changes from one iteration to the next are < MaxDeltaV
Max_dV = 0.0001*V_square

#---- Mininum time to run
MinIterations = 0.4*(W+H-R)

#---- How many iterations should go before we print an update?
UpdateFrequency = 40

#---- Over-relaxation parameter
alpha = 2*(1-1.00*np.pi/np.max([W,H]))
#alpha=1
#---- Identifies places to check for boundary condition
BC = np.zeros([W,H])
V  = np.zeros([W,H])

#---- Let's set up BC to have the potential equal V_square for  points less than R from (xc,yc)
numused = 0

##################################################################

#---- Creates grid with potential of 1 in each point
for x in range(1, W-1):
    for y in range(1, H-1):
      
        V[x,y] = V_square
        numused+=1

#---- Creates a box inside the potential. Inside box boundary condition equal 1 and outside 0.
for x in range(30, 71):
    for y in range(30, 71):
        
        BC[x,y] = 1

#---- Now we start iterating
AvgDeltaV = 2*Max_dV
i=0

#---- Iterations
while (AvgDeltaV>Max_dV) or (i<=MinIterations):
    
    #---- Before making changes, sum of changes is 0
    Sum_dV=0
    
    #---- Looping from the origin of box out to its edge
    for x in range(50,W-1):
        for y in range(50,x+1):

               #---- Everywhere outside the box
               if BC[x,y]==0:

               #---- Diagonal Points
                   if x == y:
                      dV=alpha*(0.5*(V[x+1,y]+V[x,y-1])-V[x,y])
                      V[x,y]=V[x,y]+dV
                      Sum_dV+=abs(dV)

               #---- Horizontal Points
                   elif y==0:
                      dV=alpha*(0.25*(2*V[x,y+1] + V[x-1,y] + V[x+1,y])-V[x,y])
                      V[x,y]=V[x,y]+dV
                      Sum_dV+=abs(dV)

                #---- Other points
                   else:
                      dV=alpha*(0.25*(V[x+1,y]+V[x-1,y]+V[x,y+1]+V[x,y-1])-V[x,y])
                      V[x,y]=V[x,y]+dV
                      Sum_dV+=abs(dV)

               #---- Using Symmetry on other quadrants
               V[((W)/2) + (y-((H)/2)), x]           = V[x,y]   #Quadrant I
               V[y,((H-1)/2) - (x-(W/2))]            = V[x,y]   #Lower Quadrant II
               V[x,W-1-y]                            = V[x,y]   #Upper Quadrant II
               V[-y+(H-1), ((H-1)/2) - (x-(W/2))]    = V[x,y]   #Lower Quadrant III
               V[-x+H-1,-y+W-1]                      = V[x,y]   #Upper Quadrant III
               V[-x+H-1,y]                           = V[x,y]   #Lower Quadrant IV
               V[(W/2)-(y-(H/2)),x]                  = V[x,y]   #Upper Quadrant IV

    #---- Calculate the average change
    AvgDeltaV = Sum_dV/numused

    # Print out update
    if (1.0*i)/UpdateFrequency == i/UpdateFrequency:
        print('Iteration %d.  Average voltage change= %f' % (i, AvgDeltaV))
    i +=1

##################################################################

#---- Print statements

print('System size: %d' % (H))
print('Total iterations: %d' % (i))
if i>MinIterations:
    print('This exceeded the minimum number of iterations (%d).' % (MinIterations))
else:
    print('This was the minimum number of iterations.')
print('Average voltage change in last iteration: %.3g' % (AvgDeltaV))

#---- Printing potentials on grid
print V

##################################################################

#---- Plotting Contour plot
#...Plot 1
plt.subplot(211)
plt.contour(V,10)
plt.axis('equal')
plt.title('$Equipotential$ $Line$', fontsize=20)
plt.ylabel('$y$', fontsize=20)

#...Plot 2
plt.subplot(212)
plt.imshow(BC)
plt.axis('equal')

#...Extra
plt.xlabel('$x$', fontsize=20)
plt.ylabel('$y$', fontsize=20)
plt.show()
