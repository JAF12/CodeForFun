"""
    Created:      Th/09/28/17 (Class)
    Last update:  Th/09/28/17 (Class)
    Author:       Jose Flores
    
    """

#######################################################################

#...Import key external routines:

import numpy as np

#######################################################################

"""This code multipies a 5x3 matrix with a 3x5 matrix. Simple manipulation of the code can create this function to multiply any size matrix multiplications. All that must be changed is the dummy arrays used as well as the size of the matrices. Ex: AB=C"""

#######################################################################

#---- Dummy Arrays used as rows in the Matrices
array1 = np.array([1,2,3])
array2 = np.array([2,3,4])
array3 = np.array([5,3,1])
array4 = np.array([2,4,3])
array5 = np.array([9,8,7])
array6 = np.array([1,2,3,4,5])
array7 = np.array([2,3,4,7,1])
array8 = np.array([5,3,1,4,8])

#---- Number of Rows in matrix A
Rows_A    = 5
#---- Number of Columns in matrix A
Columns_A = 3

#---- Number of Rows in matrix B
Rows_B    = 3

#---- Number of Columns in matrix B
Columns_B = 5

#---- Create Matrix A
A = np.zeros((Rows_A, Columns_A))

#---- Create Matrix B
B = np.zeros((Rows_B,Columns_B))

#---- Updating rows of Matrix A and B with the dummy arrays above

#---- A Matrix
A[0] = array1
A[1] = array2
A[2] = array3
A[3] = array4
A[4] = array5

#---- B Matrix
B[0] = array6
B[1] = array7
B[2] = array8

#---- Resultant Matrix in Matrix Multiplication (Rows_A x Columns_B)
C = np.zeros((A.shape[0],B.shape[1]))

#---- Matrix Multiplication Calculation

#---- Loops through the rows in matrix A
for rows_A in np.arange(A.shape[0]):                #Iteration from: 0-4
    
    #---- Loops through the columns of B
    for columns_B in np.arange(B.shape[1]):         #Iteration from: 0-4
        
        #---- Loops through the rows in matrix B
        for rows_B in np.arange(B.shape[0]):        #Iteration from: 0-3
            
            #---- Stores matrix multiplication results to the resultant matrix constructed above.
            C[rows_A][columns_B] += A[rows_A][rows_B]*B[rows_B][columns_B]

#---- Print Result
print "A:", A
print "B:", B
print "AB=C:",C
