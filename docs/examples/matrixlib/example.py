
import matrixlib
import numpy as np

A = np.arange(4*5*6*7).reshape(4,5,6,7)

C = matrixlib.squareRoot(A)

# check result
Cnp = np.sqrt(A)

print('Difference between Python and C++ result: ',np.linalg.norm(C-Cnp))

print('')
print('The result')
print(C)
