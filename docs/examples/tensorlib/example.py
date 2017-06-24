
import tensorlib
import numpy as np

A = np.arange(4*4*4*4).reshape(4,4,4,4)
B = np.arange(4*4    ).reshape(4,4    )

C = tensorlib.ddot42(A,B)

# check result
Cnp = np.zeros((4,4))
for i in range(4):
  for j in range(4):
    for k in range(4):
      for l in range(4):
        Cnp[i,j] += A[i,j,k,l]*B[l,k]

print('Difference between Python and C++ result: ',np.linalg.norm(C-Cnp))
