from scipy import sparse, linalg, stats
from scipy.sparse.linalg import svds, aslinearoperator, LinearOperator
import numpy as np

A = np.random.rand(5, 5)
U, S, Vt = np.linalg.svd(A)
print(U)
print(S)
print(Vt)

U, S, Vt = svds(A, 4)
print(U)
print(S)
print(Vt)
np.linalg.qr()