import numpy as np
from scipy.linalg import lu, qr, eig
from scipy.sparse import spdiags

def rSVDBKIr(A, k, i, Q=None):
    s = 5
    if Q is None:
        m, n = A.shape
        B = np.random.randn(n, k + s)
        H = np.zeros((m, (k + s) * i))
        H[:, :k + s], _ = lu(A @ B)
        for j in range(2, i + 1):
            start = (k + s) * (j - 2)
            end = (k + s) * (j - 1)
            H[:, start:end], _ = lu(A @ (A.T @ H[:, start:end]))
        Q, _ = qr(H, mode='economic')
        kn = i * (k + s)
    else:
        kn = Q.shape[1]

    T = A.T @ Q
    v, d = eig(T.T @ T)
    ss = np.sqrt(np.diag(d))
    S = spdiags(ss, 0, kn, kn)
    u = (np.linalg.solve(S, T @ v).T).T
    V = u
    x = range(kn - k, kn)
    S = ss[x]
    U = Q @ v[:, x]
    V = V[:, x]

    return U, S, V, Q
