import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import svds
from numpy.linalg import norm
# from multiprocessing import Pool
from rSVDBKIr import rSVDBKIr


def fastSVT_U(M, tol, ran, i_reuse, q_reuse, delta):
    # Initial setup
    m, n = M.shape
    Omega = M != 0
    Ns = Omega.sum()
    if delta is None:
        delta = 1.2 * m * n / Ns

    xi, yi = Omega.nonzero()
    tau = 5 * n
    l = 5
    i_max = 1000
    PM = M.copy()
    normPM2 = svds(PM, 1, return_singular_vectors=False)[0]
    normPM = norm(PM, 'fro')
    k0 = np.ceil(tau / delta / normPM2)
    Y0 = k0 * delta * PM

    dec = 0
    r = 0
    p = 2
    q = 0
    err_before = 1000

    # Main loop
    for i in range(1, i_max + 1):
        r_before = r
        r += 1
        if i % 50 == 0:
            delta /= 1.1

        if i > i_reuse and q < q_reuse:
            U, S, V, Q = rSVDBKIr(Y0, r, p, U)
            q += 1
        else:
            U, S, V, Q = rSVDBKIr(Y0, r, p)
            q = 0

        while S[0] > tau:
            r += l
            U, S, V, Q = rSVDBKIr(Y0, r, p)

        j = next((j for j, val in enumerate(S) if val > tau), r)
        r_max = r
        r = max(r_max - j + 1, r_before)
        x = range(r_max - r, r_max)
        S[list(x)] -= tau

        # Parallel processing (requires proper translation of parfor loop)
        # This section needs to be adapted for Python's multiprocessing capabilities

        # Create X
        x_now = np.zeros(Ns)
        for j in range(Ns):
            temp = U[xi[j], list(x)] @ V[yi[j], list(x)].T
            x_now[j] = np.clip(temp, ran[0], ran[1])

        X = csr_matrix((x_now, (xi, yi)), shape=(m, n))
        PX = X - PM
        err = norm(PX, 'fro') / normPM

        if err <= tol:
            X = U[:, list(x)] @ V[:, list(x)].T
            X = np.clip(X, ran[0], ran[1])
            k = r
            iters = i
            break

        if err > err_before:
            dec = 0
            p += 1
            q = 10
        else:
            if p > 5:
                dec += 1
                if dec == 10:
                    p -= 1
                    dec = 0

        err_before = err
        print([i, r, err, p])
        Y0 -= delta * PX

    return X, iters, k

# Example usage (assuming M is a sparse matrix):
# X, iters, k = fastSVT_U(M, tol, ran, i_reuse, q_reuse, delta)
