import numpy as np
import numba as nb

@nb.jit(nopython=True)
def a_one(x):
    return(x+1)

@nb.jit(nopython=True, parallel=True)
def test():
    N = int(1e3)

    a = np.random.random((N,N)).astype(np.float64)
    b = np.ones(N, dtype=np.float64)

    for i in nb.prange(N):
        a[i,i] = b[i]
    
    #assert((a.all == test.all))

if __name__ == '__main__':
    test()
