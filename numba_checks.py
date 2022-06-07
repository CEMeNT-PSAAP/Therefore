import numba as nb
import numpy as np

@nb.jit(nopython=True)
def zeros(size):
    size = int(size)
    x = np.zeros((size,size), dtype=np.float64)
    return(x)
    
if __name__ == '__main__':
    x = zeros(18)
    print(x.shape)
