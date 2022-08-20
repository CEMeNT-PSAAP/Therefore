import numpy as np
import matplotlib.pyplot as plt
import therefore

file = np.load('outputs.npz')

sim_perams = {'L': 10,
              'dt': .1,
              'offset': .1,
              'ratio': .9}

therefore.MovieMaker(file['SI'], file['OCI'], sim_perams)

#therefore.MovieMaker(scalar_flux, scalar_flux2, L)
