import numpy as np
import matplotlib.pyplot as plt
import therefore

file = np.load('outputs.npz')
therefore.MovieMaker(file['SI'], file['OCI'], file['L'])

#therefore.MovieMaker(scalar_flux, scalar_flux2, L)
