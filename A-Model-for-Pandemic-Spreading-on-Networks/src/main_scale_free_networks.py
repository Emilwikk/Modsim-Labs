import config
import utils

import numpy as np
V = config.V
# Empty adjacency matrix
A_zeros = np.zeros((V,V),dtype=int)

degrees = utils.sample_degrees_powerlaw(A_zeros, config.gamma) 
print(degrees)
weird_list = []
for i,deg in enumerate(degrees):
    [weird_list.append(i) for j in range(deg)]
