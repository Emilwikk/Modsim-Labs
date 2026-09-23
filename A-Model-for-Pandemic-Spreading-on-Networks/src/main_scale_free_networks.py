import config
import utils

import numpy as np

A = utils.powerlaw_adjacency_matrix(config.V,config.pl_gamma)
print(A)
