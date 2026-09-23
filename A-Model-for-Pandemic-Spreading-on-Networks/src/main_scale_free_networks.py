from pathlib import Path

import config
import utils

import numpy as np

mainPy_dir = Path(__file__).resolve().parent
fig_dir = mainPy_dir / "../results/SFN"
fig_dir.mkdir(parents=True, exist_ok=True)

# Empty list for Scale Free Networks
SFNs = []

# Create the SFNs
for pl_gamma in config.pl_gammas:
    SFNs.append(utils.SFN(config.V, pl_gamma))

# Plot degree distributions
utils.plot_SFNs_degrees(SFNs, fig_dir=fig_dir/"degree_distributions.png")

for sfn in SFNs:
    print(sfn.fitted_gamma)
