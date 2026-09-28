from pathlib import Path

import config
import utils
from matplotlib import pyplot as plt

mainPy_dir = Path(__file__).resolve().parent
fig_dir = mainPy_dir / "../results/diffusion/"
fig_dir.mkdir(parents=True, exist_ok=True)

# Create and Populate Network 
sfn = utils.SFN(config.V_dan,config.gamma_dan)
sfn.populate(S=config.S_dan, I=config.I_dan)

#for _ in range(config.steps):
#    sfn.move_individuals()

steady_state_reached = False
while steady_state_reached==False:
    sfn.move_individuals()
    steady_state_reached, steady_state = sfn.steady_state_reached(wdw_size=10, tol=0.03)

sfn.plot_populations(out_dir=fig_dir)
sfn.steady_state = steady_state

utils.plot_rho_vs_k(sfn)


