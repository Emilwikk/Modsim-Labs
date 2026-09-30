from pathlib import Path

import config
import utils
from matplotlib import pyplot as plt
import numpy as np

mainPy_dir = Path(__file__).resolve().parent
fig_dir = mainPy_dir / "../results/combining_it_all/"
fig_dir.mkdir(parents=True, exist_ok=True)

rhos = np.linspace(0.1,2.0,8)

for V in config.V_list:
    init_Ss = (1-config.init_proportion_infected)*rhos*V
    init_Is = config.init_proportion_infected*rhos*V

    init_conds = zip(init_Ss, init_Is)

    for init_cond in init_conds:

        # Create and Populate Network 
        sfn = utils.SFN(V,config.gamma5)
        sfn.populate(S=init_cond[0], I=init_cond[1])

        steady_state_reached = False
        step = 0
        while steady_state_reached==False:
            print(f"---Step: {step}---", end="\r", flush=True)

            utils.stochastic_state(sfn, beta=config.beta5, mu=config.mu5, timestep = 0.25)
            sfn.move_individuals()
            steady_state_reached, steady_state = sfn.steady_state_reached(wdw_size=10, tol=0.05, N = init_cond[0]+init_cond[1], V=V)
            step = step+1
        print(f"Steady State Reached After {step} Steps")

        sfn.plot_populations(out_dir=fig_dir)
        sfn.steady_state = steady_state

        utils.plot_rho_vs_k(sfn, out_dir=fig_dir)


