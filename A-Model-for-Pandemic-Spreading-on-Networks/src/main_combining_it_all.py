from pathlib import Path

import config
import utils
from matplotlib import pyplot as plt
import numpy as np

mainPy_dir = Path(__file__).resolve().parent
fig_dir = mainPy_dir / "../results/combining_it_all/"
fig_dir.mkdir(parents=True, exist_ok=True)

rhos = np.linspace(config.rho_min,config.rho_max,config.n_rhos)

n_jobs = len(config.V_list)*len(rhos)
job_n = 1

SFNs = []

for V in config.V_list:

    init_Ss = (1-config.init_proportion_infected)*rhos*V
    init_Is = config.init_proportion_infected*rhos*V

    init_conds = zip(init_Ss, init_Is)

    for init_cond in init_conds:
        print(f"Running simulation [{job_n}/{n_jobs}]")
        fig_dir = mainPy_dir / "../results/combining_it_all/"
        fig_dir = fig_dir / Path(f"V_{V}_S_{init_cond[0]:.1f}_I_{init_cond[1]:.1f}/")
        fig_dir.mkdir(parents=True, exist_ok=True)
            
        # Create and Populate Network 
        sfn = utils.SFN(V,config.gamma5)
        sfn.populate(S=init_cond[0], I=init_cond[1])

        steady_state_reached = False
        step = 0
        while steady_state_reached==False:
            print(f"---Step: {step}---", end="\r", flush=True)

            utils.stochastic_state(sfn, beta=config.beta5, mu=config.mu5, timestep = 0.25)
            sfn.move_individuals()
            steady_state_reached, steady_state = sfn.steady_state_reached(wdw_size=1000, tol=0.05, N = init_cond[0]+init_cond[1], V=V)
            step = step+1
        print(f"Steady State Reached After {step} Steps")

        sfn.plot_populations(out_dir=fig_dir)
        sfn.steady_state = steady_state

        utils.plot_rho_vs_k(sfn, out_dir=fig_dir)

        SFNs.append(sfn)
        print(f"Simulation [{job_n}/{n_jobs}] complete\n")
        job_n += 1

fig_dir = mainPy_dir / "../results/combining_it_all/"
fig_dir.mkdir(parents=True, exist_ok=True)
utils.plot_rel_infect_vs_init_rho(sfn_list=SFNs, out_dir=fig_dir)
