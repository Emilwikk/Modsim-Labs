from pathlib import Path

import config
import utils
from matplotlib import pyplot as plt
import numpy as np
import shutil
from pathlib import Path

mu_list =  [0.5]    # [2,     2,      0.5]
gamma_list = [-3.6]  # [-1.5,  -3.5,   -1.5]
beta_list =    ["beta_i"] #[1,     1,      "beta_i"]
timestep_list = [0.97]#[0.22, 0.22, 0.97]

rhos = np.linspace(config.rho_min,config.rho_max,config.n_rhos)

n_jobs = len(config.V_list)*len(rhos)*len(mu_list)
job_n = 1

for gamma, mu, beta, timestep in zip(gamma_list, mu_list, beta_list, timestep_list, strict=True):
    SFNs = []
    mainPy_dir = Path(__file__).resolve().parent
    fig_dir_real = mainPy_dir / f"../results/combining_mu_{mu}_gamma_{gamma}_beta_{beta}/"

    #Remove folder if it exists already
    if fig_dir_real.exists():
        shutil.rmtree(fig_dir_real)

    fig_dir_real.mkdir(parents=True, exist_ok=True)

    if isinstance(beta, str):
        print("Warning: beta is now a string!\n")
    else:
        print(f"Warning: beta is now {beta}\n")
    for V in config.V_list:

        init_Ss = (1-config.init_proportion_infected)*rhos*V
        init_Is = config.init_proportion_infected*rhos*V

        init_conds = zip(init_Ss, init_Is)

        for init_cond in init_conds:
            print(f"Running simulation [{job_n}/{n_jobs}]")
            #fig_dir =fig_dir_real
            #fig_dir = fig_dir / Path(f"V_{V}_S_{init_cond[0]:.1f}_I_{init_cond[1]:.1f}/")
            #fig_dir.mkdir(parents=True, exist_ok=True)
                
            # Create and Populate Network 
            sfn = utils.SFN(V,gamma)
            sfn.populate(S=init_cond[0], I=init_cond[1])

            steady_state_reached = False
            step = 0
            while steady_state_reached==False:
                print(f"---Step: {step}---", end="\r", flush=True)
                if isinstance(beta, str):
                    utils.stochastic_state(sfn, beta=None, mu=mu, timestep = timestep)
                else:
                    utils.stochastic_state(sfn, beta=beta, mu=mu, timestep = timestep)
                sfn.move_individuals()
                steady_state_reached, steady_state = sfn.infected_steady_state_reached(10, 0.01, 0.0, step)
                step = step+1
            print(f"Steady State Reached After {step} Steps")

            #sfn.plot_populations(out_dir=fig_dir)
            #sfn.steady_state = steady_state

            #utils.plot_rho_vs_k(sfn, out_dir=fig_dir)

            SFNs.append(sfn)
            print(f"Simulation [{job_n}/{n_jobs}] complete\n")
            job_n += 1


    fig_dir_real.mkdir(parents=True, exist_ok=True)
    utils.plot_rel_infect_vs_init_rho(sfn_list=SFNs, out_dir=fig_dir_real, gamma=gamma, mu=mu, beta=beta)
