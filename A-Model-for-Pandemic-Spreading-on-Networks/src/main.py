import time 
from pathlib import Path

import config
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp 
import utils

mainPy_dir = Path(__file__).resolve().parent
fig_dir = mainPy_dir / "../results/figures"
fig_dir.mkdir(parents=True, exist_ok=True)

t_span_rdIS = [0,config.days_simulated*config.seconds_per_day] # seconds
y0_rdIS = [config.initial_infected,config.initial_susceptible]

mus_seconds = np.array(config.mus)/config.seconds_per_day
betas_seconds = np.array(config.betas)/config.seconds_per_day

print(zip(mus_seconds, betas_seconds))

# Creates pairs of mus and betas, not all possible combos
for mu, beta in zip(mus_seconds, betas_seconds): 
    sol = solve_ivp(
        fun=utils.reaction_diffusion_IS,
        t_span=t_span_rdIS,
        y0=y0_rdIS,
        method='Radau', # The goat
        args=(mu,beta),
        dense_output=True
    )

    # For plotting
    t = np.linspace(t_span_rdIS[0], t_span_rdIS[1], 1000)
    y = sol.sol(t)

    utils.plot_sim1(t,y,mu,beta,spd=config.seconds_per_day, fig_dir=fig_dir, show_plots=config.SHOW_PLOTS)

