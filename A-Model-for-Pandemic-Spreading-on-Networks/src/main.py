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

sol = solve_ivp(
    fun=utils.reaction_diffusion_IS,
    t_span=t_span_rdIS,
    y0=y0_rdIS,
    method='Radau', # The goat
    args=(config.mu,config.beta),
    dense_output=True
)

# For plotting
t = np.linspace(t_span_rdIS[0], t_span_rdIS[1], 1000)
y = sol.sol(t)

plt.figure()
plt.title("Disease Spreading Simulation")
plt.plot(t/config.seconds_per_day,y[0]/(y[0]+y[1]),label='# Infected')
plt.plot(t/config.seconds_per_day,y[1]/(y[0]+y[1]), label='# Susceptible')
plt.legend()
plt.grid()
plt.xlabel("Time [Days]")
plt.ylabel("Relative Population")
plt.savefig(fig_dir / "fig.png")