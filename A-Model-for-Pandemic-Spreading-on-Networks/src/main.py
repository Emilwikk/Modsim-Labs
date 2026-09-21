import time 
from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
import utils
from scipy.integrate import solve_ivp 

mainPy_dir = Path(__file__).resolve().parent
fig_dir = mainPy_dir / "../results/figures"
fig_dir.mkdir(parents=True, exist_ok=True)
seconds_per_day = 24*60*60

days = 5
t_span_rdIS = [0,days*seconds_per_day] # seconds
y0_rdIS = [20,20] # I,S

mu = 2e-6   # Recovery rate [s^-1]
beta = 1e-6  # Infection rate [s^-1]

sol = solve_ivp(
    fun=utils.reaction_diffusion_IS,
    t_span=t_span_rdIS,
    y0=y0_rdIS,
    method='Radau', # The goat
    args=(mu,beta),
    dense_output=True
)
t = np.linspace(t_span_rdIS[0], t_span_rdIS[1], 1000)
y = sol.sol(t)

plt.figure()
plt.plot(t/seconds_per_day,y[0],label='# Infected')
plt.plot(t/seconds_per_day,y[1], label='# Susceptible')
plt.xlabel("Time [Days]")
plt.savefig(fig_dir / "fig.png")