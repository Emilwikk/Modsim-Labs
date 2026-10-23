import time 
from pathlib import Path

import config
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp 
import utils

def susceptible_to_infected(mask_susceptible, N,beta, dt, rng):

    infected = N - np.sum(mask_susceptible)
    prob = 1 - (1 - beta*dt)**infected
    new_infections = mask_susceptible * prob > rng.random(N)

    return new_infections

def infected_to_susceptible(mask_infected,  N, mu, dt, rng):

    prob = mu*dt
    new_susceptible = mask_infected * prob > rng.random(N)

    return new_susceptible

class State:
    INFECTED = 1
    SUSCEPTIBLE = 0

mainPy_dir = Path(__file__).resolve().parent
fig_dir = mainPy_dir / Path(f"../results/stochastic_simulation/")
fig_dir.mkdir(parents=True, exist_ok=True)

#Constants
tend = 5
dt = 0.0001
time_step = dt


time = (np.arange(0,tend,time_step))
#print(time)


#N = np.array([2*4*1000/32, 500, 750]).astype(int)
N = (config.initial_infected + config.initial_susceptible) * config.pop_mult

#initial_infected = N*0.1 #Number of staring infected
#initial_infected = initial_infected.astype(int)

initial_infected = config.initial_infected * config.pop_mult
initial_susceptible = config.initial_susceptible * config.pop_mult

#initial_susceptible = N - initial_infected

tol = 1e-2
ss_check = np.int64(0.4/time_step)*2


t_span_rdIS = [0,config.days_simulated] # seconds
t = np.linspace(t_span_rdIS[0], t_span_rdIS[1], 1000)


mus_days = np.array(config.mus) * (config.pop_mult)
betas_days = np.array(config.betas) 

mus_seconds = mus_days * config.seconds_per_day
betas_seconds = np.array(config.betas)  * config.seconds_per_day

y = []

# Creates pairs of mus and betas, not all possible combos
for i, (mu, beta) in enumerate(zip(mus_days, betas_days)): 
    y0_rdIS = [initial_infected, initial_susceptible]
    sol = solve_ivp(
        fun=utils.reaction_diffusion_IS,
        t_span=t_span_rdIS,
        y0=y0_rdIS,
        method='Radau', # The goat
        args=(mu,beta),
        dense_output=True
    )

    # For plotting
    y.append(sol.sol(t))


rel_inf_list = []
rel_sus_list = []
ss_list = []

for i, (mu, beta) in enumerate(zip(mus_days, betas_days)):
    #Initialize the populaiton
    rng = np.random.default_rng(42)

    individuals = np.full(N, State.SUSCEPTIBLE)
    infected_idx = rng.choice(N, initial_infected, replace = False)
    individuals[infected_idx] = State.INFECTED  

    nmb_susceptible = np.zeros((time.shape))
    nmb_infected = np.zeros((time.shape))

    steady_state = None
    for j, sloppy in enumerate(time):
        # STEADY STATE CHECKER
        if j > ss_check*2:
            mean_drift = np.mean(nmb_infected[j-ss_check:j]/N)-np.mean(nmb_infected[j-ss_check*2:j-ss_check]/N)
            if np.abs(mean_drift) < tol and steady_state == None:
                steady_state = j - ss_check
                ss_list.append(steady_state)
                print(f"Steady state was reached at time: {time[steady_state]:0.3}, at index: {steady_state}")

        #STATE CHANGER
        nmb_susceptible[j] = np.sum(individuals == State.SUSCEPTIBLE)
        nmb_infected[j] = np.sum(individuals == State.INFECTED)

        mask_susceptible = individuals == State.SUSCEPTIBLE
        mask_infected = individuals == State.INFECTED

        new_susceptible = infected_to_susceptible(mask_infected, N, mu, dt, rng)
        new_infected = susceptible_to_infected(mask_susceptible, N, beta, dt, rng)

        individuals[new_infected] = State.INFECTED
        individuals[new_susceptible] = State.SUSCEPTIBLE

    
    rel_infected = nmb_infected/N
    rel_susceptible = nmb_susceptible/N

    rel_inf_list.append(rel_infected)
    rel_sus_list.append(rel_susceptible)

    f1 = y[i]
    spd=config.seconds_per_day
    plt.figure()
    plt.plot(t,f1[0]/(f1[0]+f1[1]),label='ODE_Infected')
    plt.plot(t,f1[1]/(f1[0]+f1[1]), label='ODE_Susceptible')
    plt.plot(time, rel_infected, color = "royalblue", linestyle =(0,(1,1)), label = "Sto_Infected")
    plt.plot(time, rel_susceptible, color = "orange", linestyle =(0,(1,1)), label = "Sto_Susceptible")
    if steady_state is not None:
        plt.axvline(time[steady_state], color="gray",
                    linestyle="--", label="steady state")
    plt.xlabel("Time")
    plt.ylabel("Fraction of individuals")
    plt.title(f"Stochastic Simulation, N={N} mu={mu}, beta={beta}")
    plt.grid()
    plt.legend()
    if steady_state is not None:
        plt.xlim(0, min(time[steady_state] * 1.3, tend))
    else:
        plt.xlim(0, tend)   
    #plt.show()
    plt.savefig(fig_dir / f"N{N}_mu{mu}_beta{beta}.pdf", dpi = 300)