import numpy as np
from matplotlib import pyplot as plt


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


#Constants
tend = 5
time_step = 0.1
time_steps = np.int16((tend+1)/time_step)

time = (np.arange(0,tend,time_step))
#print(time)


N = 1000 #Number of individuals
initial_infected = 5 #Number of staring infected



tol = 1e-2
ss_check = np.int64(0.4/time_step)

mu_list = [5, 1, 10]
beta_list = [0.01, 0.02, 0.07]

rel_inf_list = []
rel_sus_list = []
ss_list = []

for mu, beta in list(zip(mu_list, beta_list)):
        #Initialize the populaiton
    rng = np.random.default_rng(42)

    individuals = np.full(N, State.SUSCEPTIBLE)
    infected_idx = rng.choice(N, initial_infected, replace = False)
    individuals[infected_idx] = State.INFECTED

    nmb_susceptible = np.zeros((time.shape))
    nmb_infected = np.zeros((time.shape))

    steady_state = None
    print(ss_check)
    for i, dt in enumerate(time):
        # STEADY STATE CHECKER
        if i > ss_check*2:
            mean_drift = np.mean(nmb_infected[i-ss_check:i]/N)-np.mean(nmb_infected[i-ss_check*2:i-ss_check]/N)
            if np.abs(mean_drift) < tol and steady_state == None:
                steady_state = i - ss_check
                ss_list.append(steady_state)
                print(f"Steady state was reached at time: {time[steady_state]:0.3}, at index: {steady_state}")

        #STATE CHANGER
        nmb_susceptible[i] = np.sum(individuals == State.SUSCEPTIBLE)
        nmb_infected[i] = np.sum(individuals == State.INFECTED)

        mask_susceptible = individuals == State.SUSCEPTIBLE
        mask_infected = individuals == State.INFECTED

        new_susceptible = infected_to_susceptible(mask_infected, N, mu, time_step, rng)
        new_infected = susceptible_to_infected(mask_susceptible, N, beta, time_step, rng)

        individuals[new_infected] = State.INFECTED
        individuals[new_susceptible] = State.SUSCEPTIBLE

        rel_infected = nmb_infected/N
        rel_susceptible = nmb_susceptible/N

        rel_inf_list.append(rel_infected)
        rel_sus_list.append(rel_susceptible)

    plt.figure()
    plt.plot(time, rel_infected, color = "royalblue", linestyle =(0,(1,1)), label = "infected")
    plt.plot(time, rel_susceptible, color = "orange", label = "susceptible")
    plt.axvline(time[steady_state], color="gray", linestyle="--", label="steady state")
    plt.xlabel("Time")
    plt.ylabel("Fraction of individuals")
    plt.title(f"Stochastic Simulation, mu={mu}, beta={beta}")
    plt.grid()
    plt.legend()
    #plt.show()
    plt.savefig(f"sim2_mu{mu}_beta{beta}.pdf", dpi = 300)