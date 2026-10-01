from pathlib import Path

from collections import Counter
from matplotlib import pyplot as plt
import numpy as np
import random

def reaction_diffusion_IS(t, y, mu, beta):
    I, S = y

    dIdt = -mu*I + beta*I*S
    dSdt = mu*I - beta*I*S
    return [dIdt, dSdt]

def plot_sim1(t, y, mu, beta, N, spd, fig_dir, show_plots):
    plt.figure()
    plt.title(f"Disease Spreading Simulation \nwith mu={mu:.2e} and beta={beta:.2e}")
    plt.plot(t/spd,y[0]/(y[0]+y[1]),label='# Infected', color="r")
    plt.plot(t/spd,y[1]/(y[0]+y[1]), label='# Susceptible', color="b")
    plt.axhline(y=(N-mu/beta)/N, label="Theoretical Steady State #I", color="r", linestyle="dotted")
    plt.axhline(y=(mu/beta)/N, label="Theoretical Steady State #S", color="b", linestyle="dotted")
    plt.legend()
    plt.grid()
    plt.xlabel("Time [Days]")
    plt.ylabel("Relative Population")
    plt.savefig(fig_dir / f"sim1_mu{mu:.1e}_beta{beta:.1e}.pdf")
    if show_plots == 1:
        plt.show()


def powerlaw_adjacency_matrix(V, gamma):

    print("Building adjacency matrix...")

    A = np.zeros((V, V), dtype=int)

    k_min = 2
    k_max = np.sqrt(V)

    target_degrees = []

    for _ in range(V):
        r = np.random.rand()

        ki = np.floor(
            (
                (k_max**(gamma + 1) - k_min**(gamma + 1)) * r
                + k_min**(gamma + 1)
            ) ** (1 / (gamma + 1))
        )

        target_degrees.append(int(ki))

    weird_list = []

    for i, deg in enumerate(target_degrees):
        weird_list.extend([i] * deg)

    failed_attempts = 0
    max_failed_attempts = 1000

    while len(weird_list) > 1:

        idx1, idx2 = random.sample(range(len(weird_list)), 2)

        i = weird_list[idx1]
        j = weird_list[idx2]

        if i != j and A[i, j] == 0:
            A[i, j] = 1
            A[j, i] = 1

            # Remove by index so we don't search through the list
            for idx in sorted((idx1, idx2), reverse=True):
                weird_list.pop(idx)

            failed_attempts = 0

        else:
            failed_attempts += 1
            print(f"failed attempt: {failed_attempts}", end="\r", flush=True) 

        if failed_attempts >= max_failed_attempts:
            print("Broke")
            break
            

    print("Adjacency matrix complete")

    return A, target_degrees

def fit_gamma(degrees):

    degree, count = np.unique(degrees,return_counts=True)
    slope,inter = np.polyfit(np.log(degree),np.log(count),deg=1)
    fitted_gamma = slope
    return fitted_gamma

class Node:
    def __init__(self,S,I):
        self.S = S
        self.I = I
        self.Pop = S+I
        self.degree = None
        self.beta = 1/(self.Pop)

        self.S_history = [S]
        self.I_history = [I]
        self.Pop_history = [S+I]

        self.steady_counter = 0
        self.steady_state_reached = False

class SFN:
    def __init__(self, V, gamma):
        self.V = V
        self.gamma = gamma
        A, target_degrees = powerlaw_adjacency_matrix(self.V,self.gamma)
        self.A = A
        self.target_degrees = target_degrees
        self.actual_degrees = np.sum(self.A, axis=0)
        self.avg_network_degree = np.mean(self.actual_degrees)
        self.fitted_gamma = fit_gamma(self.target_degrees)
        self.I_total_history = []

    def populate(self, S, I):
        self.nodes = []
        probabilities = np.ones(self.V)/self.V
        #probabilities = np.zeros(self.V)
        #probabilities[0] = 1

        S_nodes = np.random.multinomial(S,probabilities)
        I_nodes = np.random.multinomial(I,probabilities)

        for i in range(self.V):
            self.nodes.append(Node(S_nodes[i],I_nodes[i]))
            self.nodes[i].degree = np.sum(self.A[i,:])
            self.nodes[i].beta = 1/(self.nodes[i].I+self.nodes[i].S)
        self.avg_node_population = (S+I)/self.V
        self.I_total_history.append(I)

        

    def move_individuals(self):
        S_move = np.zeros(self.V, dtype=int)
        I_move = np.zeros(self.V, dtype=int)
        for n, node in enumerate(self.nodes):
            if node.degree==0:
                S_move[n] += node.S
                I_move[n] += node.I 
                continue 

            move_prob = self.A[n,:]/node.degree

            # Move susceptible individuals to new node
            S_move = S_move + np.random.multinomial(node.S, move_prob)
            I_move = I_move + np.random.multinomial(node.I, move_prob)
        I_hist = 0
        # Add changes to actual nodes
        for n, node in enumerate(self.nodes):
            node.S = S_move[n]
            node.I = I_move[n]
            node.Pop = node.S+node.I
            node.beta = 1/node.Pop

            # Save to history
            node.S_history.append(node.S)
            node.I_history.append(node.I)
            node.Pop_history.append(node.Pop)
            I_hist = I_hist + node.I
        self.I_total_history.append(I_hist)

    def plot_populations(self, out_dir:Path):
        if len(self.nodes[0].S_history)>1:
            plt.figure()
            plt.title(f"Number of Susceptible Individuals Per Node Over Time\nGamma = {self.gamma}")
            for n, node in enumerate(self.nodes):
                plt.plot(range(len(node.S_history)), node.S_history, label=f"Node {n}")
            #plt.legend()
            #plt.axis([0,config.steps+1,0,400])
            plt.xlabel("Time []")
            plt.ylabel("Quantity")
            plt.grid()
            plt.savefig(out_dir/f"susceptible_over_time_gamma_{self.gamma}.pdf")

            plt.figure()
            plt.title(f"Number of Infected Individuals Per Node Over Time\nGamma = {self.gamma}")
            for n, node in enumerate(self.nodes):
                plt.plot(range(len(node.I_history)), node.I_history, label=f"Node {n}")
            #plt.legend()
            #plt.axis([0,config.steps+1,0,400])
            plt.xlabel("Time []")
            plt.ylabel("Quantity")
            plt.grid()
            plt.savefig(out_dir/Path(f"infected_over_time_gamma_{self.gamma}.pdf"))

            plt.figure()
            plt.title(f"Total Population Per Node Over Time\nGamma = {self.gamma}")
            for n, node in enumerate(self.nodes):
                plt.plot(range(len(node.Pop_history)), node.Pop_history, label=f"Node {n}")
            #plt.legend()
            #plt.axis([0,config.steps+1,0,400])
            plt.xlabel("Time []")
            plt.ylabel("Quantity")
            plt.grid()
            plt.savefig(out_dir/Path(f"population_over_time_gamma_{self.gamma}.pdf"))

            plt.figure()
            plt.title(f"Total Infected Over Time\nGamma = {self.gamma}")
            plt.plot(range(len(self.I_total_history)), self.I_total_history, label = "Infected")
            #plt.legend()
            #plt.axis([0,config.steps+1,0,400])
            plt.xlabel("Time []")
            plt.ylabel("Quantity")
            plt.ylim(bottom=-1)
            plt.grid()
            plt.savefig(out_dir/Path(f"Tot_infected_over_time_gamma_{self.gamma}.pdf"))
        else:
            print("No movement of individuals were recorded")

    def steady_state_reached(
        self,
        window=50,
        rel_tol=0.05,
        abs_tol=0.5
    ):
        if len(self.nodes[0].Pop_history) < 2 * window:
            return False, None

        # Steady-state estimate for each individual node
        node_means = np.array([
            np.mean(node.Pop_history[-window:])
            for node in self.nodes
        ])

        # Check convergence using degree groups
        for k in np.unique(self.actual_degrees):

            nodes_k = [
                node for node in self.nodes
                if node.degree == k
            ]

            old_mean = np.mean([
                np.mean(node.Pop_history[-2*window:-window])
                for node in nodes_k
            ])

            new_mean = np.mean([
                np.mean(node.Pop_history[-window:])
                for node in nodes_k
            ])

            threshold = abs_tol + rel_tol * new_mean

            if abs(new_mean - old_mean) > threshold:
                return False, None
            

        return True, node_means

    def infected_steady_state_reached(self, window=50, I_rel_tol=0.02, I_abs_tol=1.0, step = 0):

        if len(self.nodes[0].Pop_history) < 2 * window:
            return False, None

        if step%5 == 1:
            return False, None

        # Individual node steady-state estimates
        node_means = np.array([
            np.mean(node.Pop_history[-window:])
            for node in self.nodes
        ])

        # Check total infected population
        old_I = np.mean(
            self.I_total_history[-2*window:-window]
        )

        new_I = np.mean(
            self.I_total_history[-window:]
        )

        I_threshold = I_abs_tol + I_rel_tol * new_I

        if abs(new_I - old_I) > I_threshold:
            return False, None

        return True, node_means



def plot_SFNs_degrees(SFNs,fig_dir):
    """
    Make a plot of degree distributions for the nodes in a list of Scale-Free Networks (SFNs)
    Figures are saved to fig_dir
    """
    
    plt.figure()
    for sfn in SFNs:
        degrees = sfn.target_degrees
        values, counts = np.unique(degrees,return_counts=True)
        plt.plot(values,counts, label=f"gamma={sfn.gamma}")
    plt.title(f"Degree Distributions for SFNs with \nVarying Power Law Gammas")
    plt.xlabel("Degree")
    plt.ylabel("# Occurrances")
    plt.grid()
    plt.legend()

    plt.savefig(fig_dir)

def plot_rho_vs_k(sfn, out_dir):
    unique_k = np.unique(sfn.actual_degrees)

    rho_measured = []
    rho_theory = []

    for k in unique_k:
        populations = [sfn.steady_state[i] for i, node in enumerate(sfn.nodes) if node.degree==k]
        rho_measured.append(np.mean(populations))

        rho_theory.append(k/sfn.avg_network_degree*sfn.avg_node_population)

    plt.figure()
    plt.title("Average populations per node depending on degree")
    plt.scatter(unique_k, rho_measured, label="Measured avg. populations")
    plt.plot(unique_k, rho_theory, label="Theoretical avg. populations")
    plt.xlabel("Degree k")
    plt.ylabel("Population")
    plt.legend()
    plt.grid()
    plt.savefig(out_dir/Path(f"pop_vs_degree_gamma_{sfn.gamma}.pdf"))




class State:
    INFECTED = 1
    SUSCEPTIBLE = 0


def susceptible_to_infected(mask_susceptible, N,beta, dt, rng):

    infected = N - np.sum(mask_susceptible)
    prob = 1 - (1 - beta*dt)**infected
    new_infections = mask_susceptible * prob > rng.random(N)

    return new_infections

def infected_to_susceptible(mask_infected,  N, mu, dt, rng):

    prob = mu*dt
    new_susceptible = mask_infected * prob > rng.random(N)

    return new_susceptible


#Stochastic infection
def stochastic_state(SFN: SFN = None , beta = None, mu = 0, timestep = 0):
    #STATE CHANGER
    rng = np.random.default_rng()
    for node in SFN.nodes:
        N = node.I + node.S
        individuals = np.full(N, State.SUSCEPTIBLE)
        infected_idx = rng.choice(N, node.I, replace = False)
        individuals[infected_idx] = State.INFECTED

        mask_susceptible = individuals == State.SUSCEPTIBLE
        mask_infected = individuals == State.INFECTED
        if beta == None:
            new_infected = susceptible_to_infected(mask_susceptible, N, node.beta, timestep, rng)
        else:
            new_infected = susceptible_to_infected(mask_susceptible, N, beta, timestep, rng)

        new_susceptible = infected_to_susceptible(mask_infected, N, mu, timestep, rng)

        individuals[new_infected] = State.INFECTED
        individuals[new_susceptible] = State.SUSCEPTIBLE

        nmb_susceptible = np.sum(individuals == State.SUSCEPTIBLE)
        nmb_infected = np.sum(individuals == State.INFECTED)

        node.S = nmb_susceptible
        node.I = nmb_infected
    return 

def plot_rel_infect_vs_init_rho(sfn_list, out_dir):
    plt.figure()
    plt.title(f"Steady State Infected Proportion \nDependency on Rho")
    V = []
    rho = []
    I_pop_rels = []
    for sfn in sfn_list:
        V.append(sfn.V)
        rho.append(sfn.avg_node_population)
        I_tot = 0

        for node in sfn.nodes:
            I_tot += node.I_history[-1] 

        I_pop_rels.append((I_tot/sfn.V)/sfn.avg_node_population) 

    # Group by V for plotting
    for v in np.unique(V):
        rho_group = [rho[i] for i in range(len(V)) if v==V[i]]
        I_group = [I_pop_rels[i] for i in range(len(V)) if v==V[i]]

        plt.plot(rho_group, I_group, label=f"V = {v}")
    out_file = out_dir / Path(f"rel_infect_vs_init_rho_gamma_{sfn_list[0].gamma}.pdf")
    plt.legend()
    plt.grid()
    plt.xlabel("Average Initial Population Per Node")
    plt.ylabel("Proportion of Total Population Which is Infected")

    plt.savefig(out_file)

    
        