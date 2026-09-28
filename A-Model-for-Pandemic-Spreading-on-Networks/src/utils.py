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

def plot_sim1(t, y, mu, beta,  spd, fig_dir, show_plots):
    plt.figure()
    plt.title(f"Disease Spreading Simulation \nwith mu={mu:.2e} and beta={beta:.2e}")
    plt.plot(t/spd,y[0]/(y[0]+y[1]),label='# Infected')
    plt.plot(t/spd,y[1]/(y[0]+y[1]), label='# Susceptible')
    plt.legend()
    plt.grid()
    plt.xlabel("Time [Days]")
    plt.ylabel("Relative Population")
    plt.savefig(fig_dir / f"sim1_mu{mu:.1e}_beta{beta:.1e}.pdf")
    if show_plots == 1:
        plt.show()

def powerlaw_adjacency_matrix(V, gamma):
    """
    Create an adjacency matrix A for V number of nodes using power law: P(k_i)=k^gamma.
    Returns adjacency matrix 'A' and degrees of each node in vector 'degrees'    
    """

    A = np.zeros((V,V),dtype=int)
    target_degrees = []
    k_min = 2
    k_max = np.sqrt(V)

    for i in range(V):
        # Random float between 0 and 1
        r = np.random.rand()
        ki = np.floor(((k_max**(gamma+1)-k_min**(gamma+1))*r+k_min**(gamma+1))**(1/(gamma+1)))
        target_degrees.append(int(ki))

    weird_list = []
    for i,deg in enumerate(target_degrees):
        [weird_list.append(i) for j in range(deg)]

    while weird_list and len(set(weird_list))>1:
        unique_nodes = set(weird_list)

        valid_pairs = [(i,j)
         for i in unique_nodes
         for j in unique_nodes
         if i<j and A[i,j]==0]

        if not valid_pairs:
            break

        i, j = random.choice(valid_pairs)

        A[i,j] = 1
        A[j,i] = 1

        weird_list.remove(i)
        weird_list.remove(j)

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

        self.avg_node_population = (S+I)/self.V

        

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

        # Add changes to actual nodes
        for n, node in enumerate(self.nodes):
            node.S = S_move[n]
            node.I = I_move[n]
            node.Pop = node.S+node.I

            # Save to history
            node.S_history.append(node.S)
            node.I_history.append(node.I)
            node.Pop_history.append(node.Pop)

    def plot_populations(self, out_dir:Path):
        if len(self.nodes[0].S_history)>1:
            plt.figure()
            plt.title("Number of Susceptible Individuals Per Node Over Time")
            for n, node in enumerate(self.nodes):
                plt.plot(range(len(node.S_history)), node.S_history, label=f"Node {n}")
            plt.legend()
            #plt.axis([0,config.steps+1,0,400])
            plt.xlabel("Time []")
            plt.ylabel("Quantity")
            plt.grid()
            plt.savefig(out_dir/"susceptible_over_time.pdf")

            plt.figure()
            plt.title("Number of Infected Individuals Per Node Over Time")
            for n, node in enumerate(self.nodes):
                plt.plot(range(len(node.I_history)), node.I_history, label=f"Node {n}")
            plt.legend()
            #plt.axis([0,config.steps+1,0,400])
            plt.xlabel("Time []")
            plt.ylabel("Quantity")
            plt.grid()
            plt.savefig(out_dir/"infected_over_time.pdf")

            plt.figure()
            plt.title("Total Population Per Node Over Time")
            for n, node in enumerate(self.nodes):
                plt.plot(range(len(node.Pop_history)), node.Pop_history, label=f"Node {n}")
            plt.legend()
            #plt.axis([0,config.steps+1,0,400])
            plt.xlabel("Time []")
            plt.ylabel("Quantity")
            plt.grid()
            plt.savefig(out_dir/"population_over_time.pdf")
        else:
            print("No movement of individuals were recorded")

    def steady_state_reached(self, wdw_size=10, tol=0.01, N = 0, V = 0):
        means = []

        for node in self.nodes:

            if len(node.Pop_history) < 2 * wdw_size:
                return False, None
                
            old_mean = np.mean(node.Pop_history[-2*wdw_size:-wdw_size])
            new_mean = np.mean(node.Pop_history[-wdw_size:])
            if np.abs(old_mean-new_mean)>(tol*node.degree*N)/V:
                return False, new_mean
            
            means.append(new_mean)

        return True, means



def plot_SFNs_degrees(SFNs,fig_dir):
    """
    Make a plot of degree distributions for the nodes in a list of Scale-Free Networks (SFNs)
    Figures are saved to fig_dir
    """
    
    plt.figure()
    for sfn in SFNs:
        degrees = sfn.degrees
        values, counts = np.unique(degrees,return_counts=True)
        plt.plot(values,counts, label=f"gamma={sfn.gamma}")
    plt.title(f"Degree Distributions for SFNs with \nVarying Power Law Gammas")
    plt.xlabel("Degree")
    plt.ylabel("# Occurrances")
    plt.grid()
    plt.legend()

    plt.savefig(fig_dir)

def plot_rho_vs_k(sfn):
    rho_k = []
    k = []

    for i, node in enumerate(sfn.nodes):
        rho_k.append(node.degree/sfn.avg_network_degree*sfn.steady_state[i])
        k.append(node.degree)

    plt.figure()
    plt.scatter(k,rho_k)
    plt.show()


