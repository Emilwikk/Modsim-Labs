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
    plt.savefig(fig_dir / f"sim1_mu{mu:.1e}_beta{beta:.1e}.png")
    if show_plots == 1:
        plt.show()

def powerlaw_adjacency_matrix(V, gamma):
    """
    Create an adjacency matrix A for V number of nodes using power law: P(k_i)=k^gamma.
    Returns adjacency matrix 'A' and degrees of each node in vector 'degrees'    
    """

    A = np.zeros((V,V),dtype=int)
    degrees = []
    k_min = 2
    k_max = np.sqrt(V)

    for i in range(V):
        # Random float between 0 and 1
        r = np.random.rand()
        ki = np.floor(((k_max**(gamma+1)-k_min**(gamma+1))*r+k_min**(gamma+1))**(1/(gamma+1)))
        degrees.append(int(ki))

    weird_list = []
    for i,deg in enumerate(degrees):
        [weird_list.append(i) for j in range(deg)]

    while weird_list and len(set(weird_list))>1:
        i,j = random.sample(weird_list,2)
        if i!=j and A[i,j]==0: # A[i,j]==0 because if not, i and j will be removed from weird_list even though a new connection has not been made    
            A[i,j] = 1
            A[j,i] = 1
            weird_list.remove(i)
            weird_list.remove(j)

    return A, degrees

def fit_gamma(degrees):

    degree, count = np.unique(degrees,return_counts=True)
    slope,inter = np.polyfit(np.log(degree),np.log(count),deg=1)
    fitted_gamma = slope
    return fitted_gamma

class Node:
    def __init__(self,S,I):
        self.S = S
        self.I = I

class SFN:
    def __init__(self, V, gamma):
        self.V = V
        self.gamma = gamma
        A, degrees = powerlaw_adjacency_matrix(self.V,self.gamma)
        self.A = A
        self.degrees = degrees
        self.fitted_gamma = fit_gamma(self.degrees)

    def populate(self, S, I):
        self.nodes = []
        probabilities = np.ones(self.V)/self.V

        S_nodes = np.random.multinomial(S,probabilities)
        I_nodes = np.random.multinomial(I,probabilities)

        for i in range(self.V):
            self.nodes.append(Node(S_nodes[i],I_nodes[i]))


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

