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

    return A

