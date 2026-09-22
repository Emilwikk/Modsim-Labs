from matplotlib import pyplot as plt
import numpy as np

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

def sample_degrees_powerlaw(A, gamma):
    """
    Sample a degree k for each node i in adjacency matrix A.
    Adjacency matrix "A" must be numpy square matrix. 
    Returns list of degree for each node (length is the number of rows or columns of A). 
    """

    degrees = []
    V = len(A[:,0])
    k_min = 2
    k_max = np.sqrt(V)

    for i in range(V):
        # Random float between 0 and 1
        r = np.random.rand()
        ki = np.floor(((k_max**(gamma+1)-k_min**(gamma+1))*r+k_min**(gamma+1))**(1/(gamma+1)))
        degrees.append(int(ki))
    return degrees

