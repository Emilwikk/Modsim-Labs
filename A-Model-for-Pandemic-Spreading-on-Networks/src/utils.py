from matplotlib import pyplot as plt

def reaction_diffusion_IS(t, y, mu, beta):
    I, S = y

    dIdt = -mu*I + beta*I*S
    dSdt = mu*I - beta*I*S
    return [dIdt, dSdt]

def plot_sim1(t, y, mu, beta,  spd, fig_dir, show_plots):
    plt.figure()
    plt.title(f"Disease Spreading Simulation \nwith mu={mu} and beta={beta}")
    plt.plot(t/spd,y[0]/(y[0]+y[1]),label='# Infected')
    plt.plot(t/spd,y[1]/(y[0]+y[1]), label='# Susceptible')
    plt.legend()
    plt.grid()
    plt.xlabel("Time [Days]")
    plt.ylabel("Relative Population")
    plt.savefig(fig_dir / f"sim1_mu{mu}_beta{beta}.png")
    if show_plots == 1:
        plt.show()    
