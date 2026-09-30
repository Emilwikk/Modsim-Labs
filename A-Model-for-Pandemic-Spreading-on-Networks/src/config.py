## Global ##
seconds_per_day = 24*60*60
SHOW_PLOTS = 0

## Reaction Diffusion SIS ##

mus =   [4   , 4    , 4   ]    # Recovery rate [per day]
betas = [8/32, 10/32, 6/32]      # Infection rate [per day]

initial_infected = 2
initial_susceptible = 30

days_simulated = 5


## Scale Free Networks ##

# Number of nodes
V = 100

# Power law parameter
pl_gammas = [-1.5,-2.5,-3.5]

## Diffusion Along Networks

V_dan = 1000
gamma_dan = -2.5
gammas_dan = [-1.5, -2.5, -3.5]
S_dan = 100000000
I_dan = 30000

## Combining it all ##

mu5 = 2
beta5 = 1
gamma5 = -1.5

V_list = [100, 500, 1000]
rho_min = 0.1
rho_max = 2
n_rhos = 8

init_proportion_infected = 0.1
