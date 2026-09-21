# Global
seconds_per_day = 24*60*60
SHOW_PLOTS = 0

# Simulation parameters

mus = [1e-6, 2e-6, 4e-6]    # Recovery rate [s^-1]
betas = [4e-6, 1e-6, 6e-6]  # Infection rate [s^-1]

initial_infected = 2
initial_susceptible = 30

days_simulated = 5