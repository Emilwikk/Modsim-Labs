from pathlib import Path

import config
import time
import utils

start = time.perf_counter()
simNr = 1
sims = len(config.gammas_dan)
for gamma in config.gammas_dan:

    print(f"Running simulation [{simNr}/{sims}]")
    mainPy_dir = Path(__file__).resolve().parent
    fig_dir = mainPy_dir / Path(f"../results/diffusion/gamma_{gamma}/")
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Create and Populate Network 
    sfn = utils.SFN(config.V_dan,gamma)
    sfn.populate(S=config.S_dan, I=config.I_dan)

    steady_state_reached = False
    step = 0
    while steady_state_reached==False:
        print(f"---Step: {step}---", end="\r", flush=True)
        sfn.move_individuals()
        steady_state_reached, steady_state = sfn.steady_state_reached(wdw_size=10, tol=0.05, N = config.S_dan+config.I_dan, V=config.V_dan)
        step += 1
    print(f"Steady State Reached After {step} Steps")
    sfn.plot_populations(out_dir=fig_dir)
    sfn.steady_state = steady_state

    utils.plot_rho_vs_k(sfn, out_dir=fig_dir)

    print(f"Simulation [{simNr}/{sims}] complete\n")
    simNr += 1

end = time.perf_counter()
print(f"Total Simulation Time: {end-start:.1f}s\n")

