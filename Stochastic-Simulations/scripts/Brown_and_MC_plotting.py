from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

from MD_plotting import read_efile, calc_E_stats

this_dir = Path(__file__).resolve().parent

brown_efiles = this_dir / Path("../brown/efile")

sims_list = []
# Read brownian files
for file in brown_efiles.iterdir():
    if file.is_file() and file.name.startswith("0064_r0.300_T1.000_alpha10.00_dt00"):
        sims_list.append(calc_E_stats(file))

# Read MC efile
MC_efile = this_dir / Path("../mc/efile/0064_r0.300_T1.000")
sims_list.append(calc_E_stats(MC_efile))

# In order to plot MC sim at dt/alpha=0 as per exercise description
sims_list[-1].params["dt"] = 0
sims_list[-1].params["alpha"] = 1 # (avoid division by zero)

# Sort list of simulations by dt
sims_list.sort(key=lambda sim: sim.params["dt"])

dt_over_alpha = []
exp_E = []
std_E = []
for sim in sims_list:
    dt_over_alpha.append(sim.params["dt"] / sim.params["alpha"])
    exp_E.append(sim.exp_E)
    std_E.append(sim.std_E)

outDir = this_dir / Path("../results/mc_brown")
outDir.mkdir(exist_ok=True, parents=True)

plt.figure()
plt.title("Expected value of potential energy depending on \ntimestep and damping coefficient")
plt.xlabel("dt/alpha *[1e4]")
plt.ylabel("Energy")
plt.grid()
plt.plot(1e4*np.array(dt_over_alpha), exp_E, color="red")
plt.scatter(1e4*np.array(dt_over_alpha), exp_E, color="red")
plt.plot(1e4*np.array(dt_over_alpha), np.array(exp_E)+np.array(std_E), color="black", linestyle="dashed")
plt.plot(1e4*np.array(dt_over_alpha), np.array(exp_E)-np.array(std_E), color="black", linestyle="dashed")
plt.savefig(outDir / Path("E_plot_alt_and_probably_wrong.pdf"))

# Manually extracted values from terminal prints:
exp_pot_E = [-0.8383, -0.8324, -0.8299, -0.8237, -0.8232]
SE_pot_E = [0.001542, 0.001699, 0.001423, 0.001670, 0.001579] # standard error (not standard deviation)

plt.figure()
plt.title(f"Expected value of potential energy depending on \ntimestep and damping coefficient")
plt.xlabel("dt/alpha * 1e4")
plt.ylabel("Energy")
plt.grid()
plt.plot(1e4*np.array(dt_over_alpha), exp_pot_E, color="red", label="Mean potential energy")
plt.scatter(1e4*np.array(dt_over_alpha), exp_pot_E, color="red")
plt.plot(1e4*np.array(dt_over_alpha), np.array(exp_pot_E)+np.array(SE_pot_E), color="black", linestyle="dashed", label=r"Mean $\pm$ SE")
plt.plot(1e4*np.array(dt_over_alpha), np.array(exp_pot_E)-np.array(SE_pot_E), color="black", linestyle="dashed")
plt.tight_layout()
plt.legend()
plt.savefig(outDir / Path("E_plot_correct.pdf"))