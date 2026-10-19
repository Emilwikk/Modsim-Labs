from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from MD_plotting import read_efile, calc_E_stats

this_dir = Path(__file__).resolve().parent
outDir = this_dir / Path("../results/lang/")
outDir.mkdir(parents=True, exist_ok=True)

# Initialize list of simulation objects
sims_list = []

# Loop through all files in efile directory
efiles = this_dir / Path("../lang/efile")
for efile in efiles.iterdir():
    # Only include if correct parameters
    if efile.is_file() and efile.name.startswith("0064_r0.600_T1.000_alpha") and efile.name.endswith("dt010"):
        sim = calc_E_stats(efilePath=efile)
        sims_list.append(sim)


# Plot E vs t for different alpha
for sim in sims_list:
    t = np.linspace(0,sim.params["dt"]*len(sim.E), len(sim.E))
    plt.figure()
    plt.title(f"Average energy over time for alpha={sim.params["alpha"]}")
    plt.xlabel("Time")
    plt.ylabel("Energy")
    plt.grid()
    plt.plot(t,sim.E)
    outFile = outDir / Path(f"E_vs_t_alpha{sim.params["alpha"]}.pdf")
    plt.savefig(outFile)

# Sort list by alpha
sims_list.sort(key=lambda sim: sim.params["alpha"])
std_E = []
alpha = []
[std_E.append(sim.std_E) for sim in sims_list]
[alpha.append(sim.params["alpha"]) for sim in sims_list]

outFile = outDir / Path("std_E_vs_alpha.pdf")
plt.figure()
plt.plot(alpha, std_E)
plt.xlabel("Damping coefficient alpha")
plt.ylabel("Energy")
plt.title("Energy fluctuation size [std(E)]")
plt.grid()
plt.savefig(outFile)








