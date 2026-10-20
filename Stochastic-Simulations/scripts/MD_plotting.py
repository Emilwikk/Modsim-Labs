from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

this_dir = Path(__file__).resolve().parent


def read_conf(inPath: Path):

    # Read parameters from top comment of file
    with open(inPath, "r") as f:
       values = list(map(float, f.readline().lstrip("#").split()))

    params = dict(zip(["N", "unknown1", "rho", "T", "unknown2"], values))
    # Get dt from file name
    params["dt"] = int(inPath.name.split("_dt")[-1]) / 1000

    # Get alpha from file name
    params["alpha"] = float(inPath.name.split("alpha")[-1].split("_dt")[0])
    
    df = pd.read_csv(inPath,
                sep=r"\s+", # Interpret every space (or multiple spaces) as separator
                comment="#",
                header=None,
                names=["x","y","vx","vy"])

    return df, params


def read_efile(inPath: Path):

    # Read parameters from top comment of file
    #params = dict(zip(["N", "unknown1", "rho", "T", "unknown2"], values))
    params = {}

    # Get dt from file name
    if "_dt" in inPath.name:
        params["dt"] = int(inPath.name.split("_dt")[-1]) / 1000
    else:
        params["dt"] = None

    # Get alpha from file name
    if "alpha" in inPath.name:
        params["alpha"] = float(inPath.name.split("alpha")[-1].split("_dt")[0])
    else:
        params["alpha"] = None
        
    df = pd.read_csv(inPath,
                sep=r"\s+", # Interpret every space (or multiple spaces) as separator
                comment="#",
                header=None,
                names=["E"])

    return df, params    

class simulation:   
    def __init__(self, params):
        self.params = params

def calc_E_stats(efilePath: Path) -> simulation: 

    df, params = read_efile(efilePath)
    sim = simulation(params=params)
    sim.E = df
    exp_E = np.mean(df)
    std_E = np.std(df)
    sim.exp_E = exp_E
    sim.std_E = std_E

    return sim

# Initialize list of simulation objects
sims_list = []

# Loop through all files in efile directory
efiles = this_dir / Path("../lang/efile")
for efile in efiles.iterdir():
    # Only include if correct parameters
    if efile.is_file() and efile.name.startswith("0064_r0.500_T1.000_alpha0.00") and not efile.name.endswith("16"):
        sim = calc_E_stats(efilePath=efile)
        sims_list.append(sim)


def plot_E_over_dt(sims_list, outDir):        
    # Sort list by dt
    sims_list.sort(key=lambda sim: sim.params["dt"])

    # Make vectors for plot
    dt = []
    [dt.append(sim.params["dt"]) for sim in sims_list]
    exp_E = []
    [exp_E.append(sim.exp_E) for sim in sims_list]
    std_E = []
    [std_E.append(sim.std_E) for sim in sims_list]

    plt.figure()
    plt.title("Average total energy per particle")
    plt.plot(dt, exp_E, color="r")
    plt.plot(dt, np.array(exp_E)+np.array(std_E), linestyle="dashed", color="black")
    plt.plot(dt, np.array(exp_E)-np.array(std_E), linestyle="dashed", color="black")
    plt.xlabel("dt")
    plt.ylabel("Energy")
    plt.grid()

    # Make results directory if it does not exist and save figure
    outDir.mkdir(parents=True, exist_ok=True)
    outFile = outDir / Path("dt_vs_E.pdf")
    plt.savefig(outFile)

    plt.figure()
    plt.title("Standard deviation of total energy per particle")
    plt.plot(dt, std_E, color="r")
    plt.xlabel("dt")
    plt.ylabel("Energy")
    plt.grid()
    outFile = outDir / Path("dt_vs_std_E.pdf")
    plt.savefig(outFile)

outDir = this_dir / Path("../results/MD/")
plot_E_over_dt(sims_list=sims_list, outDir=outDir)

### 3.1.2 ###
efile_instability = this_dir / Path("../lang/efile/0064_r0.500_T1.000_alpha0.00_dt016")
sim = calc_E_stats(efilePath=efile_instability)
t = np.linspace(0,sim.params["dt"] * len(sim.E), len(sim.E))
plt.figure()
plt.plot(t, np.array(sim.E).flatten())
plt.ylim([-0.4, 1])
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title(f"Energy over time using dt={sim.params["dt"]}")
plt.grid()

outFile = this_dir / Path("../results/MD/unstable_simulation.pdf")
plt.savefig(outFile)