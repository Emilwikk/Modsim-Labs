from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

this_dir = Path(__file__).resolve().parent
print(this_dir)


def read_conf(inPath: Path):

    # Read parameters from top comment of file
    params = {}
    with open(inPath, "r") as f:
       params = list(map(float, f.readline().lstrip("#").split()))

    dt = int(inPath.name.split("_dt")[-1]) / 1000
    
    df = pd.read_csv(inPath,
                sep=r"\s+", # Interpret every space (or multiple spaces) as separator
                comment="#",
                header=None,
                names=["x","y","vx","vy"])

    return df, params

inPath = this_dir / Path("../lang/conf/0064_r0.500_T1.000_alpha0.10_dt010")

df, params = read_conf(inPath)
print(df.head())
print(params)
plt.figure()
plt.plot(df["x"],np.linspace(0, 0.1, len(df["x"])))
plt.show()



