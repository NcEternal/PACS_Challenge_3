import pandas as pd
import matplotlib.pyplot as plt

import os

stats = ["calculation_time_(ms)", "L2_error"]

for file in sorted(os.listdir("./data")):
    if file.endswith("procs.csv"):
        df = pd.read_csv("./data/" + file)
        processors = file[file.find("procs.csv") - 1] + " processors"
        for stat in stats:
            plt.figure(stat)
            plt.plot(df["n"], df[stat], label = processors)

for stat in stats:
    plt.figure(stat)
    plt.xlabel("n")
    title = " ".join([w.capitalize() for w in stat.split("_")])
    plt.legend()
    plt.title(title)
    title = "_".join(title.replace("(", "").replace(")", "").split())
    plt.savefig("./data/" + "_".join(title.split()) + ".png")
    