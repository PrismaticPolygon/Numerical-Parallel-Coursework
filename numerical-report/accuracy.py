import subprocess
import os
import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

rc("text", usetex="True")   # Use LaTeX text rendering

args = [
    "./numerical-report",                       # Compiled C executable
    "0.01",                                     # tPlotDelta
    "3.0",                                      # tFinal
    "1e-8",                                     # timeStepSize
    "0", "0", "0",    "0", "0", "0",    "4",    # [x, y, z], [v_x, v_y, v_z], m
    "3", "0", "0",    "0", "0", "0",    "5",    # [x, y, z], [v_x, v_y, v_z], m
    "3", "4", "0",    "0", "0", "0",    "3"     # [x, y, z], [v_x, v_y, v_z], m
]

# Collision should occur around 2.88s.
# 7.3e-6 is the largest time step that yields collisions.

def run():

    # https://stackoverflow.com/questions/477486/how-to-use-a-decimal-range-step-value
    mantissa = np.linspace(1, 9.9, 90)

    with open("accuracy.csv", "w") as file:

        file.write("step,x_1,y_1,z_1,t_1,x_2,y_2,z_2,t_2\n")

        for e in range(6, 10):

            for m in mantissa:

                step = "{}e-{}".format(m, e)

                print("Testing {}... ".format(step), end="")

                start = time.time()

                args[3] = step

                process = subprocess.Popen(args, stdout=subprocess.PIPE) # https://stackoverflow.com/questions/6657690/python-getoutput-equivalent-in-subprocess
                out = process.communicate()[0].decode("utf-8").split("\n")

                if len(out) >= 2:

                    coordinates = ",".join(out[:-1])

                    message = coordinates

                    file.write(step + "," + message + "\n")

                else:

                    message = "NO COLLISION"

                print("{} ({:.2f}s)".format(message, time.time() - start))

def distance(row):

    return math.sqrt(
        (row["x_1"] - row["x_2"]) ** 2 +
        (row["y_1"] - row["y_2"]) ** 2 +
        (row["z_1"] - row["z_2"]) ** 2
    )

def plot():

    df = pd.read_csv("accuracy.csv")

    plt.figure(figsize=(6, 3))  # Default is (6, 4)

    df["distance"] = df.apply(distance, axis=1)
    df["step"] = np.log10(df["step"])

    plt.scatter(df["step"], df["distance"], alpha=0.5)

    plt.xlim([df["step"].min(), df["step"].max()])
    plt.ylim([df["distance"].min(), df["distance"].max()])

    plt.title(r"Observing convergence")
    plt.xlabel(r'$\log_{10} h$')
    plt.ylabel(r'Collision distance $|d|$')

    plt.tight_layout(pad=0.8)

    plt.savefig("accuracy.png")

    plt.show()


if __name__ == "__main__":

    subprocess.call(["g++", "-O3", "--std=c++11", "report.c", "-o", "report"])

    if not os.path.exists("accuracy.csv"):

        run()

    plot()