import time
import subprocess
import os
import math
import pandas as pd
import numpy as np
from scipy import optimize
from decimal import Decimal
import matplotlib.pyplot as plt

from matplotlib import rc

# I could probably merge these two. Let's do that later, though.

rc("text", usetex="True")   # Use LaTeX text rendering

args = [
    "./numerical-report",                       # Compiled C executable
    "0.01",                                     # tPlotDelta
    "2.9",                                      # tFinal
    "1e-8",                                     # timeStepSize
    "0", "0", "0",    "0", "0", "0",    "4",    # [x, y, z], [v_x, v_y, v_z], m
    "3", "0", "0",    "0", "0", "0",    "5",    # [x, y, z], [v_x, v_y, v_z], m
    "3", "4", "0",    "0", "0", "0",    "3"     # [x, y, z], [v_x, v_y, v_z], m
]

def run():

    h = 0.0000073   # Smallest timestep at which a collision occurs

    with open("convergence.csv", "w") as file:

        file.write("step,x_1,y_1,z_1,t_1,x_2,y_2,z_2,t_2\n")

        for _ in range(13):

            step = str(h)

            print("Testing {}... ".format(step), end="")

            start = time.time()

            args[3] = step

            process = subprocess.Popen(args, stdout=subprocess.PIPE) # https://stackoverflow.com/questions/6657690/python-getoutput-equivalent-in-subprocess
            out = process.communicate()[0].decode("utf-8").split("\n")

            coordinates = ",".join(out[:-1])
            file.write(step + "," + coordinates + "\n")

            print("{} ({:.2f}s)".format(coordinates, time.time() - start))

            h = h / 2

def difference(a, b):

    return math.sqrt(
        (a["x_2"] - b["x_2"]) ** 2 +
        (a["y_2"] - b["y_2"]) ** 2 +
        (a["z_2"] - b["z_2"]) ** 2
    )

def test_func(h, C, p):

    return (C * (h ** p)) - (C * ((h/2) ** p))

def format_func(value, tick_number):

    float_str = "{0:.2g}".format(value)

    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
    else:

        return float_str

def plot():

    df = pd.read_csv("convergence.csv")

    plt.figure(figsize=(6, 3))  # Default is (6, 4)

    distances = []
    timesteps = []

    for i in range(1, len(df)):

        y_h = df.iloc[i - 1]
        y_h_2 = df.iloc[i]

        distance = difference(y_h, y_h_2)
        h = y_h["step"]

        distances.append(distance)
        timesteps.append(h)

    timesteps = np.array(timesteps)
    distances = np.array(distances)

    plt.scatter(timesteps, distances, label="Data", alpha=0.5)

    params, params_covariance = optimize.curve_fit(test_func, timesteps, distances)

    print("C:", params[0])
    print("p:", params[1])

    plt.plot(timesteps, test_func(timesteps, params[0], params[1]), label='Fitted function', color="red")

    plt.xlim([min(timesteps), max(timesteps)])
    plt.ylim([min(distances), max(distances)])

    axes = plt.gca()

    axes.xaxis.set_major_formatter(plt.FuncFormatter(format_func))

    plt.legend(loc='best')

    plt.title(r"Calculating the order of convergence")
    plt.xlabel(r'Timestep $h$')
    plt.ylabel(r'Collision distance $|d|$')

    plt.tight_layout(pad=0.8)

    plt.savefig("convergence.png")

    plt.show()


if __name__ == "__main__":

    subprocess.call(["g++", "-O3", "--std=c++11", "report.c", "-o", "report"])

    #if not os.path.exists("convergence.csv"):    # Run convergence if no output exists
run()
plot()
