import time
import subprocess
import os
import math
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

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

def plot():

    df = pd.read_csv("convergence.csv")

    distances = []
    timesteps = []

    a = df.iloc[0] # Largest timestep: 7.3e-6.

    for i in range(1, len(df)):

        b = df.iloc[i]

        distance = difference(a, b)
        h = b["step"]

        distances.append(distance)
        timesteps.append(h)

    timesteps = np.array(timesteps)
    distances = np.array(distances)

    plt.scatter(timesteps, distances, label="Data")

    params, params_covariance = optimize.curve_fit(test_func, timesteps, distances)

    print("C:", params[0])
    print("p:", params[1])

    plt.plot(timesteps, test_func(timesteps, params[0], params[1]), label='Fitted function')

    plt.xlim([min(timesteps), max(timesteps)])
    plt.ylim([min(distances), max(distances)])

    plt.legend(loc='best')

    plt.title("Order of convergence")
    plt.xlabel('Timestep h')
    plt.ylabel('Distance between second collision points |d|')

    plt.tight_layout()

    plt.savefig("convergence.png")

    plt.show()


if __name__ == "__main__":

    subprocess.call(["g++", "-O3", "--std=c++11", "numerical-report.c", "-o", "numerical-report"])

    if not os.path.exists("convergence.csv"):    # Run convergence if no output exists

        run()

    plot()
