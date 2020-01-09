import subprocess
from scipy import optimize
import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# How does the position of the points depend on the time step chosen, i.e. can you observe that the position
# becomes more accurate (i.e. converges) if you use a finer time step size? Track the collision points for
# various time steps.

# Compute the distances between any two collision points for different time step sizes. I.e. pick two random collision points
# and compute the distance between them. If I have close to... 300 points, that gives me a ridiculous number of points!

# Can you derive the converge order of the method experimentally using the collision points?

# https://stackoverflow.com/questions/477486/how-to-use-a-decimal-range-step-value
mantissa = np.linspace(1, 9.9, 90)

def run():

    with open("numerical-report.csv", "w") as file:

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

def difference(a, b):

    return math.sqrt(
        (a["x_2"] - b["x_2"]) ** 2 +
        (a["y_2"] - b["y_2"]) ** 2 +
        (a["z_2"] - b["z_2"]) ** 2
    )

def test_func(h, C, p):

    return (C * (h ** p)) - (C * ((h/2) ** p))

# Rate of convergence is also known as the global order of accuracy. Describes the decrease in the maximum
# numerical error one can expect for a given decrease in time step h in the limit h -> 0.
# A numerical method has global order of accuracy p if the MAXIMUM
# v ** n - u ** n, where u  is the
# This maximum is spooking me.

# Equivalent, but different forms of notation, methinks.
# Convergent if the maximum error tends to 0 as the timestep tends to 0
# As we shrink the time step smaller and smaller, the largest absolute error between the numerical
# and exact solution will also get smaller and smaller.
# Not quite useful. Oh well; I've done my best.
# D'oh! We have the change the base. Silly boy.
# So it does have to be divisible.
# Bollocks. Well - I've just divided it by 10 each time. Right?

def convergence():

    df = pd.read_csv("numerical-report.csv")

    distances = []
    timesteps = []

    a = df.iloc[333] # Smallest timestep: 9.9e-9

    for i in range(len(df)):

        if i != 333:

            b = df.iloc[i]

            distance = difference(a, b)
            h = b["step"]

            distances.append(distance)
            timesteps.append(h)

    timesteps = np.array(timesteps)
    distances = np.array(distances)

    plt.scatter(timesteps, distances, label="Data")

    params, params_covariance = optimize.curve_fit(test_func, timesteps, distances, p0=[11000, 1.3])

    print("C:", params[0])
    print("p:", params[1])

    plt.plot(timesteps, test_func(timesteps, params[0], params[1]), label='Fitted function')

    plt.xlim([min(timesteps), max(timesteps)])
    plt.ylim([min(distances), max(distances)])

    plt.legend(loc='best')

    plt.title("Order of convergence")
    plt.xlabel('Timestep h')
    plt.ylabel('Distance between second collision point |d|')

    plt.tight_layout()

    plt.savefig("convergence.png")

    plt.show()

def log_plot():

    df = pd.read_csv("numerical-report.csv")

    df["distance"] = np.log10(df.apply(distance, axis=1))
    df["step"] = np.log10(df["step"])

    plt.scatter(df["step"], df["distance"], alpha=0.5)

    plt.xlim([df["step"].min(), df["step"].max()])
    plt.ylim([df["distance"].min(), df["distance"].max()])

    plt.title("Log10-log10 time step against collision distance")
    plt.xlabel('Log10 time step')
    plt.ylabel('Log10 collision distance')

    plt.tight_layout()

    plt.savefig("log-log.png")

    plt.show()


def plot():

    df = pd.read_csv("numerical-report.csv")

    df["distance"] = df.apply(distance, axis=1)

    plt.scatter(df["step"], df["distance"], alpha=0.5)

    plt.xlim([df["step"].min(), df["step"].max()])
    plt.ylim([df["distance"].min(), df["distance"].max()])

    plt.title("Time step against collision distance")
    plt.xlabel('Time step')
    plt.ylabel('Collision distance')

    plt.tight_layout()

    plt.savefig("linear.png")

    plt.show()



if __name__ == "__main__":

    # subprocess.call(["g++", "-O3", "--std=c++11", "numerical-report.c", "-o", "numerical-report"])
    #
    # if not os.path.exists("numerical-report.csv"):    # Run solution-step2 if no output exists
    #
    #     run()
    #
    # plot()
    #
    # log_plot()

    convergence()