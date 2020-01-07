import subprocess
import os
import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

args = [
    "./solution-step2",                         # Compiled C executable
    "0.01",                                     # tPlotDelta
    "10",                                       # tFinal
    "1e-8",                                     # timeStepSize
    "0", "0", "0",    "0", "0", "0",    "4",    # [x, y, z], [v_x, v_y, v_z], m
    "3", "0", "0",    "0", "0", "0",    "5",    # [x, y, z], [v_x, v_y, v_z], m
    "3", "4", "0",    "0", "0", "0",    "3"     # [x, y, z], [v_x, v_y, v_z], m
]

# 7e-6 is the largest time step that yields collisions.

# How does the position of the points depend on the time step chosen, i.e. can you observe that the position
# becomes more accurate (i.e. converges) if you use a finer time step size? Track the collision points for
# various time steps.

# Compute the distances between any two collision points for different time step sizes. I.e. pick two random collision points
# and compute the distance between them. If I have close to... 300 points, that gives me a ridiculous number of points!

# Can you derive the converge order of the method experimentally using the collision points?

# https://stackoverflow.com/questions/477486/how-to-use-a-decimal-range-step-value
mantissa = np.linspace(1, 9.9, 90)

def run():

    with open("solution-step2.csv", "w") as file:

        file.write("step,x,y,z\n")

        for e in range(7, 10):

            for m in mantissa:

                step = "{}e-{}".format(m, e)

                print("Testing {}... ".format(step), end="")

                start = time.time()

                args[3] = step

                process = subprocess.Popen(args, stdout=subprocess.PIPE) # https://stackoverflow.com/questions/6657690/python-getoutput-equivalent-in-subprocess
                out = process.communicate()[0].decode("utf-8").split("\n")

                if len(out) >= 2:

                    coordinates = out[-2]

                    message = coordinates

                    file.write(step + "," + coordinates.replace(" ", "") + "\n")

                else:

                    message = "NO COLLISION"

                print("{} ({:.2f}s)".format(message, time.time() - start))

def plot(num_points=100, type=None):

    step_differences = []
    collision_differences = []

    min_step_difference = math.inf
    max_step_difference = -math.inf

    min_collision_difference = math.inf
    max_collision_difference = -math.inf

    df = pd.read_csv("solution-step2.csv", dtype={"step": np.float64, "x": np.float64, "y": np.float64, "z": np.float64})

    for i in range(num_points):

        sample = df.sample(n=2, random_state=i) # Set random_state for reproducibility

        a = sample.iloc[0]
        b = sample.iloc[1]

        step_difference = abs(a["step"] - b["step"])

        collision_difference = abs(math.sqrt(
            (a["x"] - b["x"]) ** 2 +
            (a["y"] - b["y"]) ** 2 +
            (a["z"] - b["z"]) ** 2
        ))

        if type == "log":

            step_difference = math.log(step_difference)
            collision_difference = math.log(collision_difference)

        step_differences.append(step_difference)
        collision_differences.append(collision_difference)

        if step_difference < min_step_difference:

            min_step_difference = step_difference

        if step_difference > max_step_difference:

            max_step_difference = step_difference

        if collision_difference < min_collision_difference:

            min_collision_difference = collision_difference

        if collision_difference > max_collision_difference:

            max_collision_difference = collision_difference

    plt.scatter(step_differences, collision_differences, alpha=0.5)

    plt.xlim([min_step_difference, max_step_difference])
    plt.ylim([min_collision_difference, max_collision_difference])


    plt.title("Collision difference against step difference \n"
              "for {} randomly chosen pairs of collisions".format(num_points))
    plt.xlabel('Step difference')
    plt.ylabel('Collision difference')

    if type == "log":

        plt.title("Log-log collision difference against step difference \n"
                  "for {} randomly chosen pairs of collisions".format(num_points))
        plt.xlabel('Log step difference')
        plt.ylabel('Log collision difference')


    plt.tight_layout()

    plt.show()

if __name__ == "__main__":

    if not os.path.exists("solution-step2"):        # Compile solution-step2 if it doesn't exist

        subprocess.call(["g++", "-O3", "--std=c++11", "solution-step2.c", "-o", "solution-step2"])

    if not os.path.exists("solution-step2.csv"):    # Run solution-step2 if no output exists

        run()

    plot(type="log", num_points=300)
