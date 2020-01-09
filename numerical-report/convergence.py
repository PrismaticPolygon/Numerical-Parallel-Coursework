import time
import subprocess
import os

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

    h = 0.0000073

    with open("convergence.csv", "w") as file:

        file.write("step,x_1,y_1,z_1,t_1,x_2,y_2,z_2,t_2\n")

        for _ in range(13):

            step = str(h)

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

            h = h / 2

if __name__ == "__main__":

    subprocess.call(["g++", "-O3", "--std=c++11", "numerical-report.c", "-o", "numerical-report"])

    if not os.path.exists("convergence.csv"):    # Run solution-step2 if no output exists

        run()
