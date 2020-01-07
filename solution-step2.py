import subprocess
import numpy as np

args = [
    "solution-step2"                            # Compiled C executable
    "0.01",                                     # tPlotDelta
    "10",                                       # tFinal
    "1e-8",                                     # timeStepSize
    "0", "0", "0",    "0", "0", "0",    "4",    # [x, y, z], [v_x, v_y, v_z], m
    "3", "0", "0",    "0", "0", "0",    "5",
    "3", "4", "0",    "0", "0", "0",    "3"
]

# 7e-6 is the largest time step that yields collisions (I've been told)

# So what does this file need to do? In short, we need to iterate over a bunch of different time steps,
# compute and tabulate and the results, and then go from there.

# How does the position of the points depend on the time step chosen, i.e. can you observe that the position
# becomes more accurate (i.e. converges) if you use a finer time step size? Track the collision points for
# various time steps. Compute the distances between any two collision points for different time step sizes.

# Does that mean 'compute the distance between two locations where collisions have occurred'? Or
# 'compute the distance between two points that are about to collide'?

# Can you derive the converge order of the method experimentally using the collision points?

# https://stackoverflow.com/questions/477486/how-to-use-a-decimal-range-step-value
mantissa = np.linspace(1, 9.9, 90)

# Ah, I see. Peter generates both at once. I, on the other hand, will have two

def run():

    with open("output.csv", "w") as file:

        file.write("step,x,y,z\n")

        for e in range(7, 10):

            for m in mantissa:

                time = "{}e-{}".format(m, e)

                print("Testing {}".format(time))

                file.write(time + ",")

                args[3] = time

                process = subprocess.Popen(args, stdout=subprocess.PIPE) # https://stackoverflow.com/questions/6657690/python-getoutput-equivalent-in-subprocess
                out = process.communicate()

                print(out)

def plot():

    pass


# if not os.path.exists('foo'):
#
#     subprocess.call(["gcc", "foo.c", "-ofoo", "-std=c99", '-w', '-Ofast'])
#
# subprocess.call(["./foo"], stdin = sys.stdin)

if __name__ == "__main__":

    # if not os.path.exists("output.csv"):

    run()

    plot()

