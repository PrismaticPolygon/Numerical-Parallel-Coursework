import subprocess
import os
import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

args = [
    "./numerical-report",                       # Compiled C executable
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

# We could simply compute the absolute value from the origin.
# We don't even need to, in fact.
# Deriving the convergence order of the method.
# An algorithm converges iff it is consistent and stable.
# An algorithm converges with order p if |y_h(t) - y(t)| <= C dot h^p.

def convergence():

    df = pd.read_csv("numerical-report.csv")

    # This ISN'T what we want. We want the distance between their collision points.

    # distances = list(df.apply(distance, axis=1))  # So this is the distances between collision points.
    # steps = list(df["step"])

    # Ah. We also need timestep, naturally.

    # print(distances)
    # print(steps)
    # Is it a question of too many points? No, that's never the case.

    errors = []
    hs = []

    # I'm not seeing any correlation whatsoever.
    # And
    # If the convergence order CAN be derived, it's eluding me.
    # This seems a good point to stop.
    # Actually, I'll have a look into number 5.

    # Yes.
    # nums = [1, 2, 4, 8, 16, 32, 64, 128, 256]

    # That's more interesting. In fact
    # Now we can do it with EVERY pair. That's 300(300 - 1) / 2. Definitely possible
    # Let's do it.

    # Interesting. Is this what I was expecting?
    # I'm not sure.
    # No. In his, every point is on the same line.
    # I feel I'm just making this up!
    # At least there's no negative values.
    # Am I even calculating y properly? Yeah, the only way that's meaningful anyway.

    # num = len(df)
    #
    # for i in range(num - 1):
    #
    #     a = df.iloc[i]
    #
    #     for j in range(i, num):
    #
    #         b = df.iloc[j]
    #
    #         error = difference(a, b)
    #         h = a["step"]
    #
    #         errors.append(error)
    #         hs.append(h)
    # That's not right either! It's DIVIDED.
    # For logarithmic purposes. So how will I get the HALF value? We start with the biggest step.

    # print(df.iloc[0]["step"])   # 1e-6, 5e-7, 2.5e-7, 1.25e7

    # Okay.

    a = df.iloc[0]
    # i = h / 2

    # I guess a straight line of 0 means it converges to 0...
    # Or do I choose different values of h?
    # I'm convinced that must be what I'm supposed to do.
    # When I did that I got my unholy other graph, however. It's possible that
    # it's not possible to
    # That's not right either. Every point we be plotted on h if it were.

    for i in range(1, 333):

        b = df.iloc[i]

        error = difference(a, b)
        h = b["step"]

        errors.append(error)
        hs.append(h)



    # while i > df.iloc[-1]["step"]:
    #
    #     i = i / 2
    #
    #     print(i)

        # I don't have those values! I could fabricate them...
        # Yeah. Let's do that. And perhaps change our smallest one too.

    # for i in range(9):
    #
    #     j = pow(2, i) - 1
    #
    #     print(j)

    #     a = df.iloc[nums[0]]
    #     b = df.iloc[nums[i]]
    #
    #     error = difference(a, b)
    #     h = (steps[nums[0]] + steps[nums[i]]) / 2
    #
    #     # But... what do I say the h is?
    #
    #     errors.append(error)
    #     hs.append(h)
    #
    # print(errors)
    # print(hs)

    # What the fuck does this mean?
    # I'm MAYBE seeing peaks and troughs.
    # Could very well be random noise, though.
    # And it's so much denser in the other one. Ugh.

    # Okay. Let's make them lists.

    # So asides fro the binary powers thing...



    # def error(row):
    #
    #     return math.sqrt(
    #         (actual["x_2"] - row["x_2"]) ** 2 +
    #         (actual["y_2"] - row["y_2"]) ** 2 +
    #         (actual["z_2"] - row["z_2"]) ** 2
    #     )
    #
    # def h(row):
    #
    #     return abs(actual["step"] - row["step"])


    # df["e"] = df.apply(error, axis=1)
    # df["h"] = df.apply(h, axis=1)
    #
    plt.scatter(hs, errors, alpha=0.5)

    plt.xlim([min(hs), max(hs)])
    plt.ylim([min(errors), max(errors)])

    plt.title("Order of convergence")
    plt.xlabel('h')
    plt.ylabel('error |e|')

    plt.tight_layout()

    plt.savefig("convergence.png")

    plt.show()

    # Time stepping for ODEs (ordinary differential equations). We know the equations.
    # The hard part is just figuring out how many times they're applied.
    # Time slicing: cut the solution into snapshots, with fixed time steps sizes delta_t.
    # We want the numerical solution f_h at fixed time steps.
    # f(0) = f_h(0).
    # Assume that f_h(t) is known. Use Taylor: f_h(t + delta_t) = f(t) + ....
    # Cut after the second term, and insert the ODE. Produces explicit Euler:

    # Ah. Each line is for a different timestep. Each line is composed of different points
    # computed AT that time step.
    # The smaller h, the better our scheme approximates the real solution.
    # The error equals the truncated terms. Makes sense.

    # A scheme is consistent if, as the timestep tends to 0, our approximation tends to the actual value.

    # f(t + delta_t) = f(t) + delta_t . f'(t).
    # f(t + h) = f(t) + h . f'(t).
    # A simple for loop. A single step method: a new solution only needs the previous one.
    # delta_t == h / dt
    # A finite difference method, as we use the difference quotient instead of the derivative.
    # For more complex problems, we don't know the analytical solution, and so cannot directly compute
    # the error.
    # Interesting. But the overshooting is not really relevant to us.
    # An algorithm converge iff it is consistent and stable. It is stable if errors do not blow up.
    # So... perhaps we actually need some lower values. At a shitty approximation, it will definitely blow up.
    # Local error: the error that we make in one step: e(t) = y_h(t) - y(t).
    # At the first time step, this is trivially e(h) = y_h(h) - y(h)

    # Second order because that's where we truncate from.

    # Global error: the error between the numerical solution and true solution (with multiple local steps inbetween).
    # NOT the sum of the local errors. An additional error is propagated, reducing scheme (globally) to first order.

    # An algorithm converges with order p if |y_h(t) - y(t)| <= C dot h^p.

    # An algorithm converges with order p if the absolute difference between the numerical solution
    # and true solution at time t is smaller than some constant multiplied by the size of the timestep
    # raised to p.

    # In other words: an algorithm converges if, for any time t, the absolute difference between the
    # numerical solution and true solution tends to 0 as the timestep size h tends to 0.

    # The error for one h is e_h(t) = (y_h - y)(t)
    # Drop the t, because we always measure things at one point in time.
    # The error for one h is e_h = y_h - y.
    # The error for a more accuration simulation is e_(h/2) = y_(h/2) - y
    # We can combine the two via y: e_h - e_(h/2) = y_h - y_h/2.

    # Okay. I mostly get it. How does this lend itself to MY problem?
    # The error at timestep 2e-8 minus the error at timestep 4e-8 equals y at 2e-8 - y at 4e-8.

    # e is a member of h^p, i.e. e = C dot e^p.

    # We can do this for every single y. Ok.
    # I don't suppose it matters what the difference between the timesteps is, mind you.
    # I wasn't far off!. The trickier part is e.
    # Wait. How do I make the jump?
    # At. Gotcha. Prepare to get dicked, you fuckers.
    # We want to calculate




    # However. We don't HAVE the true solution at any point.

    # Error for one h: e_h() = t(y_h - y)
    # Plug subsequent errors into each other, and we can see how the error has changed...
    # Evaluation of v is stable, i.e. v_m = v(1 + C_e)
    # x_m(t + dt) = x(t) + dt . v(t) + (2 + C) . dt . v(t) . e + x(t) . e
    # Rounding error scales with the solution, and therefore the relative error remains bounded.
    # If we say that it converges to our most accurate result - which seems reasonable - that's a start.
    # We want to know how quickly the error, e_n = x_n - x*, converges to 0.
    # as n tends to infinity (i.e. the size of the representation increases0
    # Hm. This isn't quite what we want, though it does sound interesting.
    # Actually, it is.
    # Let's take our final row as the ground truth of where each point should collide.
    # That will give us a series.
    # And Taylor's theorem somehow fits in, of course.
    # The order estimates the rate of convergence in terms of polynomial behaviour.
    # Isn't this rather easier? If we drew the graph... what would happen?
    # Let's try.


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