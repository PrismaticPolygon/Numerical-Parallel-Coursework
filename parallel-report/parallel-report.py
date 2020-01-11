# Study the scalability of your code and compare it to both a weak scaling and a strong scaling model.

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import optimize

rc("text", usetex="True")   # Use LaTeX text rendering

def amdahl(N, s, p):

    return 1 / (s + (p / N))

def gustafson(N, s):

    return s + (p * N)

def format_func(value, tick_number):

    float_str = "{0:.2g}".format(value)

    if "e" in float_str:

        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))

    else:

        return float_str

def strong():

    df = pd.read_csv("results.csv")

    t_1 = df.iloc[0]["runtime"]    # Time taken on 1 core

    df["speedup"] = df["runtime"].map(lambda x: t_1 / x)

    df = df[df["type"] == "strong"]

    plt.figure(figsize=(6, 3))  # Default is (6, 4)

    plt.scatter(df["threads"], df["speedup"], label="Data", alpha=0.5)

    params, params_covariance = optimize.curve_fit(amdahl, df["threads"], df["speedup"])

    print("s:", params[0]) # s = 0.03 for Amdahl's law
    print("p:", params[1])

    plt.plot(df["threads"], amdahl(df["threads"], params[0], params[1]), label='Fitted function', color="red")

    plt.xlabel(r'Number of threads $n$')
    plt.ylabel(r'Speedup')

    plt.title(r"Amdahl's law")
    plt.legend(loc='best')

    plt.tight_layout()

    plt.savefig("strong.png")

    plt.show()

def weak():

    df = pd.read_csv("dummy.csv")

    # Hm. I don't think that we can equate speedup in the same way.
    # Instead, we should be plotting
    # Instead, we calculate scaled speedup, and plot that. Right?
    # That's all that we can plot.

    t_1 = df.iloc[0]["runtime"]    # Time taken on 1 core

    df["speedup"] = df["runtime"].map(lambda x: t_1 / x)
    # We put in x, then y.
    # x is the number of threads
    # y is the speedup.

    df = df[df["type"] == "weak"]

    plt.figure(figsize=(6, 3))  # Default is (6, 4)

    plt.scatter(df["threads"], df["speedup"], label="Data", alpha=0.5)

    params, params_covariance = optimize.curve_fit(gustafson, df["threads"], df["speedup"])

    print("s:", params[0]) # s = 0.1 for Gustafson's law
    print("p:", params[1])

    plt.plot(df["threads"], gustafson(df["threads"], params[0], params[1]), label='Fitted function', color="red")

    plt.xlabel(r'Number of threads $n$')
    plt.ylabel(r'Scaled speedup')

    plt.title(r"Gustafson's law")
    plt.legend(loc='best')

    plt.tight_layout()

    plt.savefig("weak.png")

    plt.show()


if __name__ == "__main__":

    strong()

    weak()

# It's... flat. Why? What the fuck does that mean? It's not at all parallelised?


# https://www.kth.se/blogs/pdc/2018/11/scalability-strong-and-weak-scaling/

# Strong scaling: how the solution time varies with the number of processors
# for a fixed total problem size.

# Weak scaling: how the solution time varies with the number of processors
# for a fixed problem size per processor.

# The speedup in parallel computing can be straightforwardly defined as:
# t_1 / t_N, where t_1 is the computational time running on one processor,
# and t_N is the computational time running on N processors.
# Ideally, software should have a linear speedup.

# Amdahl's law: speedup is limited b the fraction of the serial part of
# the software not amenable to parallelisation:

# speedup = 1 / (s + p / N), where s is the proportion of execution time
# spent on the serial part, p is the proportion spent on the parallalised part,
# and N is the number of processors. This means that, for a fixed problem,
# the upper limit of speedup is determined by the serial fraction of the code.
# This is STRONG scaling.

# Gustafson's law is based on the approximations that:
# a) the parallel part p scales linearly with the amount of resources
# b) the serial part does not increas with respect to the size of the problem.

# To test strong scaling: test how the overall computational time of the job
# scales with the number of processing elements.

# To test weak scaling, increase both the job size and the number
# of processing elements.

# For the strong scaling, we'll want to run with identical parameters
# on a different number of threads (1 - 24).

# For the weak scaling, we'll want to increase both the thread count and
# the size of the thing.

# schedule(dynamic) provides better workload distribution at the cost of some extra overhead.

# This gives the formula:
# scaled speedup = s + p x N.
# The scaled speedup increases linearly with respect to the number of processors

# This is going to be difficult to do until I get my supercomputer access.
# Yikes.

# Study the scalability of your code. Compare it to be a weak scaling
# and a strong scaling model.
# Switch OFF file output for your experiments.
# Assess the quality of my parallelisation using tool support (correctness
# and analysis).