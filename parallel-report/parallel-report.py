# Study the scalability of your code and compare it to both a weak scaling and a strong scaling model.

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import optimize

rc("text", usetex="True")   # Use LaTeX text rendering

def amdahl(N, s, p):

    return 1 / (s + (p / N))

def gustafson(N, s, p):

    return s + (p * N)

def format_func(value, tick_number):

    float_str = "{0:.2g}".format(value)

    if "e" in float_str:

        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))

    else:

        return float_str

# Need to figure out graphing the weak. Then I'm good for the day. Am I? It's not very late.
# I'll start thinking about my game. Most likely, too my dissertation script has finished running.

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

    plt.xlabel(r'Number of threads $N$')
    plt.ylabel(r'Speedup')

    plt.title(r"Amdahl's law")
    plt.legend(loc='best')

    plt.tight_layout()

    plt.savefig("strong.png")

    plt.show()

global t_1, s_1

def scaled_speedup(row):

    return (t_1 / row["runtime"]) * (row["size"] / s_1)

def weak():

    global t_1, s_1

    df = pd.read_csv("results_2.csv")
    df = df[df["type"] == "weak"]

    t_1 = df.iloc[0]["runtime"]    # Time taken on 1 core
    s_1 = df.iloc[0]["size"]       # Input size of first problem

    df["speedup"] = df.apply(scaled_speedup, axis=1)

    plt.figure(figsize=(6, 3))  # Default is (6, 4)

    plt.scatter(df["threads"], df["speedup"], label="Data", alpha=0.5)

    params, params_covariance = optimize.curve_fit(gustafson, df["threads"], df["speedup"])

    print("s:", params[0])
    print("p:", params[1])

    plt.plot(df["threads"], gustafson(df["threads"], params[0], params[1]), label='Fitted function', color="red")

    plt.xlabel(r'Number of threads $N$')
    plt.ylabel(r'Scaled speedup')

    plt.title(r"Gustafson's law")
    plt.legend(loc='best')

    plt.tight_layout()

    plt.savefig("weak.png")

    plt.show()


if __name__ == "__main__":

    # strong()

    weak()
