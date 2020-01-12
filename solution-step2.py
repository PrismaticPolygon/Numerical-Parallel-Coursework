import subprocess

if __name__ == "__main__":

    subprocess.call(["g++", "-O3", "--std=c++11", "solution-step2.c", "-o", "solution-step2"])

    args = [
        "./solution-step2",                 # Compiled C executable
        "0.01",                             # tPlotDelta
        "3.0",                               # tFinal
        "1e-6",                             # timeStepSize
        "0", "0", "0", "0", "0", "0", "4.0",  # [x, y, z], [v_x, v_y, v_z], m
        "3.0", "0", "0", "0", "0", "0", "5.0",  # [x, y, z], [v_x, v_y, v_z], m
        "3.0", "4.0", "0", "0", "0", "0", "3.0"   # [x, y, z], [v_x, v_y, v_z], m
    ]

    process = subprocess.Popen(args)
