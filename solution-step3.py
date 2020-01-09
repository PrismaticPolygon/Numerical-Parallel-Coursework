import subprocess
from random import random

def build_particle_string(num_particles=100, min_mass=1, max_mass=10):

    particles_per_axis = int(num_particles ** (1.0 / 3.0))
    string = ""
    h = 1.0

    for x in range(0, particles_per_axis):

        for y in range(0, particles_per_axis):

            for z in range(0, particles_per_axis):

                x_pos = (random() - 0.5) * 0.9 * h + x * h
                y_pos = (random() - 0.5) * 0.9 * h + y * h
                z_pos = (random() - 0.5) * 0.9 * h + z * h

                mass = random() * (max_mass - min_mass) + min_mass

                string += " {} {} {} 0 0 0 {}".format(x_pos, y_pos, z_pos, mass)

    return string

if __name__ == "__main__":

    subprocess.call(["g++", "-O3", "--std=c++11", "solution-step3.c", "-o", "solution-step3"])

    # args = [
    #    "./solution-step3",                    # Compiled C executable
    #    "0.01",                                     # tPlotDelta
    #    "1.0",                                       # tFinal
    #     "0.1",                                     # timeStepSize
    #     "0", "0", "0",    "0", "0", "0",    "4",    # [x, y, z], [v_x, v_y, v_z], m
    #     "3", "0", "0",    "0", "0", "0",    "5",    # [x, y, z], [v_x, v_y, v_z], m
    #     "3", "4", "0",    "0", "0", "0",    "3"     # [x, y, z], [v_x, v_y, v_z], m
    # ]
    #
    # string = " ".join(args)

    args = [
        "./solution-step3",     # Compiled C executable
       "0.01",                  # tPlotDelta
       "1.01",                   # tFinal
       "1",                  # timeStepSize
    ]

    string = " ".join(args) + build_particle_string(num_particles=8, min_mass=1, max_mass=10)
    
#print(string)    
    
process = subprocess.Popen(string, shell=True)
