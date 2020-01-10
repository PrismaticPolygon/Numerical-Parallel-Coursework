import random

# Strong scaling: test with identical parameters but a different number of threads
# Weak scaling: test with different parameters AND different number of threads

template = """#!/bin/csh
#SBATCH --job-name="{}"
#SBATCH -o {}.%A.out
#SBATCH -e {}.%A.err
#SBATCH -p test.q
#SBATCH -t 00:05:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --cpus-per-task={}
#SBATCH --mail-user=ffgt86@durham.ac.uk
#SBATCH --mail-type=ALL
source /etc/profile.d/modules.sh
module load intel/xe_2017.2
export OMP_NUM THREADS={}
./a.out {}"""

threads = [1, 2, 4, 8, 12, 16, 20, 24]
args = [
       "0.01",                  # tPlotDelta
       "10.0",                  # tFinal
       "1e-7",                  # timeStepSize
]
num_particles = 500

random.seed(0)

def build_particle_string(num_particles, min_mass=1, max_mass=10):

    particles_per_axis = int(num_particles ** (1.0 / 3.0))
    string = ""
    h = 1.0

    for x in range(0, particles_per_axis):

        for y in range(0, particles_per_axis):

            for z in range(0, particles_per_axis):

                x_pos = (random.random() - 0.5) * 0.9 * h + x * h
                y_pos = (random.random() - 0.5) * 0.9 * h + y * h
                z_pos = (random.random() - 0.5) * 0.9 * h + z * h

                mass = random.random() * (max_mass - min_mass) + min_mass

                string += " {} {} {} 0 0 0 {}".format(x_pos, y_pos, z_pos, mass)

    return string

# It seems that we can BATCH a bunch of different jobs. It

for i, thread in enumerate(threads):

    strong_name = "strong_{}".format(thread)

    with open(strong_name + ".slurm-script", "w") as strong_file:

        string = " ".join(args) + build_particle_string(num_particles)

        strong_file.write(template.format(strong_name, strong_name, strong_name, thread, thread, string))

    k = num_particles + (i * 100)

    weak_name = "weak_{}_{}".format(thread, k)

    with open(weak_name + ".slurm-script", "w") as weak_file:

        string = " ".join(args) + build_particle_string(k)

        weak_file.write(template.format(weak_name, weak_name, weak_name, thread, thread, string))





