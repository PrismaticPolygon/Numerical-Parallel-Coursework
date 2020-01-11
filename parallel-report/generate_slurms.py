import random

# Strong scaling: test with identical parameters but a different number of threads
# Weak scaling: test with different parameters AND different number of threads

template = """#!/bin/sh
#SBATCH --job-name="{}"
#SBATCH -o {}.%A.out
#SBATCH -e {}.%A.err
#SBATCH -p par7.q
#SBATCH -t 02:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=ffgt86@durham.ac.uk
#SBATCH --mail-type=ALL
source /etc/profile.d/modules.sh
module load intel/xe_2017.2
export OMP_NUM THREADS={}
./solution-step5 0.01 5.0 1e-6 $(cat {})"""

# Didn't occur to me that each core could have multiple threads.
# Then again, I'm not 'allowed' to share anyway, soo....
# Maybe I should have some backup runs on MIRA. Yeah, let's do that.
# Just to be safe.

threads = [1, 2, 4, 8, 12, 16, 20, 24]
num_particles = 100

random.seed(0)

def build_particle_string(num_particles, min_mass=1, max_mass=10):

    string = ""

    for i in range(num_particles):

        x_pos = random.random() * 2 - 1
        y_pos = random.random() * 2 - 1
        z_pos = random.random() * 2 - 1

        mass = random.random() * (max_mass - min_mass) + min_mass

        string += " {} {} {} 0 0 0 {}".format(x_pos, y_pos, z_pos, mass)

    return string

for i, thread in enumerate(threads):

    k = num_particles * thread
    conditions = "conditions_{}.txt".format(k)

    with open("conditions_{}.txt".format(k), "w") as conditions_file:

        conditions_file.write(build_particle_string(k))

    strong_name = "strong_{}".format(thread)

    with open(strong_name + ".slurm-script", "w") as strong_file:

        strong_file.write(template.format(strong_name, strong_name, strong_name, thread, "conditions_{}.txt".format(num_particles)))

    weak_name = "weak_{}_{}".format(thread, k)

    with open(weak_name + ".slurm-script", "w") as weak_file:

        weak_file.write(template.format(weak_name, weak_name, weak_name, thread, conditions))





