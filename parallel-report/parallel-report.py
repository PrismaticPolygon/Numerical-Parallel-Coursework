# Study the scalability of your code and compare it to both a weak scaling and a strong scaling model.

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