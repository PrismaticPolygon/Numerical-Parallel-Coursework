# Numerical-Parallel-Coursework

Steps 1 through 3 can be compiled using:

`g++ -O3 --std=c++11 solution-step2.c -o solution-step2`

and run using:

`./solution-step2 0.01 10 1e-8 0 0 0 0 0 0 4 3 0 0 0 0 0 5 3 4 0 0 0 0 3`

with suitable argument replacements. 

Steps 4 through 6 can be compiled using

`g++ -fopenmp -O3 --std=c++11 solution-step4.c -o solution-step4` or, on Hamilton, `icpc -O3 -fopenmp --std=c++11 solution-step4.c -o solution-step4`.

I would like to note here that despite my best efforts I was unable to get satisfactory results on Hamilton for step 5,
and I could not figure out why. None of my experiments produced the results I obtained experimentally on Mira, and I'm 
fairly sure that my parallelisation is sound. I have briefly documented my method in the report. If you can figure
out why, I'd love to know in my submission mark comments!