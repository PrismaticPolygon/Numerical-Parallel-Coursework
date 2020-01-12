// icpc -O3 -fopenmp --std=c++11 solution-step5.c -o solution-step5

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>
#include <omp.h>
#include <chrono>

using namespace std::chrono;

auto start = high_resolution_clock::now();

double t            = 0;
double tFinal       = 0;
double timeStepSize = 0.0; // Global time step size used.

int NumberOfBodies = 0;

double** x;	        // A 2D array of particle coordinates: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
double** v;         // A 2D array of particle velocities: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
double** forces;    // A 2D array of particle forces: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
double*  mass;      // An array of particle masses [m1, m2, m3]

double   maxV;  // Maximum velocity of all particles
double   minDx; // Minimum distance between two particles

void setUp(int argc, char** argv) {

  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i = 0; i < NumberOfBodies; i++) {

    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

  }

  std::cout << "Running... " << std::endl;

}

void updateBody() {

  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  double* forces0 = new double[NumberOfBodies];
  double* forces1 = new double[NumberOfBodies];
  double* forces2 = new double[NumberOfBodies];

  #pragma omp parallel for reduction(std::min:minDx) reduction(+:forces0[:NumberOfBodies]) reduction(+:forces1[:NumberOfBodies]) reduction(+:forces2[:NumberOfBodies])
  for (int k = 0; k < NumberOfBodies * (NumberOfBodies - 1) / 2; k++) {

	size_t i = k / NumberOfBodies, j = k % NumberOfBodies;

	if (j <= i) {

		i = NumberOfBodies - i - 2;
		j = NumberOfBodies - j - 1;

	}

    double distance = sqrt(
      (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
      (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
      (x[i][2] - x[j][2]) * (x[i][2] - x[j][2])
    );

    if (distance != 0) { // Just in case

      double force0 = (x[j][0] - x[i][0]) * mass[i] * mass[j] / distance / distance / distance;
      double force1 = (x[j][1] - x[i][1]) * mass[i] * mass[j] / distance / distance / distance;
      double force2 = (x[j][2] - x[i][2]) * mass[i] * mass[j] / distance / distance / distance;

      forces0[i] += force0;
      forces0[j] += -force0;

      forces1[i] += force1;
      forces1[j] += -force1;

      forces2[i] += force2;
      forces2[j] += -force2;

    }

    minDx = distance;

  }

  #pragma omp parallel for reduction(std::max:maxV)
  for (int i = 0; i < NumberOfBodies; i++) {
	
	x[i][0] = x[i][0] + timeStepSize * v[i][0];
	x[i][1] = x[i][1] + timeStepSize * v[i][1];
	x[i][2] = x[i][2] + timeStepSize * v[i][2];  

	v[i][0] = v[i][0] + timeStepSize * forces0[i] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * forces1[i] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * forces2[i] / mass[i];

	maxV = sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

  }

  t += timeStepSize;

  delete[] forces0;
  delete[] forces1;
  delete[] forces2;

}

/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {

  setUp(argc, argv);

  while (t <= tFinal) {

    updateBody();

  }

  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);

  std::cout << duration.count() << std::endl;

  return 0;

}
