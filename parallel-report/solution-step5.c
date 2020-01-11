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
  forces = new double*[NumberOfBodies];

  #pragma omp parallel for
  for (int i = 0; i < NumberOfBodies; i++) {

	 forces[i] = new double[3]{0.0, 0.0, 0.0};

  }

  #pragma omp parallel for reduction(min:minDx)
  for (int k = 0; k < NumberOfBodies * (NumberOfBodies - 1) / 2; k++) {

	size_t i = k / NumberOfBodies, j = k % NumberOfBodies;

	if (j <= i) {

		i = NumberOfBodies - i - 2;
		j = NumberOfBodies - j - 1;

	}

	// Calculate the distance from particle i to particle j
    double distance = sqrt(
      (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
      (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
      (x[i][2] - x[j][2]) * (x[i][2] - x[j][2])
    );

    if (distance != 0) { // Just in case

      double force0 = (x[j][0] - x[i][0]) * mass[i] * mass[j] / distance / distance / distance;
      double force1 = (x[j][1] - x[i][1]) * mass[i] * mass[j] / distance / distance / distance;
      double force2 = (x[j][2] - x[i][2]) * mass[i] * mass[j] / distance / distance / distance;

      forces[i][0] += force0;
      forces[j][0] -= force0;

      forces[i][1] += force1;
      forces[j][1] -= force1;

      forces[i][2] += force2;
      forces[j][2] -= force2;

    }

    if (distance < minDx) {

      minDx = distance;

    }

  }

  #pragma omp parallel for reduction(max:maxV)
  for (int i = 0; i < NumberOfBodies; i++) {

	x[i][0] = x[i][0] + timeStepSize * v[i][0];
	x[i][1] = x[i][1] + timeStepSize * v[i][1];
	x[i][2] = x[i][2] + timeStepSize * v[i][2];

	v[i][0] = v[i][0] + timeStepSize * forces[i][0] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * forces[i][1] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * forces[i][2] / mass[i];

	double totalV = sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

	if (totalV > maxV) {

		maxV = totalV;
    }

  }

  t += timeStepSize;

  for (int i = 0; i < NumberOfBodies; i+) { // Free up memory.

        delete[] forces[i];
  }

  delete[] forces;

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
