// g++ -O3 --std=c++11 solution-step3.c -o solution-step3
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>


double t             = 0;   // Start time.
double tFinal        = 0;   // End time.
double tPlot         = 0;
double tPlotDelta    = 0;
double timeStepSize  = 0.0; // Global time step size used.

int NumberOfBodies = 0;     // Number of particles.

double** x;	                // A 2D array of particle coordinates: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
double** v;                 // A 2D array of particle velocities: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]
double*  mass;              // An array of particle masses [m1, m2, m3]
double** forces;            // A 2D array of particle forces: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]


double   maxV     = 0.0;    // Maximum velocity of all particles.
double   minDx;             // Minimum distance between two particles.
double   diameter = 0.01;   // Diameter below which particles merge.

/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {

  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
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

    if (mass[i] <= 0.0) {

      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);

    }

  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

  if (tPlotDelta <= 0.0) {

    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;

  } else {

    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;

  }

}

std::ofstream videoFile;

/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}

/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "results/result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}

// Tasks can communicate through shared memory.
// Critical sections remain available.
// I should probably check if step3 actually works first.
// Each thread places n/p elements in each local buckets, then thread i gathers the contents
// of bucket i from each processor and sorts.
// We, obviously, don't need to sort.
// It looks like it works, with just a few hitches.
// I'm really very unhappy with my solution.
// But I've just realised that I don't need a 2D array, just a contiguous 1-D array, and I'm good.
// So let's work on that: it should simplify things.
// Aha, that must be the sorting that he means!
// So I have my array of particles. I can map each one to a velocity, and I sort those.
// But I CAN'T lose the index that I start with.

// Parallelise. In particular: parallelise the sorting. We SORT the particles into buckets
// after each time step has passed. It doesn't really matter when, actually.

// Parallel sorting in OpenMP: https://homepages.math.uic.edu/~jan/mcs572/parallelsorting.pdf
// Complex, may be useful: https://www.smaizys.com/programing/bucket-sort-parallel-algorithm-using-c-openmpi/

/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {

  minDx  		  = std::numeric_limits<double>::max();	// The minimum distance between particles
  forces          = new double*[NumberOfBodies];	    // A 2D array of the forces on each molecule
  int numBuckets  = 2;									// The number of buckets

  for (int i = 0; i < NumberOfBodies; i++) {            // Initialise forces on each particle to 0

	forces[i] = new double[3]{0.0, 0.0, 0.0};

  }

  /*
  * 1. Make a 1D array, buckets1D. Each i in bucketsID corresponds to particle i. buckets1D[i] is the bucket particle i is in.
  * 2. Make a 1D array. bucketCounts. The number of particles in each bucket k.
  * 3. Make a 1D array, bucketCounts2. The number of particles currently in each bucket k.
  */

  int* buckets1D = new int[NumberOfBodies];  // Each element i of buckets will map to a particle i. The value will be the bucket that particle i belongs to.
  int* bucketCounts = new int[numBuckets]{0, 0};	 // The number of particles to go in each bucket i
  int* bucketCounts2 = new int[numBuckets]{0, 0};	 // The number of particles in each bucket i

  double vBucket = maxV / (numBuckets - 1);		// The partition

  for (int i = 0; i < NumberOfBodies; i++) {	// For each particle, calculate which bucket it should go into.

    if (maxV == 0) {

      buckets1D[i] = 0;
      bucketCounts[0]++;

    } else {

      double totalV = 0;

      for (int k = 0; k < 3; k++) {	// Loop through each dimension

        totalV += v[i][k] * v[i][k];

      }

      totalV = std::sqrt(totalV);
	  int bucket = round(totalV / vBucket);

      buckets1D[i] = bucket;
      bucketCounts[bucket]++;

  	}

  }

//  std::cout << "Buckets: ";
//
//  for (int i = NumberOfBodies - 1; i >= 0; i--) // Correct.
//
//    std::cout << buckets1D[i];
//
//  std::cout << std::endl;
//
//  std::cout << "Counts: ";

//  for (int i = 0; i < numBuckets; i++) // Correct.
//
//    std::cout << i << ": " << bucketCounts[i] << ", ";
//
//  std::cout << std::endl;

  int** buckets2D = new int*[numBuckets];	// Create a 2D list of buckets using bucketCounts.

  for (int i = 0; i < numBuckets; i++) {

	buckets2D[i] = new int[bucketCounts[i]];

  }

  for (int i = 0; i < NumberOfBodies; i++) {    // Set the current element of bucket k to particle i

    int k = buckets1D[i];

//	std::cout << i << ", " << k << ", " << bucketCounts2[k] << std::endl;

    buckets2D[k][bucketCounts2[k]] = i;
    bucketCounts2[k]++;

  }

  //for (int i = 0; i < numBuckets; i++) { // 000 as expected.

	//std::cout << "Bucket " << i << ": ";

    //for (int j = bucketCounts[i] - 1; j >= 0; j--) {

		//std::cout << buckets2D[i][j];	// Has two three elements in 0, 0 and 2. NOT 1. Why?

    //}

    //std::cout << std::endl;

  //}

  for (int k = 0; k < numBuckets; k++) {	// Iterate through buckets

    int bucketSize = bucketCounts[k];                           // The number of particles in the bucket
    int timeSteps = pow(2, k);									// The number of timesteps to run bucket k for
	double timeStepSizeEuler = timeStepSize / timeSteps;        // The size of the timestep for bucket k

//    std::cout << k << ", " << bucketSize << ", " << timeSteps << ", " << timeStepSize << ", " << timeStepSizeEuler << std::endl;

	for (int q = 0; q < timeSteps; q++) {

      for (int a = 0; a < bucketSize - 1; a++) {

        int i = buckets2D[k][a];

        for (int b = a + 1; b < bucketSize; b++) {

          int j = buckets2D[k][b];

		  std::cout << "(" << i << ", " << j << ")" << std::endl;

          const double distance = sqrt(
            (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
            (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
            (x[i][2] - x[j][2]) * (x[i][2] - x[j][2])
          );

          minDx = std::min( minDx,distance );

          if (distance < diameter) {

            for (int z = 0; z < 3; z++) {

                v[i][z] = (mass[i] / (mass[i] + mass[j])) * v[i][z] + (mass[j] / (mass[i] + mass[j])) * v[j][z];

            }

            mass[i] += mass[j]; // Merge masses

            for (int c = j; c < NumberOfBodies; c++) {	// Remove particle from global arrays

                x[c] = x[c+1];					// Co-ordinates
                mass[c] = mass[c+1];			// Mass
                v[c] = v[c+1];					// Velocity
                buckets1D[c] = buckets1D[c+1];  // Probably unnecessary

            }

            for (int c = b; c < bucketSize; c++) {	    // Remove particle from bucket array

                buckets2D[c] = buckets2D[c+1];

            }

            NumberOfBodies--;
            bucketSize--;
            b--;	// Decrement b as the "old" b has been deleted


          } else {

            for (int z = 0; z < 3; z++) {

                double force = (x[j][z] - x[i][z]) * mass[i] * mass[j] / distance / distance / distance ;

                forces[i][z] += force;
                forces[j][z] -= force;

            }

          }

        }

      }

      for (int a = 0; a < bucketSize - 1; a++) {

        int i = buckets2D[k][a];
        double totalV = 0;

        for (int z = 0; z < 3; z++) {

          x[i][z] = x[i][z] + timeStepSizeEuler * v[i][z];	              // Update coordinates in dimension z
          v[i][z] = v[i][z] + timeStepSizeEuler * forces[i][z] / mass[i]; // Update velocity in dimension z

          totalV += v[i][z] * v[i][z];

        }

        maxV = std::max(maxV, std::sqrt(totalV));

  	  }

    }

  }

  if (NumberOfBodies == 1) {	// Terminate

	tFinal = t;

	std::cout << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

  }

  t += timeStepSize;

  delete[] forces;

}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;

  while (t <= tFinal) {

    updateBody();
    timeStepCounter++;

    if (t >= tPlot) {

      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }

  }

  closeParaviewVideoFile();

  return 0;
}
