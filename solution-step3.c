// Translate this file with
//
// g++ -O3 --std=c++11 solution-step3.c -o solution-step3
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

/*
 * All objects should move freely through space. Ensure that the global statistics
 * (minimal distance between all objects and maximum velocity) are still computed correctly
 * Marks are given for correctness and efficiency (try to spot redundant computations).
 * Worth 25 marks. 
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>


double t           = 0;
double tFinal      = 0;
double tPlot       = 0;
double tPlotDelta  = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;	// A 2D array of molecule coordinates: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]

/**
 * Equivalent to x storing the velocities.
 */
double** v; // A 2D array of molecule velocities: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]

/**
 * One mass entry per molecule/particle.
 */
double*  mass; // An array of molecule masses [m1, m2, m3]

/**
 * Equivalent to x storing the forces.
 */
double** forces; // A 2D array of molecule forces: [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV = 0.0;

/**
 * Minimum distance between two particles.
 */
double   minDx;

/*
* Diameter below which particles merge.
*/
double   diameter = 0.01;

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
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
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

/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {

  minDx  		  = std::numeric_limits<double>::max();	// The minimum distance between particles
  forces          = new double*[NumberOfBodies];	    // A 2D array of the forces on each molecule
  int numBuckets  = 10;									// The number of buckets

  for (int i = 0; i < NumberOfBodies; i++) {            // Initialise forces on each particle to 0
	  
	forces[i] = new double[3]{0.0, 0.0, 0.0};
	  
  }

  /*
  * 1. Make a 1D array, buckets1D. Each i in bucketsID corresponds to particle i. buckets1D[i] is the bucket particle i is in.
  * 2. Make a 1D array. bucketCounts. The number of particles in each bucket k.
  * 3. Make a 1D array, bucketCounts2. The number of particles currently in each bucket k.
  */

  int* buckets1D = new int[NumberOfBodies];  // Each element i of buckets will map to a particle i. The value will be the bucket that particle i belongs to.
  int* bucketCounts = new int[numBuckets]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};	 // The number of particles in each bucket i
  int* bucketCounts2 = new int[numBuckets]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0};	 // The number of particles in each bucket i

  double vBucket = maxV / 10;						// The partition

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

      buckets1D[i] = round(totalV / vBucket);
      bucketCounts[i]++;

  	}

  }
  
  int** buckets2D = new int*[numBuckets];	// Create a 2D list of buckets using bucketCounts.
  
  for (int i = 0; i < numBuckets; i++) {
  
	buckets2D = new int[bucketCounts[i]];
  
  }

  for (int i = 0; i < NumberOfBodies; i++) {    // Set the next element of bucket k to i

    int k = buckets1D[i]

    buckets2D[k][bucketCounts2[i]] = i;
    bucketCounts2[i]++;

  }

  for (int k = 0; k < numBuckets; k++) {	// Iterate through buckets

    int bucketSize = bucketCounts[k];                           // The number of particles in the bucket
    int timeSteps = pow(2, k);									// The number of timesteps to run bucket k for
	double timeStepSizeEuler = timeStepSize * (1 / timeSteps);  // The size of the timestep for bucket k

	for (int q = 0; q < timeSteps; q++) {

      for (int a = 0; a < bucketSize - 1; a++) {

        int i = buckets2D[k][a]

        for (int b = a + 1; b < bucketSize; b++) {

          int j == buckets2D[k][b]

          const double distance = sqrt(
            (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
            (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
            (x[i][2] - x[j][2]) * (x[i][2] - x[j][2])
          );

          minDx = std::min(minDx, distance);

          if (distance < diameter) {

            for (int z = 0; z < 3; z++) {

                v[i][z] = (mass[i] / (mass[i] + mass[j])) * v[i][z] + (mass[j] / (mass[i] + mass[j])) * v[j][z];

            }

            mass[i] += mass[j]; // Merge masses

            for (int c = j; c < NumberOfBodies; c++) {	// Remove particle from global arrays

                x[c] = x[c+1];					// Co-ordinates
                mass[c] = mass[c+1];			// Mass
                v[c] = v[c+1];					// Velocity
                buckets1D[c] = buckets1D[c+1]   // Probably unnecessary

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

        int i = buckets2D[k][a]
        double totalV = 0;

        for (int z = 0; z < 3; z++) {

          x[i][z] = x[i][z] + timeStepSizeEuler * v[i][z];	                // Update particle i coordinates in dimension z
          v[i][z] = v[i][z] + timeStepSizeEuler * forces[i][z] / mass[i];	// Update particle i velocity in dimension z

          totalV += v[i][k] * v[i][k];

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
