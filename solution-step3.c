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
#include <vector>


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

// Not quite what I had in mind.

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

    //std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;

  } else {
  
    //std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
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

// New approach. Create our list of buckets, then SORT IT.
// Nope. It's the indices that we care about.
// General idea: slow moving particles shouldn't be close to fast ones. This is wrong. We're trying to avoid
// truncation errors for particles with high velocities by limiting how fast they move and updating them at smaller increments.
// Perhaps I should write it so that they are ALL in the first bucket.
// That way I don't have to do any of the sorting at the moment. Seems a bit pointless. 
// Check final lecture for shit on buckets.
// Come ON, Dom! Think.
// What is the problem with this solution? It's just crappy in general. It's overcomplex and thus prone to bugs. It's overcomplex
// by Pythonic standards, not by C. A segfault is a very common error, so we should be able to find it. Should? I have little experience
// with C so I don't know how to find segfaults.

/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {

  minDx  		   = std::numeric_limits<double>::max();	// The minimum distance between particles
  forces          = new double*[NumberOfBodies];	    // A 2D array of the forces on each molecule
  int numBuckets  = 3;									// The number of buckets

  for (int i = 0; i < NumberOfBodies; i++) {            // Initialise forces on each particle to 0
	  
	forces[i] = new double[3]{0.0, 0.0, 0.0};
	  
  }
  
  // It shouldn't be the case the these guys can only operate on particles in their own bucket either, hence the sorting. 
  // And it follows that I don't actually need some of my data structures.
  // Right? Where did my brainwave just go?
  // Am I allowed to use other libraries? Is vector a standard libray?

  /*
  * 1. Make a 1D array, buckets1D. Each i in bucketsID corresponds to particle i. buckets1D[i] is the bucket particle i is in.
  * 2. Make a 1D array. bucketCounts. The number of particles in each bucket k.
  * 3. Make a 1D array, bucketCounts2. The number of particles currently in each bucket k.
  */

  int* buckets1D = new int[NumberOfBodies];  // Each element i of buckets will map to a particle i. The value will be the bucket that particle i belongs to.
  int* bucketCounts = new int[numBuckets]{0, 0};	 // The number of particles to go in each bucket i
  int* bucketCounts2 = new int[numBuckets]{0, 0};	 // The number of particles currently in each bucket i

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

	   int bucket = round(std::sqrt(totalV) / vBucket);	

      buckets1D[i] = bucket;	// Particle i goes in bucket
      bucketCounts[bucket]++;	// bucket should have 1 more particle in.

  	}
  	
  	// So that's wrong on the first iteration.
  	// No it isn't. We're saying that every 
  	// Yes it is.

  }
  
  maxV = 0.0;	// Reset maxV for next time around

  std::cout << "buckets1D: ";	// Every particle is in bucket 0, as expected.

  for (int i = NumberOfBodies - 1; i >= 0; i--) { // Correct.

    std::cout << buckets1D[i] << " ";
    
 }

 std::cout << std::endl;
//
  std::cout << "bucketCounts: ";

  for (int i = 0; i < numBuckets; i++) // Correct.

    std::cout << i << ": " << bucketCounts[i] << ", ";

  std::cout << std::endl;
  
  int** buckets2D = new int*[numBuckets];	// Create a 2D list of buckets using bucketCounts.
  
  for (int i = 0; i < numBuckets; i++) {
  
	buckets2D[i] = new int[bucketCounts[i]];	// Set the size of each to the number of particles that should go into it.
  
  }

  for (int i = 0; i < NumberOfBodies; i++) {    // Set the current element of bucket k to particle i

    int k = buckets1D[i];	// Get the bucket of the particle

    int* bucket = buckets2D[k];	// The bucket k that particle i belongs to
    int cur = bucketCounts2[k];  // The current number of particles in bucket k

	std::cout << "Particle " << i << " placed in bucket " << k << " at " << cur << std::endl;
    

    bucket[cur] = i;	// Put particle i in bucket k
    bucketCounts2[k]++;				// Increment the number of particles in bucket k

  }
  
  // This
  // Instant seg fault. Perhaps if my stuff were better... 

  for (int i = 0; i < numBuckets; i++) { // 000 as expected. Yeah, this is wrong.

	std::cout << "Bucket " << i << ": ";

    for (int j = 0; j < bucketCounts[i]; j++) { // This should be the element of each.

		std::cout << buckets2D[i][j] << " ";	// Has two three elements in 0, 0 and 2. NOT 1. Why?

    }

    std::cout << std::endl;

  }

  for (int k = 0; k < numBuckets; k++) {	// Iterate through buckets

    int bucketSize = bucketCounts[k];                           // The number of particles in the bucket
    int timeSteps = pow(2, k);									// The number of timesteps to run bucket k for
	double timeStepSizeEuler = timeStepSize / timeSteps;        // The size of the timestep for bucket k

    std::cout << k << ", " << bucketSize << ", " << timeSteps << ", " << timeStepSize << ", " << timeStepSizeEuler << std::endl;

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
          	
          	std::cout << "Merging" << std::endl;

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

        maxV = std::max( maxV,std::sqrt(totalV) );

  	  }

    }

  }
  
  std::cout << "NumberOfBodies: " << NumberOfBodies << std::endl;  
  
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
      //std::cout << "plot next snapshot"
    	//	    << ",\t time step=" << timeStepCounter
    	std::cout   << "t="         << t;
		//		<< ",\t dt="        << timeStepSize
		//		<< ",\t v_max="     << maxV
		//		<< ",\t dx_min="    << minDx
		//		<< std::endl;

      tPlot += tPlotDelta;
    }
    
  }

  closeParaviewVideoFile();

  return 0;
}
