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

// A large RHS (that is, force) induces fast acceleration, so Euler becomes a poor approximation
// For these setups, we should reduce the time-step size h
// Otherwise, we use as large an h as possible.
// So the larger the velocity, the smaller the timestep we should use.
// That's currently implemented.
// New approach. Create our list of buckets, then SORT IT.
// Nope. It's the indices that we care about.
// General idea: slow moving particles shouldn't be close to fast ones.
// Check final lecture for shit on buckets.

// A few, light, fast particles mess up all the time step sizes.
// Sort all particles into buckets according to their speed.
// Assign each bucket characteristic time step (h, h / 2, h / 4...
// And loop over buckets.
// That's it. Let's not forget this sorting teaser, though.

// A stiff problem is a problem that makes (not unconditionally stable) algorithms
// such as explicit schemes using extremely small time step sizes / discretisation lengths.
// Explicit Euler is conditionally stable. We can't say that it IS stable for some h
// unless we fix delta.

// Sort the objects according to their velocity after every timestep.
// Easy.
// It's implied I should be using a sorting algorithm.
// I have an array of velocities, and so I could sort and slice
// those.
// But that is functionally equivalent to what I'm doing at the moment.

// Once we've used maxV to set the buckets, we reset it.

/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {

  minDx  		  = std::numeric_limits<double>::max();	// The minimum distance between particles
  forces          = new double*[NumberOfBodies];	    // A 2D array of the forces on each molecule
  int numBuckets  = 10;
  int* buckets    = new int[NumberOfBodies];            // The number of buckets
  double vBucket  = maxV / (numBuckets - 1);		    // The partition

  for (int i = 0; i < NumberOfBodies; i++) {
	  
	forces[i] = new double[3]{0.0, 0.0, 0.0};           // Initialise forces on each particle to 0

    if (maxV == 0) {

      buckets[i] = 0;

    } else {

      double totalV = 0;

      for (int k = 0; k < 3; k++) {	// Loop through each dimension

        totalV += v[i][k] * v[i][k];

      }

      buckets[i] = round(std::sqrt(totalV) / vBucket);

      //std::cout << "Particle " << i << " in bucket " << buckets[i] << std::endl;

    }
	  
  }

  std::cout << std::endl;

  for (int k = 0; k < numBuckets; k++) {    // Iterate through buckets

    int timeSteps = pow(2, k);									// The number of timesteps to run bucket k for
	double timeStepSizeEuler = timeStepSize / timeSteps;        // The size of the timestep for bucket k

    //std::cout << "Simulating bucket " << k << " (" << timeSteps << " time steps)" << std::endl;


    // Almost. I need to find some way to not compare particles in the same bucket twice.
    // e.g. (0, 2) and (2, 0).
    // I can rely on them being sorted. Right?
    // We don't WANT particles to be recomputed in that way.

    for (int i = 0; i < NumberOfBodies - 1; i++) {  // Iterate through particles

      if (buckets[i] != k) {   // If it is not in bucket k, skip

        continue;

      }

      for (int j = i + 1; j < NumberOfBodies; j++) {

        //std::cout << "Comparing particles " << i << " (bucket " << buckets[i] << ") and " << j << " (bucket " << buckets[j] << ")" << std::endl;

        const double distance = sqrt(
          (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
          (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
          (x[i][2] - x[j][2]) * (x[i][2] - x[j][2])
        );

        minDx = std::min( minDx,distance );

        if (distance < diameter) {

          //std::cout << "Merging particles " << i << " (bucket " << buckets[i] << ") and " << j << " (bucket " << buckets[j] << ")" << std::endl;

          for (int z = 0; z < 3; z++) { // Merge velocities

            v[i][z] = (mass[i] / (mass[i] + mass[j])) * v[i][z] + (mass[j] / (mass[i] + mass[j])) * v[j][z];

          }

          mass[i] += mass[j]; // Merge masses

          for (int c = j; c < NumberOfBodies; c++) {	// Remove particle from global arrays

            x[c] = x[c+1];					// Co-ordinates
            mass[c] = mass[c+1];			// Mass
            v[c] = v[c+1];					// Velocity
            buckets[c] = buckets[c+1];      // Buckets

          }

          NumberOfBodies--;
          j--;	// Decrement b as the "old" b has been deleted

        } else {

          for (int z = 0; z < 3; z++) {

              double force = (x[j][z] - x[i][z]) * mass[i] * mass[j] / distance / distance / distance ;

              forces[i][z] += force;
              forces[j][z] -= force;

          }

        }

      }

    }

  }

  for (int i = 0; i < NumberOfBodies; i++) {	// Update positions and velocity of particles

	double totalV = 0;

	for (int k = 0; k < 3; k++) {

	  x[i][k] = x[i][k] + timeStepSize * v[i][k];	// Update particle i coordinates in dimension k

	  v[i][k] = v[i][k] + timeStepSize * forces[i][k] / mass[i];	// Update particle i velocity in dimension k

	  totalV += v[i][k] * v[i][k];

	}

	maxV = std::max(maxV, std::sqrt(totalV));

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
