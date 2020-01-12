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
  filename << "result-" << counter <<  ".vtp";
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

  double* forces0 = new double[NumberOfBodies]{0};
  double* forces1 = new double[NumberOfBodies]{0};
  double* forces2 = new double[NumberOfBodies]{0};
  int* buckets    = new int[NumberOfBodies]{0};             // The number of buckets

  int numBuckets  = 10;
  double vBucket  = maxV / (numBuckets - 1);		         // The partition

  if (maxV != 0.0) {    // Buckets are already at 0.

    for (int i = 0; i < NumberOfBodies; i++) {

      double totalV = sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

      buckets[i] = round(totalV / vBucket);

    }

  }

  maxV = 0.0;

  for (int i = 0; i < NumberOfBodies - 1; i++) {  // Iterate through particles

    int timeSteps = pow(2, buckets[i]);									// The number of timesteps to run bucket k for
    double timeStepSizeEuler = timeStepSize / timeSteps;                // The size of the timestep for bucket k

    for (int q = 0; q < timeSteps; q++) { // Iterate through timesteps

      for (int j = i + 1; j < NumberOfBodies; j++) {  // Iterate through other bodies

          // std::cout << "Comparing particles " << i << " (bucket " << buckets[i] << ") and " << j << " (bucket " << buckets[j] << ")" << std::endl;

          double distance = (
            (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
            (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
            (x[i][2] - x[j][2]) * (x[i][2] - x[j][2])
          );

          if (distance <= 0.0001) {

            //std::cout << "Merging particles " << i << " (bucket " << buckets[i] << ") and " << j << " (bucket " << buckets[j] << ")" << std::endl;

            double newWeight = mass[i] + mass[j];
		    double weight_i_over = mass[i] / newWeight;
		    double weight_j_over = mass[j] / newWeight;

		    v[i][0] = weight_i_over * v[i][0] + weight_j_over * v[j][0];
		    v[i][1] = weight_i_over * v[i][1] + weight_j_over * v[j][1];
		    v[i][2] = weight_i_over * v[i][2] + weight_j_over * v[j][2];

		    x[i][0] += weight_j_over * (x[j][0] - x[i][0]);
		    x[i][1] += weight_j_over * (x[j][1] - x[i][1]);
		    x[i][2] += weight_j_over * (x[j][2] - x[i][2]);

		    mass[i] = newWeight;

            for (int c = j; c < NumberOfBodies; c++) {	// Remove particle from global arrays

              x[c]       = x[c + 1];					// Co-ordinates
              mass[c]    = mass[c + 1];   			    // Mass
              v[c]       = v[c + 1];					// Velocity
              buckets[c] = buckets[c + 1];              // Buckets

            }

            NumberOfBodies--;
            j--;	// Decrement b as the "old" j has been deleted
            distance = sqrt(distance);

          } else {

            distance = sqrt(distance);

            double weighted_cubed_distance = mass[i] * mass[j] / distance / distance / distance;

            double force0 = (x[j][0] - x[i][0]) * weighted_cubed_distance;
            double force1 = (x[j][1] - x[i][1]) * weighted_cubed_distance;
            double force2 = (x[j][2] - x[i][2]) * weighted_cubed_distance;

            forces0[i] += force0;	// Force caused by j on i in axis 0
            //forces0[j] += -force0;

            forces1[i] += force1;
            //forces1[j] += -force1;

            forces2[i] += force2;
            //forces2[j] += -force2;

          }

          minDx = std::min( minDx,distance );

        }

      x[i][0] = x[i][0] + timeStepSizeEuler * v[i][0];
	  x[i][1] = x[i][1] + timeStepSizeEuler * v[i][1];
	  x[i][2] = x[i][2] + timeStepSizeEuler * v[i][2];

	  v[i][0] = v[i][0] + timeStepSizeEuler * forces0[i] / mass[i];
      v[i][1] = v[i][1] + timeStepSizeEuler * forces1[i] / mass[i];
      v[i][2] = v[i][2] + timeStepSizeEuler * forces2[i] / mass[i];

    }

    double totalV = sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

    maxV = std::max( maxV,totalV );

  }

  if (NumberOfBodies == 1) {	// Terminate

	tFinal = t;

	std::cout << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

  }

  t += timeStepSize;
  
  delete[] forces0;
  delete[] forces1;
  delete[] forces2;

  delete[] buckets;

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
