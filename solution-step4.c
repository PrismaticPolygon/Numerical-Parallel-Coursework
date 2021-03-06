// Translate this file with
//
// icpc -O3 -fopenmp --std=c++11 solution-step4.c -o solution-step4
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>
#include <omp.h>
//#include <chrono>

//using namespace std::chrono;

//auto start = high_resolution_clock::now();


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

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

double** forces;

/**
 * One mass entry per molecule/particle.
 */
double*  mass; // An array of molecule masses [m1, m2, m3]

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


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

  for (int i=0; i<NumberOfBodies; i++) {

    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {

      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);

    }

  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

  if (tPlotDelta<=0.0) {

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

  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  double* forces0 = new double[NumberOfBodies]{0};
  double* forces1 = new double[NumberOfBodies]{0};
  double* forces2 = new double[NumberOfBodies]{0};

  #pragma omp parallel for reduction(std::min:minDx) reduction(+:forces0[:NumberOfBodies]) reduction(+:forces1[:NumberOfBodies]) reduction(+:forces2[:NumberOfBodies])
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
  while (t<=tFinal) {
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

  //auto stop = high_resolution_clock::now();

  //auto duration = duration_cast<microseconds>(stop - start);
  //std::cout << duration.count() << std::endl;

  return 0;
}
