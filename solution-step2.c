// Translate this file with
//
// g++ -O3 --std=c++11 solution-step2.c -o solution-step2
//
// ./solution-step2 0.01 10 1e-8 0 0 0 0 0 0 4 3 0 0 0 0 0 5 3 4 0 0 0 0 3
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
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two particles.
 */
double   minDx;

/*
* Diameter below which particles merge.
*/
double diameter = 0.01;


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
	
  maxV   = 0.0;	// The highest velocity
  minDx  = std::numeric_limits<double>::max();	// The minimum distance between particles
  
  double** forces = new double*[NumberOfBodies];	// A 2D array of the forces on each molecule [[fx1, fy1, fz1], [fx2, fy2, fz2], [fx3, fy3, fz3]]. Here because we can't modify other methods.
  
  for (int i = 0; i < NumberOfBodies; i++) {
	  
	  forces[i] = new double[3]{0.0, 0.0, 0.0};
	  
  }
  
  // Iterate through the particles, from 0 to 1
  for (int i = 0; i < NumberOfBodies - 1; i++) {
	  
	  // Iterate through the particles, from i to 2
	  for (int j = i + 1;  j < NumberOfBodies; j++) {
		  
		  // Calculate the distance from particle i to particle j
		  const double distance = sqrt(
		        (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
		        (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
		        (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
		      );
	  
		  minDx = std::min( minDx,distance );
		  
		  if (distance < diameter) {
			  
			  for (int k = 0; k < 3; k++) {
			  				
				  v[i][k] = (mass[i] / (mass[i] + mass[j])) * v[i][k] + (mass[j] / (mass[i] + mass[j])) * v[j][k];
					
					
			  }
			    			
			  mass[i] += mass[j]; // Merge masses
			
			for (int c = j; c < NumberOfBodies; c++) {	// Remove particle j by left-shifting subsequent particles 
			
				x[c] = x[c+1];			// Co-ordinates
				mass[c] = mass[c+1];	// Mass
				v[c] = v[c+1];			//	Velocity
			
			}	
			
			 NumberOfBodies--;  
			 j--;	// Decrement j as the "old" j has been deleted
			  
			  
		  }	else {
		
			  for (int k = 0; k < 3; k++) {
			  
				  double force = (x[j][k] - x[i][k]) * mass[i] * mass[j] / distance / distance / distance ;
				
				  forces[i][k] += force;
				  forces[j][k] -= force;
				  
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
