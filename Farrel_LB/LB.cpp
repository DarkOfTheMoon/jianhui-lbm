#if _WIN32				//Windows headers
#include <Windows.h>
#include <process.h>
#else					//Linux headers
#include <pthread.h>
#include <sys/time.h>
#define _snprintf snprintf
#define Sleep usleep
#endif
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

////Input Variables/////////////////////////////////////////////////////////////////////////////////////////////

char InputFileDefault[1024] = "./Input.txt";			//Default input file if not specified in command line

int LocalThreadNum = 1;								//Number of threads to run

char GeometryFile[1024] = "./Sphere.txt";

bool GeometryFileBin = false;						//Is solids file binary format
bool GeometryFileVTKHeader = false;					//Input file has a VTK header
int NLinesGeometryVTKHeader = 10;					//Number of lines for VTK header

int NLatticeX = 100;								//Lattice dimensions
int NLatticeY = 100;
int NLatticeZ = 100;

double Resolution = 1;								//Voxel physical size in microns

double Viscosity = 0.16667;							//Fluid viscosity

double InitialVx = 0.0;								//Lattice initial velocity
double InitialVy = 0.0;
double InitialVz = 0.0;

int BoundaryConditionX = 0;							//Boundary conditions
int BoundaryConditionY = 0;							//0 = Loop Boundary
int BoundaryConditionZ = 0;							//1 = Solid Boundary

double BodyForceX = 3.0e-6;							//Body force on fluid
double BodyForceY = 0.0;
double BodyForceZ = 0.0;

int TimeStepsMax = 100;								//End simulation after this many timesteps (0 = no limit)
int WriteInterval = 50;								//Write simulation information to console after so many timesteps
int OutputInterval = 1000;							//Output velocity field after so many timesteps

char OutputFilesFolder[1024] = "./";

char OutputVelocityFile[1024] = "./VOut.vtk";

bool OutputVelocityFileBin = false;					//Is velocities file binary format
bool OutputVelocityFileSparse = false;				//Sparse velocity files
bool OutputVelocityFileVTKHeader = true;			//Output file has a VTK header (9 lines)

char OutputGeometryFile[1024] = "./GOut.vtk";

bool OutputGeometryFileVTKHeader = true;
bool OutputGeometryFileBin = false;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define DataType_Bool		0
#define DataType_Short		1
#define DataType_Int		2
#define DataType_Double		3
#define DataType_String		4

struct InputValueStruct{
	char* DataName;
	void* Ptr;
	int PtrLength;		//Max data length
	int DataType;
	bool Required;		//Must be read in
};

InputValueStruct InputFileData[] = {		//Defines entries in the input file to be read in

	//	{ [Property Name], [Variable Ptr], [Data Length], [Data Type (eg DataType_Int)], [must be specified?] },

		{ "NumberOfThreads", &LocalThreadNum, sizeof(LocalThreadNum), DataType_Int, false },

		{ "NLatticeX", &NLatticeX, sizeof(NLatticeX), DataType_Int, true },
		{ "NLatticeY", &NLatticeY, sizeof(NLatticeY), DataType_Int, true },
		{ "NLatticeZ", &NLatticeZ, sizeof(NLatticeZ), DataType_Int, true },

		{ "ResolutionMicron", &Resolution, sizeof(Resolution), DataType_Double, true },
		{ "Viscosity", &Viscosity, sizeof(Viscosity), DataType_Double, true },
		
		{ "BoundaryConditionX", &BoundaryConditionX, sizeof(BoundaryConditionX), DataType_Int, false },
		{ "BoundaryConditionY", &BoundaryConditionY, sizeof(BoundaryConditionY), DataType_Int, false },
		{ "BoundaryConditionZ", &BoundaryConditionZ, sizeof(BoundaryConditionZ), DataType_Int, false },

		{ "InitialVx", &InitialVx, sizeof(InitialVx), DataType_Double, true },
		{ "InitialVy", &InitialVy, sizeof(InitialVy), DataType_Double, true },
		{ "InitialVz", &InitialVz, sizeof(InitialVz), DataType_Double, true },

		{ "BodyForceX", &BodyForceX, sizeof(BodyForceX), DataType_Double, true },
		{ "BodyForceY", &BodyForceY, sizeof(BodyForceY), DataType_Double, true },
		{ "BodyForceZ", &BodyForceZ, sizeof(BodyForceZ), DataType_Double, true },

		{ "TimeStepsMax", &TimeStepsMax, sizeof(TimeStepsMax), DataType_Int, true },
		{ "WriteInterval", &WriteInterval, sizeof(WriteInterval), DataType_Int, true },
		
		{ "InputGeometryFile", GeometryFile, sizeof(GeometryFile), DataType_String, true },
			{ "GeometryFileBin", &GeometryFileBin, sizeof(GeometryFileBin), DataType_Bool, true },
			{ "GeometryFileVTKHeader", &GeometryFileVTKHeader, sizeof(GeometryFileVTKHeader), DataType_Bool, true },
			{ "NLinesGeometryVTKHeader", &NLinesGeometryVTKHeader, sizeof(NLinesGeometryVTKHeader), DataType_Int, true },

		{ "OutputGeometryFile", OutputGeometryFile, sizeof(OutputGeometryFile), DataType_String, true },
			{ "OutputGeometryFileBin", &OutputGeometryFileBin, sizeof(OutputGeometryFileBin), DataType_Bool, true },
			{ "OutputGeometryFileVTKHeader", &OutputGeometryFileVTKHeader, sizeof(OutputGeometryFileVTKHeader), DataType_Bool, true },
		{ "OutputVelocityFile", OutputVelocityFile, sizeof(OutputVelocityFile), DataType_String, true },
			{ "OutputVelocityFileBin", &OutputVelocityFileBin, sizeof(OutputVelocityFileBin), DataType_Bool, true },
			{ "OutputVelocityFileSparse", &OutputVelocityFileSparse, sizeof(OutputVelocityFileSparse), DataType_Bool, true },
			{ "OutputVelocityFileVTKHeader", &OutputVelocityFileVTKHeader, sizeof(OutputVelocityFileVTKHeader), DataType_Bool, true },
		

	};

class SimulationTimer{			//Timer class to record elapsed time

private:

#if _WIN32	//Windows
	LARGE_INTEGER frequency;        
	LARGE_INTEGER SimulationStart;
	LARGE_INTEGER LastStepTime;
#else
	timeval SimulationStart;
	timeval LastStepTime;
#endif

public:

	SimulationTimer(){

#if _WIN32	//Windows
		QueryPerformanceFrequency(&frequency);
		QueryPerformanceCounter(&SimulationStart);

		LastStepTime = SimulationStart;
#else
		gettimeofday(&SimulationStart, 0);

		LastStepTime.tv_sec = SimulationStart.tv_sec;
		LastStepTime.tv_usec = SimulationStart.tv_usec;
#endif
		
	}

	double GetTimeSinceLastStep(){		//Time since this function last called in seconds

#if _WIN32		//Windows
		LARGE_INTEGER StepTime;
		QueryPerformanceCounter(&StepTime);

		double dt = (double)(StepTime.QuadPart - LastStepTime.QuadPart) / (double)frequency.QuadPart;

		LastStepTime = StepTime;
#else
		timeval StepTime;
		gettimeofday(&StepTime, 0);

		double dt = (double)(StepTime.tv_sec - LastStepTime.tv_sec) + (double)(StepTime.tv_usec - LastStepTime.tv_usec)/1000000.0;

		LastStepTime.tv_sec = StepTime.tv_sec;
		LastStepTime.tv_usec = StepTime.tv_usec;
#endif

		return dt;		//dt in seconds
	}

	double GetTimeSinceLastStepms(){	//Time since this function last called in milliseconds

		return (GetTimeSinceLastStep() * 1000);
	}

	double GetSimulationTime(){			//Time since this SimulationTimer created

#if _WIN32		//Windows
		LARGE_INTEGER Time;
		QueryPerformanceCounter(&Time);

		double dt = (double)(Time.QuadPart - SimulationStart.QuadPart) / (double)frequency.QuadPart;
#else
		timeval Time;
		gettimeofday(&Time, 0);

		double dt = (double)(Time.tv_sec - SimulationStart.tv_sec) + (double)(Time.tv_usec - SimulationStart.tv_usec)/1000000.0;
#endif

		return dt;		//dt in seconds

	}

	void Reset(){						//Set timer start time to now

#if _WIN32	//Windows
		QueryPerformanceFrequency(&frequency);
		QueryPerformanceCounter(&SimulationStart);

		LastStepTime = SimulationStart;
#else
		gettimeofday(&SimulationStart, 0);

		LastStepTime.tv_sec = SimulationStart.tv_sec;
		LastStepTime.tv_usec = SimulationStart.tv_usec;
#endif

	}
	
};

struct ThreadStruct{
	int ThreadID;		//Thread ID
	int Status;			//Synchronisation variable
	int HoldStatus;		//Synchronisation variable
	long long Index0;	//Start index in node array
	long long Index1;	//End index in node array
	double USqTotal;	//Average velocity variable for the nodes under this thread
};

ThreadStruct* Threads;	//Array of thread data

////////////////////


//Calculation variables


long long nVoxels;		//Total number of voxels (= NLatticeX * NLatticeY * NLatticeZ)
long long nNonSolid;		//Number of non-solid nodes

long long* Solid;	//Array, nX*nY*nZ, -1 = Solid, 0 -> NVoxels-1 = index for voxel in data arrays

inline void SetSolid(int x, int y, int z, long long Value){
	Solid[NLatticeX*(NLatticeY*z + y) + x] = Value;
}

inline long long GetSolid(int x, int y, int z){
	return Solid[NLatticeX*(NLatticeY*z + y) + x];
}

struct Node{	
	double f[19];					//Distribution function over basis vectors
	double fCalc[19];				//Temporary storage before streaming step
	double Velocity[3];				//Node velocity
	double Density;					//Node density
	Node* Neighbour[19];			//Ptr to neighbouring node for each vector
};

Node* Nodes;						//Node data. Node[0] = Solid

const int BasisVector[19][3] = {	{ 0, 0, 0},		//Centre node
									{ 1, 0, 0},		//+X
									{-1, 0, 0},		//-X
									{ 0, 1, 0},		//   +Y
									{ 0,-1, 0},		//   -Y
									{ 0, 0, 1},		//      +Z
									{ 0, 0,-1},		//      -Z
									{ 1, 1, 0},		//+X +Y
									{-1, 1, 0},		//-X +Y
									{ 1,-1, 0},		//+X -Y
									{-1,-1, 0},		//-X -Y
									{ 0, 1, 1},		//   +Y +Z
									{ 0,-1, 1},		//   -Y +Z
									{ 0, 1,-1},		//   +Y -Z
									{ 0,-1,-1},		//   -Y -Z
									{ 1, 0, 1},		//+X    +Z
									{-1, 0, 1},		//-X    +Z
									{ 1, 0,-1},		//+X    -Z
									{-1, 0,-1}		//-X    -Z
								};

const double Weighting[19] = {	1.0/3.0,			//Centre node
								1.0/18.0,			//+X
								1.0/18.0,			//-X
								1.0/18.0,			//   +Y
								1.0/18.0,			//   -Y
								1.0/18.0,			//      +Z
								1.0/18.0,			//      -Z
								1.0/36.0,			//+X +Y
								1.0/36.0,			//-X +Y
								1.0/36.0,			//+X -Y
								1.0/36.0,			//-X -Y
								1.0/36.0,			//   +Y +Z
								1.0/36.0,			//   -Y +Z
								1.0/36.0,			//   +Y -Z
								1.0/36.0,			//   -Y -Z
								1.0/36.0,			//+X    +Z
								1.0/36.0,			//-X    +Z
								1.0/36.0,			//+X    -Z
								1.0/36.0			//-X    -Z
							};

const double WeightedBasis[19][3] = {	{ 0, 0, 0},						//Centre node
										{ 1.0/18.0, 0, 0},				//+X
										{-1.0/18.0, 0, 0},				//-X
										{ 0, 1.0/18.0, 0},				//   +Y
										{ 0,-1.0/18.0, 0},				//   -Y
										{ 0, 0, 1.0/18.0},				//      +Z
										{ 0, 0,-1.0/18.0},				//      -Z
										{ 1.0/36.0, 1.0/36.0, 0},		//+X +Y
										{-1.0/36.0, 1.0/36.0, 0},		//-X +Y
										{ 1.0/36.0,-1.0/36.0, 0},		//+X -Y
										{-1.0/36.0,-1.0/36.0, 0},		//-X -Y
										{ 0, 1.0/36.0, 1.0/36.0},		//   +Y +Z
										{ 0,-1.0/36.0, 1.0/36.0},		//   -Y +Z
										{ 0, 1.0/36.0,-1.0/36.0},		//   +Y -Z
										{ 0,-1.0/36.0,-1.0/36.0},		//   -Y -Z
										{ 1.0/36.0, 0, 1.0/36.0},		//+X    +Z
										{-1.0/36.0, 0, 1.0/36.0},		//-X    +Z
										{ 1.0/36.0, 0,-1.0/36.0},		//+X    -Z
										{-1.0/36.0, 0,-1.0/36.0}		//-X    -Z
									};

								// 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
const double Matrix[19][19] = {	{  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},	//0		//Collision matrix
								{-30,-11,-11,-11,-11,-11,-11,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8},	//1
								{ 12, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},	//2
								{  0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1},	//3
								{  0, -4,  4,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1},	//4
								{  0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  1, -1,  1, -1,  0,  0,  0,  0},	//5
								{  0,  0,  0, -4,  4,  0,  0,  1,  1, -1, -1,  1, -1,  1, -1,  0,  0,  0,  0},	//6
								{  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1}, 	//7
								{  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1},	//8
								{  0,  2,  2, -1, -1, -1, -1,  1,  1,  1,  1, -2, -2, -2, -2,  1,  1,  1,  1},	//9
								{  0, -4, -4,  2,  2,  2,  2,  1,  1,  1,  1, -2, -2, -2, -2,  1,  1,  1,  1},	//10
								{  0,  0,  0,  1,  1, -1, -1,  1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1},	//11
								{  0,  0,  0, -2, -2,  2,  2,  1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1},	//12
								{  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0},	//13
								{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0},	//14
								{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  -1,-1,  1},	//15
								{  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0, -1,  1, -1,  1},	//16
								{  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  1, -1,  1, -1,  0,  0,  0,  0},	//17
								{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  1,  1, -1, -1}	//18
							};

const double MatrixInv[19][19] = {		//Inverse of collision matrix
		//			0							1						2							3						4							5						6							7						8							9						10							11							12						13						14							15							16						17							18
		{  5.2631578947368418e-002, -1.2531328320802001e-002,  4.7619047619047623e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000 },
		{  5.2631578947368460e-002, -4.5948203842940691e-003, -1.5873015873015876e-002,  1.0000000000000001e-001, -1.0000000000000001e-001,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  5.5555555555555552e-002, -5.5555555555555546e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000 },
		{  5.2631578947368460e-002, -4.5948203842940691e-003, -1.5873015873015872e-002, -1.0000000000000001e-001,  1.0000000000000001e-001,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  5.5555555555555546e-002, -5.5555555555555559e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000 },
		{  5.2631578947368432e-002, -4.5948203842940691e-003, -1.5873015873015872e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  9.9999999999999992e-002, -1.0000000000000001e-001,  0.0000000000000000e+000,  0.0000000000000000e+000, -2.7777777777777776e-002,  2.7777777777777776e-002,  8.3333333333333329e-002, -8.3333333333333343e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000 },
		{  5.2631578947368432e-002, -4.5948203842940691e-003, -1.5873015873015872e-002,  0.0000000000000000e+000,  0.0000000000000000e+000, -9.9999999999999992e-002,  1.0000000000000001e-001,  0.0000000000000000e+000,  0.0000000000000000e+000, -2.7777777777777766e-002,  2.7777777777777783e-002,  8.3333333333333329e-002, -8.3333333333333343e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000 },
		{  5.2631578947368425e-002, -4.5948203842940691e-003, -1.5873015873015872e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.0000000000000001e-001, -1.0000000000000001e-001, -2.7777777777777773e-002,  2.7777777777777780e-002, -8.3333333333333329e-002,  8.3333333333333343e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000 },
		{  5.2631578947368432e-002, -4.5948203842940691e-003, -1.5873015873015872e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000, -9.9999999999999992e-002,  1.0000000000000001e-001, -2.7777777777777759e-002,  2.7777777777777783e-002, -8.3333333333333329e-002,  8.3333333333333343e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000 },
		{  5.2631578947368432e-002,  3.3416875522138691e-003,  3.9682539682539698e-003,  1.0000000000000001e-001,  2.5000000000000001e-002,  1.0000000000000001e-001,  2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.7777777777777780e-002,  1.3888888888888890e-002,  8.3333333333333329e-002,  4.1666666666666664e-002,  2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.2500000000000000e-001, -1.2500000000000000e-001,  0.0000000000000000e+000 },
		{  5.2631578947368432e-002,  3.3416875522138700e-003,  3.9682539682539698e-003, -1.0000000000000001e-001, -2.5000000000000001e-002,  1.0000000000000001e-001,  2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.7777777777777780e-002,  1.3888888888888890e-002,  8.3333333333333329e-002,  4.1666666666666664e-002, -2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.2500000000000000e-001, -1.2500000000000000e-001,  0.0000000000000000e+000 },
		{  5.2631578947368432e-002,  3.3416875522138691e-003,  3.9682539682539698e-003,  1.0000000000000001e-001,  2.5000000000000001e-002, -1.0000000000000001e-001, -2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.7777777777777780e-002,  1.3888888888888890e-002,  8.3333333333333329e-002,  4.1666666666666664e-002, -2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.2500000000000000e-001,  1.2500000000000000e-001,  0.0000000000000000e+000 },
		{  5.2631578947368432e-002,  3.3416875522138700e-003,  3.9682539682539698e-003, -1.0000000000000001e-001, -2.5000000000000001e-002, -1.0000000000000001e-001, -2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.7777777777777780e-002,  1.3888888888888890e-002,  8.3333333333333329e-002,  4.1666666666666664e-002,  2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.2500000000000000e-001,  1.2500000000000000e-001,  0.0000000000000000e+000 },
		{  5.2631578947368411e-002,  3.3416875522138678e-003,  3.9682539682539680e-003,  0.0000000000000000e+000,  0.0000000000000000e+000,  9.9999999999999978e-002,  2.4999999999999994e-002,  9.9999999999999978e-002,  2.4999999999999994e-002, -5.5555555555555552e-002, -2.7777777777777776e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.2500000000000003e-001, -1.2500000000000003e-001 },
		{  5.2631578947368446e-002,  3.3416875522138687e-003,  3.9682539682539698e-003,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.0000000000000001e-001, -2.5000000000000001e-002,  1.0000000000000001e-001,  2.5000000000000001e-002, -5.5555555555555552e-002, -2.7777777777777776e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000, -2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.2500000000000000e-001, -1.2500000000000000e-001 },
		{  5.2631578947368418e-002,  3.3416875522138687e-003,  3.9682539682539698e-003,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.0000000000000001e-001,  2.5000000000000001e-002, -1.0000000000000001e-001, -2.5000000000000001e-002, -5.5555555555555559e-002, -2.7777777777777780e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000, -2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.2500000000000000e-001,  1.2500000000000000e-001 },
		{  5.2631578947368446e-002,  3.3416875522138687e-003,  3.9682539682539698e-003,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.0000000000000001e-001, -2.5000000000000001e-002, -1.0000000000000001e-001, -2.5000000000000001e-002, -5.5555555555555546e-002, -2.7777777777777773e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.5000000000000000e-001,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.2500000000000000e-001,  1.2500000000000000e-001 },
		{  5.2631578947368432e-002,  3.3416875522138683e-003,  3.9682539682539698e-003,  1.0000000000000001e-001,  2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.0000000000000001e-001,  2.5000000000000001e-002,  2.7777777777777780e-002,  1.3888888888888890e-002, -8.3333333333333329e-002, -4.1666666666666664e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.5000000000000000e-001, -1.2500000000000000e-001,  0.0000000000000000e+000,  1.2500000000000000e-001 },
		{  5.2631578947368432e-002,  3.3416875522138691e-003,  3.9682539682539698e-003, -1.0000000000000001e-001, -2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  1.0000000000000001e-001,  2.5000000000000001e-002,  2.7777777777777780e-002,  1.3888888888888890e-002, -8.3333333333333329e-002, -4.1666666666666664e-002,  0.0000000000000000e+000,  0.0000000000000000e+000, -2.5000000000000000e-001,  1.2500000000000000e-001,  0.0000000000000000e+000,  1.2500000000000000e-001 },
		{  5.2631578947368432e-002,  3.3416875522138683e-003,  3.9682539682539698e-003,  1.0000000000000001e-001,  2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.0000000000000001e-001, -2.5000000000000001e-002,  2.7777777777777780e-002,  1.3888888888888890e-002, -8.3333333333333329e-002, -4.1666666666666664e-002,  0.0000000000000000e+000,  0.0000000000000000e+000, -2.5000000000000000e-001, -1.2500000000000000e-001,  0.0000000000000000e+000, -1.2500000000000000e-001 },
		{  5.2631578947368432e-002,  3.3416875522138691e-003,  3.9682539682539698e-003, -1.0000000000000001e-001, -2.5000000000000001e-002,  0.0000000000000000e+000,  0.0000000000000000e+000, -1.0000000000000001e-001, -2.5000000000000001e-002,  2.7777777777777780e-002,  1.3888888888888890e-002, -8.3333333333333329e-002, -4.1666666666666664e-002,  0.0000000000000000e+000,  0.0000000000000000e+000,  2.5000000000000000e-001,  1.2500000000000000e-001,  0.0000000000000000e+000, -1.2500000000000000e-001 }

	};

double SDiag[19];		//Diagonal of the relaxation-time matrix S, calculated in Initialise
double dtI_SDiag[19];	//dt * (I - 0.5S) Matrix (diagonal), calculated in Initialise

//Quantities related to the lattice constant 'c' in the equilibrium distribution function, precalculated for speed in Initialise()
double cEq;				//c = sqrt(3RT); R=8.31, T=Temp OR dx/dt 'microscopic lattice velocity'
double cEqSq;			//c^2
double Inv_cEqSq;		//1/(c^2)
double cEq_3;			//c/3
//

const int ReflMap[19] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15};	//For bounce-back condition, maps streaming vector onto its opposite

//Function declarations
	
	//Input File reading
int ReadInputFile(char* InputFileName);
bool ReadBool(char* Value);
bool CompareStr(char* DataName, char* Property);
int ReadProperty(ifstream& File, char* Property, char* Value);

	//Program
bool Initialise();											//Initialise arrays and read in from files
void Finalise();											//Delete arrays
bool ReadGeometry(long long* NumberOfNonSolids);			//Read in geometry file
bool ReadGeometryBin(long long* NumberOfNonSolids);		//Read in geometry file binary format
void OutputVelocities(int nTimeSteps);						//Output velocity field file
void OutputGeometry();										//Output geometry file
class SimulationTimer;										//Timer class
struct ThreadStruct;										//Thread data variable
void SyncThreads(int ThreadID, int ThreadNum, bool Fast);
#if _WIN32													
unsigned int __stdcall ThreadMain(void* data);				//Windows multi-threading
#else														
void* ThreadMain(void* data);								//Linux multi-threading
#endif

	//Calculation
void Collision(long long Index0, long long Index1);
void ComputeMacroscopicVariables(long long Index0, long long Index1, double* UAverage);
double fEquilibrium  (double Density, int Basis, double Velocity[3]);		//Equilibrium function for any basis vector
double fEquilibrium0 (double Density, double uSq, double Velocity[3]);		//Equilibrium function for each basis vector
double fEquilibrium1 (double Density, double uSq, double Velocity[3]);
double fEquilibrium2 (double Density, double uSq, double Velocity[3]);
double fEquilibrium3 (double Density, double uSq, double Velocity[3]);
double fEquilibrium4 (double Density, double uSq, double Velocity[3]);
double fEquilibrium5 (double Density, double uSq, double Velocity[3]);
double fEquilibrium6 (double Density, double uSq, double Velocity[3]);
double fEquilibrium7 (double Density, double uSq, double Velocity[3]);
double fEquilibrium8 (double Density, double uSq, double Velocity[3]);
double fEquilibrium9 (double Density, double uSq, double Velocity[3]);
double fEquilibrium10(double Density, double uSq, double Velocity[3]);
double fEquilibrium11(double Density, double uSq, double Velocity[3]);
double fEquilibrium12(double Density, double uSq, double Velocity[3]);
double fEquilibrium13(double Density, double uSq, double Velocity[3]);
double fEquilibrium14(double Density, double uSq, double Velocity[3]);
double fEquilibrium15(double Density, double uSq, double Velocity[3]);
double fEquilibrium16(double Density, double uSq, double Velocity[3]);
double fEquilibrium17(double Density, double uSq, double Velocity[3]);
double fEquilibrium18(double Density, double uSq, double Velocity[3]);

//

int main(int argc, char *argv[]){

	//Input file

	char* InputFile = InputFileDefault;

	if(argc >= 2){		//Assume file specified as second command argument
		InputFile = argv[1];
	}else{
		cout << "No input file specified: defaulting to \"" << InputFileDefault << "\"" << endl;
	}

	int r = ReadInputFile(InputFile);

	if(r == 1){
		cout << "The input file \"" << InputFile << "\" could not be opened" << endl;
		return 0;
	}else if(r == 2){
		cout << "The file did not specify all required input values" << endl;
		return 0;
	}

	int ThreadVxN;


	//Initialisation

	bool ret = Initialise();		//Initialise variables and arrays

	if(!ret){						//Memory allocation error
		goto End;
	}

	OutputGeometry();				//Output geometry file


	//Prepare threads

	Threads = new ThreadStruct[LocalThreadNum];

	ThreadVxN = (int)floor((double)nNonSolid / (double)LocalThreadNum);

	for(int i=0; i!=LocalThreadNum; i++){
		Threads[i].ThreadID = i;
		Threads[i].Status = 0;
		Threads[i].HoldStatus = 0;
		Threads[i].Index0 = i * ThreadVxN + min( nNonSolid%ThreadVxN, (long long)i) + 1;							//Thread's first index in node array (starts from 1)
		Threads[i].Index1 = Threads[i].Index0 + ThreadVxN + ( (i < nNonSolid%ThreadVxN) ? 1 : 0 ) - 1;	//Index of final node in thread's array (ends at, including, nNonSolid)
		
		if(i!=0){
#if _WIN32															
			_beginthreadex(0, 0, ThreadMain, &Threads[i], 0, 0);		//Start thread (Windows)
#else
			pthread_t tid;
			pthread_create(&tid, 0, &ThreadMain, (void*)(&Threads[i]));	//Start thread (Linux)
#endif
		}

	}

	ThreadMain(&Threads[0]);	//Simulation carried out in ThreadMain

	//
	
	delete[] Threads;

End:
	Finalise();				//Delete variables
	return 0;
}

		#if _WIN32													
unsigned int __stdcall ThreadMain(void* data){				//Windows multi-threading
		#else														
void* ThreadMain(void* data){								//Linux multi-threading
		#endif

	int ThreadID = ((ThreadStruct*)data)->ThreadID;

	SyncThreads(ThreadID, LocalThreadNum, true);

	//Simulation

	SimulationTimer* Timer = new SimulationTimer();		//Calculates calculation time

	//cout << "Starting Calculation" << endl;

	double UAverage = 0;		//Average velocity
	double UAverageLast = 0;	//Previous average velocity

	long long i0 = ((ThreadStruct*)data)->Index0;
	long long i1 = ((ThreadStruct*)data)->Index1;

	for(int n=1; n<=TimeStepsMax; n++){

		bool WriteSync = false;


		Collision(i0, i1);

			SyncThreads(ThreadID, LocalThreadNum, true);	//Ensure all threads have completed their calculation before continuing


		ComputeMacroscopicVariables(i0, i1, &(((ThreadStruct*)data)->USqTotal));

			SyncThreads(ThreadID, LocalThreadNum, true);


		if(n % WriteInterval == 0){								//Write to console

			WriteSync = true;

		if(ThreadID==0){

			UAverage = 0;
			for(int i=0; i!=LocalThreadNum; i++){
				UAverage += Threads[i].USqTotal;
			}
			UAverage = UAverage/nNonSolid;

			cout << "After " << n << " time-steps:" << endl;

			double UChange = 100 * (UAverage - UAverageLast)/UAverageLast;	//Percentage change in average velocity since last write interval

			if(n == WriteInterval){		//First output

				cout << "Average Velocity Magnitude = " << UAverage << endl;

			}else{

				cout << "Average Velocity Magnitude = " << UAverage << "(";
			
				if(UChange>0) 
					cout << '+';
			
				cout << UChange << "%)" << endl;

			}

			UAverageLast = UAverage;

			cout << "Elapsed time: " << Timer->GetSimulationTime() << endl;

			cout << endl;
		}
		}

		if(OutputInterval!=0 && n % OutputInterval == 0 && n!=TimeStepsMax){		//Output velocity field file

			WriteSync = true;

		if(ThreadID==0){

			OutputVelocities(n);

		}
		}

		if(WriteSync){										//Wait for thread 0 to finish writing to console/output file
			SyncThreads(ThreadID, LocalThreadNum, false);
			WriteSync = false;
		}


	}

	if(ThreadID!=0){
		return 0;
	}

	
	//Final outputs

	cout << "Simulation end after " << Timer->GetSimulationTime() << " seconds" << endl;

	OutputVelocities(TimeStepsMax);

	delete Timer;
	

	//

	return 0;
}

void SyncThreads(int ThreadID, int ThreadNum, bool Fast){

	Threads[ThreadID].Status++;			//No longer hold status

	while(true){
		bool Proceed = true;
		for(int i=0; i!=ThreadNum; i++){
			if(Threads[i].Status == Threads[ThreadID].HoldStatus){ Proceed = false; }
		}
		if(Proceed){ break; }
		if(!Fast){ Sleep(1); }			//Immediate continuation not critical
	}

	Threads[ThreadID].HoldStatus++;

}

void ComputeMacroscopicVariables(long long Index0, long long Index1, double* USqTotal){

	double VAv[3] = {0, 0, 0};		//Average velocity
	double UAv    = 0;				//Magnitude of average velocity

	const double dt = 1.0;

	double dtF_2[3];
	
	dtF_2[0] = 0.5*dt*BodyForceX;
	dtF_2[1] = 0.5*dt*BodyForceY;
	dtF_2[2] = 0.5*dt*BodyForceZ;

	for(long long i=Index0; i<=Index1; i++){

		double* f = Nodes[i].f;
		double* F = Nodes[i].fCalc;
		double* u = Nodes[i].Velocity;

		double density = 0;

		for(int k=0;k<19;k++){
			f[k] = F[k];
			density += f[k];
		}

		double u0 = (1.0/18.0) * (f[1] - f[2] + 0.5*(f[7]  - f[8]  + f[9]  - f[10] + f[15] - f[16] + f[17] - f[18]));
		double u1 = (1.0/18.0) * (f[3] - f[4] + 0.5*(f[7]  + f[8]  - f[9]  - f[10] + f[11] - f[12] + f[13] - f[14]));
		double u2 = (1.0/18.0) * (f[5] - f[6] + 0.5*(f[11] + f[12] - f[13] - f[14] + f[15] + f[16] - f[17] - f[18]));

		double InvDensity = 1/density;

		u[0] = (u0 + dtF_2[0])*InvDensity;
		u[1] = (u1 + dtF_2[1])*InvDensity;
		u[2] = (u2 + dtF_2[2])*InvDensity;

		Nodes[i].Density = density;

		VAv[0] += u[0];
		VAv[1] += u[1];
		VAv[2] += u[2];

	}

	UAv = sqrt(VAv[0]*VAv[0] + VAv[1]*VAv[1] + VAv[2]*VAv[2]);	//Add up USq for these nodes to average velocity calculation

	*USqTotal = UAv;

}

void Collision(long long Index0, long long Index1){		//Do collision step on nodes Index0 to Index1 inclusive

	for(long long i=Index0; i<=Index1; i++){

		double* f = Nodes[i].f;				//Node distribution
		double* F = Nodes[i].fCalc;			//Temporary node distribution
		double* u = Nodes[i].Velocity;		//Node velocity

		double Density = Nodes[i].Density;	//Node density
		double USq = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];		//Magnitude of node velocity


		double GuoF[19];					//Body force vector

		//Calculate body force vector GuoF[19]

			//Pre calc frequently-used constants (16 multiplications)

		double FU[3][3];					//Apologise for expletive name: tensor product of body force vector and node velocity vector
		
		FU[0][0] = BodyForceX * u[0];	FU[0][1] = BodyForceX * u[1];	FU[0][2] = BodyForceX * u[2];
		FU[1][0] = BodyForceY * u[0];	FU[1][1] = BodyForceY * u[1];	FU[1][2] = BodyForceY * u[2];
		FU[2][0] = BodyForceZ * u[0];	FU[2][1] = BodyForceZ * u[1];	FU[2][2] = BodyForceZ * u[2];

		double FC_3[3];						//Body force vector multiplied by c/3

		FC_3[0] = BodyForceX * cEq_3;
		FC_3[1] = BodyForceY * cEq_3;
		FC_3[2] = BodyForceZ * cEq_3;

		double FdotU_3 = (1.0/3.0) * (FU[0][0] + FU[1][1] + FU[2][2]);			//Force dot Velocity / 3 

		double WConst[3];					//The three different constants, dependant on basis vector weighting
		
		WConst[0] = 9.0 * Weighting[0] * Inv_cEqSq;		// = 9W/c^2 for vector  0
		WConst[1] = 9.0 * Weighting[1] * Inv_cEqSq;		// = 9W/c^2 for vectors 1 -> 6
		WConst[2] = 9.0 * Weighting[7] * Inv_cEqSq;		// = 9W/c^2 for vectors 7 -> 18

		//

		GuoF[0]  = WConst[0] * ( - FdotU_3);
		GuoF[1]  = WConst[1] * (FU[0][0] + FC_3[0] - FdotU_3);
		GuoF[2]  = WConst[1] * (FU[0][0] - FC_3[0] - FdotU_3);
		GuoF[3]  = WConst[1] * (FU[1][1] + FC_3[1] - FdotU_3);
		GuoF[4]  = WConst[1] * (FU[1][1] - FC_3[1] - FdotU_3);
		GuoF[5]  = WConst[1] * (FU[2][2] + FC_3[2] - FdotU_3);
		GuoF[6]  = WConst[1] * (FU[2][2] - FC_3[2] - FdotU_3);
		GuoF[7]  = WConst[2] * (FU[0][0] + FU[0][1] + FC_3[0] + FU[1][0] + FU[1][1] + FC_3[1] - FdotU_3);
		GuoF[8]  = WConst[2] * (FU[0][0] - FU[0][1] - FC_3[0] - FU[1][0] + FU[1][1] + FC_3[1] - FdotU_3);
		GuoF[9]  = WConst[2] * (FU[0][0] - FU[0][1] + FC_3[0] - FU[1][0] + FU[1][1] - FC_3[1] - FdotU_3);
		GuoF[10] = WConst[2] * (FU[0][0] + FU[0][1] - FC_3[0] + FU[1][0] + FU[1][1] - FC_3[1] - FdotU_3);
		GuoF[11] = WConst[2] * (FU[1][1] + FU[1][2] + FC_3[1] + FU[2][1] + FU[2][2] + FC_3[2] - FdotU_3);
		GuoF[12] = WConst[2] * (FU[1][1] - FU[1][2] - FC_3[1] - FU[2][1] + FU[2][2] + FC_3[2] - FdotU_3);
		GuoF[13] = WConst[2] * (FU[1][1] - FU[1][2] + FC_3[1] - FU[2][1] + FU[2][2] - FC_3[2] - FdotU_3);
		GuoF[14] = WConst[2] * (FU[1][1] + FU[1][2] - FC_3[1] + FU[2][1] + FU[2][2] - FC_3[2] - FdotU_3);
		GuoF[15] = WConst[2] * (FU[0][0] + FU[0][2] + FC_3[0] + FU[2][0] + FU[2][2] + FC_3[2] - FdotU_3);
		GuoF[16] = WConst[2] * (FU[0][0] - FU[0][2] - FC_3[0] - FU[2][0] + FU[2][2] + FC_3[2] - FdotU_3);
		GuoF[17] = WConst[2] * (FU[0][0] - FU[0][2] + FC_3[0] - FU[2][0] + FU[2][2] - FC_3[2] - FdotU_3);
		GuoF[18] = WConst[2] * (FU[0][0] + FU[0][2] - FC_3[0] + FU[2][0] + FU[2][2] - FC_3[2] - FdotU_3);

		//

		//Calculate (fEq - f) = equilibrium distribution minus node distribution

		double fEq_f[19];

		fEq_f[0]  = fEquilibrium0 (Density, USq, u) - f[0];
		fEq_f[1]  = fEquilibrium1 (Density, USq, u) - f[1];
		fEq_f[2]  = fEquilibrium2 (Density, USq, u) - f[2];
		fEq_f[3]  = fEquilibrium3 (Density, USq, u) - f[3];
		fEq_f[4]  = fEquilibrium4 (Density, USq, u) - f[4];
		fEq_f[5]  = fEquilibrium5 (Density, USq, u) - f[5];
		fEq_f[6]  = fEquilibrium6 (Density, USq, u) - f[6];
		fEq_f[7]  = fEquilibrium7 (Density, USq, u) - f[7];
		fEq_f[8]  = fEquilibrium8 (Density, USq, u) - f[8];
		fEq_f[9]  = fEquilibrium9 (Density, USq, u) - f[9];
		fEq_f[10] = fEquilibrium10(Density, USq, u) - f[10];
		fEq_f[11] = fEquilibrium11(Density, USq, u) - f[11];
		fEq_f[12] = fEquilibrium12(Density, USq, u) - f[12];
		fEq_f[13] = fEquilibrium13(Density, USq, u) - f[13];
		fEq_f[14] = fEquilibrium14(Density, USq, u) - f[14];
		fEq_f[15] = fEquilibrium15(Density, USq, u) - f[15];
		fEq_f[16] = fEquilibrium16(Density, USq, u) - f[16];
		fEq_f[17] = fEquilibrium17(Density, USq, u) - f[17];
		fEq_f[18] = fEquilibrium18(Density, USq, u) - f[18];

		//

		//Calculate Collision Matrix * (fEq - f) and Collision Matrix * GuoF

			//Precalculate sums for entries commonly appearing with the same factor

		double fEq_f_3_6   = fEq_f[3]  + fEq_f[4]  + fEq_f[5]  + fEq_f[6];		//3 -> 6		(5x)
		double fEq_f_1_6   = fEq_f[1]  + fEq_f[2]  + fEq_f_3_6;					//1 -> 6		(3x)
		double fEq_f_7_10  = fEq_f[7]  + fEq_f[8]  + fEq_f[9]  + fEq_f[10];		//7 -> 10		(7x)
		double fEq_f_11_14 = fEq_f[11] + fEq_f[12] + fEq_f[13] + fEq_f[14];		//11 -> 14		(5x)
		double fEq_f_15_18 = fEq_f[15] + fEq_f[16] + fEq_f[17] + fEq_f[18];		//15 -> 18		(5x)
		double fEq_f_11_18 = fEq_f_11_14 + fEq_f_15_18;							//11 -> 18		(3x)
		double fEq_f_7_18  = fEq_f_7_10  + fEq_f_11_18;							//7 -> 18		(3x)


		double MfEq_f[19];		//Matrix * (fEq - f)

		MfEq_f[0]  =       fEq_f[0] +      fEq_f_1_6 +     fEq_f_7_18;
		MfEq_f[1]  = -30 * fEq_f[0] - 11 * fEq_f_1_6 + 8 * fEq_f_7_18;
		MfEq_f[2]  =  12 * fEq_f[0] - 4  * fEq_f_1_6 +     fEq_f_7_18;
		MfEq_f[3]  =       fEq_f[1]  - fEq_f[2]  + fEq_f[7]  - fEq_f[8]  + fEq_f[9]  - fEq_f[10] + fEq_f[15] - fEq_f[16] + fEq_f[17] - fEq_f[18];
		MfEq_f[4]  = 4 * (-fEq_f[1]  + fEq_f[2]) + fEq_f[7]  - fEq_f[8]  + fEq_f[9]  - fEq_f[10] + fEq_f[15] - fEq_f[16] + fEq_f[17] - fEq_f[18];
		MfEq_f[5]  =       fEq_f[3]  - fEq_f[4]  + fEq_f[7]  + fEq_f[8]  - fEq_f[9]  - fEq_f[10] + fEq_f[11] - fEq_f[12] + fEq_f[13] - fEq_f[14];
		MfEq_f[6]  = 4 * (-fEq_f[3]  + fEq_f[4]) + fEq_f[7]  + fEq_f[8]  - fEq_f[9]  - fEq_f[10] + fEq_f[11] - fEq_f[12] + fEq_f[13] - fEq_f[14];
		MfEq_f[7]  =       fEq_f[5]  - fEq_f[6]  + fEq_f[11] + fEq_f[12] - fEq_f[13] - fEq_f[14] + fEq_f[15] + fEq_f[16] - fEq_f[17] - fEq_f[18];
		MfEq_f[8]  = 4 * (-fEq_f[5]  + fEq_f[6]) + fEq_f[11] + fEq_f[12] - fEq_f[13] - fEq_f[14] + fEq_f[15] + fEq_f[16] - fEq_f[17] - fEq_f[18];
		MfEq_f[9]  = 2 * ( fEq_f[1]  + fEq_f[2]) -   fEq_f_3_6 + fEq_f_7_10 - 2*fEq_f_11_14 + fEq_f_15_18;
		MfEq_f[10] =-4 * ( fEq_f[1]  + fEq_f[2]) + 2*fEq_f_3_6 + fEq_f_7_10 - 2*fEq_f_11_14 + fEq_f_15_18;
		MfEq_f[11] =       fEq_f[3]  + fEq_f[4]  -    fEq_f[5] - fEq_f[6]  + fEq_f_7_10 - fEq_f_15_18;
		MfEq_f[12] =-2 * ( fEq_f[3]  + fEq_f[4]) + 2*(fEq_f[5] + fEq_f[6]) + fEq_f_7_10 - fEq_f_15_18;
		MfEq_f[13] =       fEq_f[7]  - fEq_f[8]  - fEq_f[9]  + fEq_f[10];
		MfEq_f[14] =       fEq_f[11] - fEq_f[12] - fEq_f[13] + fEq_f[14];
		MfEq_f[15] =       fEq_f[15] - fEq_f[16] - fEq_f[17] + fEq_f[18];
		MfEq_f[16] =       fEq_f[7]  - fEq_f[8]  + fEq_f[9]  - fEq_f[10] - fEq_f[15] + fEq_f[16] - fEq_f[17] + fEq_f[18];
		MfEq_f[17] =      -fEq_f[7]  - fEq_f[8]  + fEq_f[9]  + fEq_f[10] + fEq_f[11] - fEq_f[12] + fEq_f[13] - fEq_f[14];
		MfEq_f[18] =      -fEq_f[11] - fEq_f[12] + fEq_f[13] + fEq_f[14] + fEq_f[15] + fEq_f[16] - fEq_f[17] - fEq_f[18];
	
			//Precalculate sums for entries commonly appearing with the same factor

		double GuoF_3_6   = GuoF[3]  + GuoF[4]  + GuoF[5]  + GuoF[6];		//3 -> 6		(5x)
		double GuoF_1_6   = GuoF[1]  + GuoF[2]  + GuoF_3_6;					//1 -> 6		(3x)
		double GuoF_7_10  = GuoF[7]  + GuoF[8]  + GuoF[9]  + GuoF[10];		//7 -> 10		(7x)
		double GuoF_11_14 = GuoF[11] + GuoF[12] + GuoF[13] + GuoF[14];		//11 -> 14		(5x)
		double GuoF_15_18 = GuoF[15] + GuoF[16] + GuoF[17] + GuoF[18];		//15 -> 18		(5x)
		double GuoF_11_18 = GuoF_11_14 + GuoF_15_18;						//11 -> 18		(3x)
		double GuoF_7_18  = GuoF_7_10  + GuoF_11_18;						//7 -> 18		(3x)


		double MGuoF[19];		//Matrix * GuoF

		MGuoF[0]  =       GuoF[0] +      GuoF_1_6 +     GuoF_7_18;
		MGuoF[1]  = -30 * GuoF[0] - 11 * GuoF_1_6 + 8 * GuoF_7_18;
		MGuoF[2]  =  12 * GuoF[0] - 4  * GuoF_1_6 +     GuoF_7_18;
		MGuoF[3]  =       GuoF[1]  - GuoF[2]  + GuoF[7]  - GuoF[8]  + GuoF[9]  - GuoF[10] + GuoF[15] - GuoF[16] + GuoF[17] - GuoF[18];
		MGuoF[4]  = 4 * (-GuoF[1]  + GuoF[2]) + GuoF[7]  - GuoF[8]  + GuoF[9]  - GuoF[10] + GuoF[15] - GuoF[16] + GuoF[17] - GuoF[18];
		MGuoF[5]  =       GuoF[3]  - GuoF[4]  + GuoF[7]  + GuoF[8]  - GuoF[9]  - GuoF[10] + GuoF[11] - GuoF[12] + GuoF[13] - GuoF[14];
		MGuoF[6]  = 4 * (-GuoF[3]  + GuoF[4]) + GuoF[7]  + GuoF[8]  - GuoF[ 9]  - GuoF[10] + GuoF[11] - GuoF[12] + GuoF[13] - GuoF[14];
		MGuoF[7]  =       GuoF[5]  - GuoF[6]  + GuoF[11] + GuoF[12] - GuoF[13] - GuoF[14] + GuoF[15] + GuoF[16] - GuoF[17] - GuoF[18];
		MGuoF[8]  = 4 * (-GuoF[5]  + GuoF[6]) + GuoF[11] + GuoF[12] - GuoF[13] - GuoF[14] + GuoF[15] + GuoF[16] - GuoF[17] - GuoF[18];
		MGuoF[9]  = 2 * ( GuoF[1]  + GuoF[2]) -   GuoF_3_6 + GuoF_7_10 - 2*GuoF_11_14 + GuoF_15_18;
		MGuoF[10] =-4 * ( GuoF[1]  + GuoF[2]) + 2*GuoF_3_6 + GuoF_7_10 - 2*GuoF_11_14 + GuoF_15_18;
		MGuoF[11] =       GuoF[3]  + GuoF[4]  -    GuoF[5] - GuoF[6]  + GuoF_7_10 - GuoF_15_18;
		MGuoF[12] =-2 * ( GuoF[3]  + GuoF[4]) + 2*(GuoF[5] + GuoF[6]) + GuoF_7_10 - GuoF_15_18;
		MGuoF[13] =       GuoF[7]  - GuoF[8]  - GuoF[9]  + GuoF[10];
		MGuoF[14] =       GuoF[11] - GuoF[12] - GuoF[13] + GuoF[14];
		MGuoF[15] =       GuoF[15] - GuoF[16] - GuoF[17] + GuoF[18];
		MGuoF[16] =       GuoF[7]  - GuoF[8]  + GuoF[9]  - GuoF[10] - GuoF[15] + GuoF[16] - GuoF[17] + GuoF[18];
		MGuoF[17] =      -GuoF[7]  - GuoF[8]  + GuoF[9]  + GuoF[10] + GuoF[11] - GuoF[12] + GuoF[13] - GuoF[14];
		MGuoF[18] =      -GuoF[11] - GuoF[12] + GuoF[13] + GuoF[14] + GuoF[15] + GuoF[16] - GuoF[17] - GuoF[18];
	
		//

		//Calculate SDiag*M*(fEq - f) + dtI_SDiag*M*GuoF where dtI_SDiag = dt*(I - 0.5*SDiag)

		double SMfEqG[19];

		SMfEqG[0]  = SDiag[0] *MfEq_f[0]  + dtI_SDiag[0] *MGuoF[0];
		SMfEqG[1]  = SDiag[1] *MfEq_f[1]  + dtI_SDiag[1] *MGuoF[1];
		SMfEqG[2]  = SDiag[2] *MfEq_f[2]  + dtI_SDiag[2] *MGuoF[2];
		SMfEqG[3]  = SDiag[3] *MfEq_f[3]  + dtI_SDiag[3] *MGuoF[3];
		SMfEqG[4]  = SDiag[4] *MfEq_f[4]  + dtI_SDiag[4] *MGuoF[4];
		SMfEqG[5]  = SDiag[5] *MfEq_f[5]  + dtI_SDiag[5] *MGuoF[5];
		SMfEqG[6]  = SDiag[6] *MfEq_f[6]  + dtI_SDiag[6] *MGuoF[6];
		SMfEqG[7]  = SDiag[7] *MfEq_f[7]  + dtI_SDiag[7] *MGuoF[7];
		SMfEqG[8]  = SDiag[8] *MfEq_f[8]  + dtI_SDiag[8] *MGuoF[8];
		SMfEqG[9]  = SDiag[9] *MfEq_f[9]  + dtI_SDiag[9] *MGuoF[9];
		SMfEqG[10] = SDiag[10]*MfEq_f[10] + dtI_SDiag[10]*MGuoF[10];
		SMfEqG[11] = SDiag[11]*MfEq_f[11] + dtI_SDiag[11]*MGuoF[11];
		SMfEqG[12] = SDiag[12]*MfEq_f[12] + dtI_SDiag[12]*MGuoF[12];
		SMfEqG[13] = SDiag[13]*MfEq_f[13] + dtI_SDiag[13]*MGuoF[13];
		SMfEqG[14] = SDiag[14]*MfEq_f[14] + dtI_SDiag[14]*MGuoF[14];
		SMfEqG[15] = SDiag[15]*MfEq_f[15] + dtI_SDiag[15]*MGuoF[15];
		SMfEqG[16] = SDiag[16]*MfEq_f[16] + dtI_SDiag[16]*MGuoF[16];
		SMfEqG[17] = SDiag[17]*MfEq_f[17] + dtI_SDiag[17]*MGuoF[17];
		SMfEqG[18] = SDiag[18]*MfEq_f[18] + dtI_SDiag[18]*MGuoF[18];


		//Calculate Collision operator + Force term as MInv * SMfEqG;

		double MInvC[19];

			//Few unique values in the inverse matrix, precalculate (30 multiplications)

		double MInvCV[19][3];			//No more than 3 unique values per column, ignoring 0

		MInvCV[0][0]  = ( 5.2631578947368432e-002 ) * SMfEqG[0];	//1 unique value in column  0
		MInvCV[1][0]  = ( 1.2531328320802001e-002 ) * SMfEqG[1];	//3 unique values in column 1
		MInvCV[1][1]  = ( 4.5948203842940691e-003 ) * SMfEqG[1];	
		MInvCV[1][2]  = ( 3.3416875522138691e-003 ) * SMfEqG[1];
		MInvCV[2][0]  = ( 4.7619047619047623e-002 ) * SMfEqG[2];	//3 unique values in column 2
		MInvCV[2][1]  = ( 1.5873015873015872e-002 ) * SMfEqG[2];	
		MInvCV[2][2]  = ( 3.9682539682539698e-003 ) * SMfEqG[2];
		MInvCV[3][0]  = ( 1.0000000000000001e-001 ) * SMfEqG[3];	//1 unique value in column  3
		MInvCV[4][0]  = ( 1.0000000000000001e-001 ) * SMfEqG[4];	//2 unique values in column 4
		MInvCV[4][1]  = ( 2.5000000000000001e-002 ) * SMfEqG[4];	
		MInvCV[5][0]  = ( 1.0000000000000001e-001 ) * SMfEqG[5];	//1 unique value in column  5
		MInvCV[6][0]  = ( 1.0000000000000001e-001 ) * SMfEqG[6];	//2 unique values in column 6
		MInvCV[6][1]  = ( 2.5000000000000001e-002 ) * SMfEqG[6];	
		MInvCV[7][0]  = ( 1.0000000000000001e-001 ) * SMfEqG[7];	//1 unique value in column  7
		MInvCV[8][0]  = ( 1.0000000000000001e-001 ) * SMfEqG[8];	//2 unique values in column 8
		MInvCV[8][1]  = ( 2.5000000000000001e-002 ) * SMfEqG[8];
		MInvCV[9][0]  = ( 5.5555555555555552e-002 ) * SMfEqG[9];	//2 unique values in column 9
		MInvCV[9][1]  = ( 2.7777777777777780e-002 ) * SMfEqG[9];
		MInvCV[10][0] = ( 5.5555555555555552e-002 ) * SMfEqG[10];	//3 unique values in column 10
		MInvCV[10][1] = ( 2.7777777777777780e-002 ) * SMfEqG[10];	
		MInvCV[10][2] = ( 1.3888888888888890e-002 ) * SMfEqG[10];
		MInvCV[11][0] = ( 8.3333333333333329e-002 ) * SMfEqG[11];	//1 unique value in column  11
		MInvCV[12][0] = ( 8.3333333333333343e-002 ) * SMfEqG[12];	//2 unique values in column 12
		MInvCV[12][1] = ( 4.1666666666666664e-002 ) * SMfEqG[12];
		MInvCV[13][0] = ( 2.5000000000000000e-001 ) * SMfEqG[13];	//1 unique value in column  13
		MInvCV[14][0] = ( 2.5000000000000000e-001 ) * SMfEqG[14];	//1 unique value in column  14
		MInvCV[15][0] = ( 2.5000000000000000e-001 ) * SMfEqG[15];	//1 unique value in column  15
		MInvCV[16][0] = ( 1.2500000000000000e-001 ) * SMfEqG[16];	//1 unique value in column  16
		MInvCV[17][0] = ( 1.2500000000000000e-001 ) * SMfEqG[17];	//1 unique value in column  17
		MInvCV[18][0] = ( 1.2500000000000000e-001 ) * SMfEqG[18];	//1 unique value in column  18

		//MInvC

		double MInvC_0_2[2];	//Two common sums across columns 0 to 2

		MInvC_0_2[0] =  MInvCV[0][0] - MInvCV[1][1] - MInvCV[2][1];		//Used by rows 1 -> 6
		MInvC_0_2[1] =  MInvCV[0][0] + MInvCV[1][2] + MInvCV[2][2];		//Used by rows 7 -> 18

		MInvC[0]  = MInvCV[0][0] - MInvCV[1][0] + MInvCV[2][0];
		MInvC[1]  = MInvC_0_2[0] + MInvCV[3][0] - MInvCV[4][0]                                                             + MInvCV[9][0] - MInvCV[10][0];
		MInvC[2]  = MInvC_0_2[0] - MInvCV[3][0] + MInvCV[4][0]                                                             + MInvCV[9][0] - MInvCV[10][0];
		MInvC[3]  = MInvC_0_2[0]                               + MInvCV[5][0] - MInvCV[6][0]                               - MInvCV[9][1] + MInvCV[10][1] + MInvCV[11][0] - MInvCV[12][0];
		MInvC[4]  = MInvC_0_2[0]                               - MInvCV[5][0] + MInvCV[6][0]                               - MInvCV[9][1] + MInvCV[10][1] + MInvCV[11][0] - MInvCV[12][0];
		MInvC[5]  = MInvC_0_2[0]                                                             + MInvCV[7][0] - MInvCV[8][0] - MInvCV[9][1] + MInvCV[10][1] - MInvCV[11][0] + MInvCV[12][0];
		MInvC[6]  = MInvC_0_2[0]                                                             - MInvCV[7][0] + MInvCV[8][0] - MInvCV[9][1] + MInvCV[10][1] - MInvCV[11][0] + MInvCV[12][0];
		MInvC[7]  = MInvC_0_2[1] + MInvCV[3][0] + MInvCV[4][1] + MInvCV[5][0] + MInvCV[6][1]                               + MInvCV[9][1] + MInvCV[10][2] + MInvCV[11][0] + MInvCV[12][1] + MInvCV[13][0]                                 + MInvCV[16][0] - MInvCV[17][0];
		MInvC[8]  = MInvC_0_2[1] - MInvCV[3][0] - MInvCV[4][1] + MInvCV[5][0] + MInvCV[6][1]                               + MInvCV[9][1] + MInvCV[10][2] + MInvCV[11][0] + MInvCV[12][1] - MInvCV[13][0]                                 - MInvCV[16][0] - MInvCV[17][0];
		MInvC[9]  = MInvC_0_2[1] + MInvCV[3][0] + MInvCV[4][1] - MInvCV[5][0] - MInvCV[6][1]                               + MInvCV[9][1] + MInvCV[10][2] + MInvCV[11][0] + MInvCV[12][1] - MInvCV[13][0]                                 + MInvCV[16][0] + MInvCV[17][0];
		MInvC[10] = MInvC_0_2[1] - MInvCV[3][0] - MInvCV[4][1] - MInvCV[5][0] - MInvCV[6][1]                               + MInvCV[9][1] + MInvCV[10][2] + MInvCV[11][0] + MInvCV[12][1] + MInvCV[13][0]                                 - MInvCV[16][0] + MInvCV[17][0];
		MInvC[11] = MInvC_0_2[1]                               + MInvCV[5][0] + MInvCV[6][1] + MInvCV[7][0] + MInvCV[8][1] - MInvCV[9][0] - MInvCV[10][1]                                                 + MInvCV[14][0]                                 + MInvCV[17][0] - MInvCV[18][0];
		MInvC[12] = MInvC_0_2[1]                               - MInvCV[5][0] - MInvCV[6][1] + MInvCV[7][0] + MInvCV[8][1] - MInvCV[9][0] - MInvCV[10][1]                                                 - MInvCV[14][0]                                 - MInvCV[17][0] - MInvCV[18][0];
		MInvC[13] = MInvC_0_2[1]                               + MInvCV[5][0] + MInvCV[6][1] - MInvCV[7][0] - MInvCV[8][1] - MInvCV[9][0] - MInvCV[10][1]                                                 - MInvCV[14][0]                                 + MInvCV[17][0] + MInvCV[18][0];
		MInvC[14] = MInvC_0_2[1]                               - MInvCV[5][0] - MInvCV[6][1] - MInvCV[7][0] - MInvCV[8][1] - MInvCV[9][0] - MInvCV[10][1]                                                 + MInvCV[14][0]                                 - MInvCV[17][0] + MInvCV[18][0];
		MInvC[15] = MInvC_0_2[1] + MInvCV[3][0] + MInvCV[4][1]                               + MInvCV[7][0] + MInvCV[8][1] + MInvCV[9][1] + MInvCV[10][2] - MInvCV[11][0] - MInvCV[12][1]                                 + MInvCV[15][0] - MInvCV[16][0]                 + MInvCV[18][0];
		MInvC[16] = MInvC_0_2[1] - MInvCV[3][0] - MInvCV[4][1]                               + MInvCV[7][0] + MInvCV[8][1] + MInvCV[9][1] + MInvCV[10][2] - MInvCV[11][0] - MInvCV[12][1]                                 - MInvCV[15][0] + MInvCV[16][0]                 + MInvCV[18][0];
		MInvC[17] = MInvC_0_2[1] + MInvCV[3][0] + MInvCV[4][1]                               - MInvCV[7][0] - MInvCV[8][1] + MInvCV[9][1] + MInvCV[10][2] - MInvCV[11][0] - MInvCV[12][1]                                 - MInvCV[15][0] - MInvCV[16][0]                 - MInvCV[18][0];
		MInvC[18] = MInvC_0_2[1] - MInvCV[3][0] - MInvCV[4][1]                               - MInvCV[7][0] - MInvCV[8][1] + MInvCV[9][1] + MInvCV[10][2] - MInvCV[11][0] - MInvCV[12][1]                                 + MInvCV[15][0] + MInvCV[16][0]                 - MInvCV[18][0];

		//

		//Finally, add original distribution to collision operator and force term to get new distribution
		//Thus new distribution = MInvC[i] + f[i]

		//Streaming and bounce-back

		F[0] = MInvC[0] + f[0];			//No streaming or bounce-back for centre node

		for (int mi=1; mi<19; mi++){

			Node* NeighbourNode = Nodes[i].Neighbour[mi];

			if(NeighbourNode == 0){									//Bounce-back condition	
				Nodes[i].fCalc[ReflMap[mi]] = MInvC[mi] + f[mi];
			}else{													//Stream to neighbour node
				NeighbourNode->fCalc[mi] = MInvC[mi] + f[mi];
			}

		}

	}

}

	//Equilibrium function for all components
inline double fEquilibrium(double Density, int Basis, double Velocity[3]){	//~12 double multiplications

	/*
	double cEq;				//c = dx/dt lattice velocity
	double cEqSq;			//c^2
	double Inv_cEqSq;		//1/(c^2)
	*/

	double uSq = Velocity[0]*Velocity[0] + Velocity[1]*Velocity[1] + Velocity[2]*Velocity[2];
	double eu = (BasisVector[Basis][0]*Velocity[0] + BasisVector[Basis][1]*Velocity[1] + BasisVector[Basis][2]*Velocity[2]);

	double f = Weighting[Basis] * Density * (1.0 + Inv_cEqSq*(eu*(3.0 + (9.0/2.0)*Inv_cEqSq*eu) - (3.0/2.0)*uSq));	//6 multiplications
	
	return f;
}

	//Individual equilibrium function for each component (for speed)
inline double fEquilibrium0(double Density, double uSq, double Velocity[3]){	//Centre ~3 double multiplications
	return  Density * ((Weighting[0]*1.0) - (Weighting[0]*3.0/2.0)*uSq*Inv_cEqSq);	
}
inline double fEquilibrium1(double Density, double uSq, double Velocity[3]){	//+X ~6 double multiplications
	double eu = Velocity[0];
	return Density * ((Weighting[1]*1.0) + Inv_cEqSq*(eu*((Weighting[1]*3.0) + (Weighting[1]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[1]*3.0/2.0)*uSq));	
}
inline double fEquilibrium2(double Density, double uSq, double Velocity[3]){	//-X ~6 double multiplications
	double eu = -Velocity[0];
	return Density * ((Weighting[2]*1.0) + Inv_cEqSq*(eu*((Weighting[2]*3.0) + (Weighting[2]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[2]*3.0/2.0)*uSq));	
}
inline double fEquilibrium3(double Density, double uSq, double Velocity[3]){	//+Y ~6 double multiplications
	double eu = Velocity[1];
	return Density * ((Weighting[3]*1.0) + Inv_cEqSq*(eu*((Weighting[3]*3.0) + (Weighting[3]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[3]*3.0/2.0)*uSq));	
}
inline double fEquilibrium4(double Density, double uSq, double Velocity[3]){	//-Y ~6 double multiplications
	double eu = -Velocity[1];
	return Density * ((Weighting[4]*1.0) + Inv_cEqSq*(eu*((Weighting[4]*3.0) + (Weighting[4]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[4]*3.0/2.0)*uSq));	
}
inline double fEquilibrium5(double Density, double uSq, double Velocity[3]){	//+Z ~6 double multiplications
	double eu = Velocity[2];
	return Density * ((Weighting[5]*1.0) + Inv_cEqSq*(eu*((Weighting[5]*3.0) + (Weighting[5]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[5]*3.0/2.0)*uSq));	
}
inline double fEquilibrium6(double Density, double uSq, double Velocity[3]){	//-Z ~6 double multiplications
	double eu = -Velocity[2];
	return Density * ((Weighting[6]*1.0) + Inv_cEqSq*(eu*((Weighting[6]*3.0) + (Weighting[6]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[6]*3.0/2.0)*uSq));	
}
inline double fEquilibrium7(double Density, double uSq, double Velocity[3]){	//+X +Y ~6 double multiplications
	double eu = Velocity[0] + Velocity[1];
	return Density * ((Weighting[7]*1.0) + Inv_cEqSq*(eu*((Weighting[7]*3.0) + (Weighting[7]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[7]*3.0/2.0)*uSq));	
}
inline double fEquilibrium8(double Density, double uSq, double Velocity[3]){	//-X +Y ~6 double multiplications
	double eu = -Velocity[0] + Velocity[1];
	return Density * ((Weighting[8]*1.0) + Inv_cEqSq*(eu*((Weighting[8]*3.0) + (Weighting[8]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[8]*3.0/2.0)*uSq));	
}
inline double fEquilibrium9(double Density, double uSq, double Velocity[3]){	//+X -Y ~6 double multiplications
	double eu = Velocity[0] - Velocity[1];
	return Density * ((Weighting[9]*1.0) + Inv_cEqSq*(eu*((Weighting[9]*3.0) + (Weighting[9]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[9]*3.0/2.0)*uSq));	
}
inline double fEquilibrium10(double Density, double uSq, double Velocity[3]){	//-X -Y ~6 double multiplications
	double eu = -Velocity[0] - Velocity[1];
	return Density * ((Weighting[10]*1.0) + Inv_cEqSq*(eu*((Weighting[10]*3.0) + (Weighting[10]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[10]*3.0/2.0)*uSq));	
}
inline double fEquilibrium11(double Density, double uSq, double Velocity[3]){	//+Y +Z ~6 double multiplications
	double eu = Velocity[1] + Velocity[2];
	return Density * ((Weighting[11]*1.0) + Inv_cEqSq*(eu*((Weighting[11]*3.0) + (Weighting[11]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[11]*3.0/2.0)*uSq));	
}
inline double fEquilibrium12(double Density, double uSq, double Velocity[3]){	//-Y +Z ~6 double multiplications
	double eu = -Velocity[1] + Velocity[2];
	return Density * ((Weighting[12]*1.0) + Inv_cEqSq*(eu*((Weighting[12]*3.0) + (Weighting[12]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[12]*3.0/2.0)*uSq));	
}
inline double fEquilibrium13(double Density, double uSq, double Velocity[3]){	//+Y -Z ~6 double multiplications
	double eu = Velocity[1] - Velocity[2];
	return Density * ((Weighting[13]*1.0) + Inv_cEqSq*(eu*((Weighting[13]*3.0) + (Weighting[13]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[13]*3.0/2.0)*uSq));	
}
inline double fEquilibrium14(double Density, double uSq, double Velocity[3]){	//-Y -Z ~6 double multiplications
	double eu = -Velocity[1] - Velocity[2];
	return Density * ((Weighting[14]*1.0) + Inv_cEqSq*(eu*((Weighting[14]*3.0) + (Weighting[14]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[14]*3.0/2.0)*uSq));	
}
inline double fEquilibrium15(double Density, double uSq, double Velocity[3]){	//+X +Z ~6 double multiplications
	double eu = Velocity[0] + Velocity[2];
	return Density * ((Weighting[15]*1.0) + Inv_cEqSq*(eu*((Weighting[15]*3.0) + (Weighting[15]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[15]*3.0/2.0)*uSq));	
}
inline double fEquilibrium16(double Density, double uSq, double Velocity[3]){	//-X +Z ~6 double multiplications
	double eu = -Velocity[0] + Velocity[2];
	return Density * ((Weighting[16]*1.0) + Inv_cEqSq*(eu*((Weighting[16]*3.0) + (Weighting[16]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[16]*3.0/2.0)*uSq));	
}
inline double fEquilibrium17(double Density, double uSq, double Velocity[3]){	//+X -Z ~6 double multiplications
	double eu = Velocity[0] - Velocity[2];
	return Density * ((Weighting[17]*1.0) + Inv_cEqSq*(eu*((Weighting[17]*3.0) + (Weighting[17]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[17]*3.0/2.0)*uSq));	
}
inline double fEquilibrium18(double Density, double uSq, double Velocity[3]){	//-X -Z ~6 double multiplications
	double eu = -Velocity[0] - Velocity[2];
	return Density * ((Weighting[18]*1.0) + Inv_cEqSq*(eu*((Weighting[18]*3.0) + (Weighting[18]*9.0/2.0)*Inv_cEqSq*eu) - (Weighting[18]*3.0/2.0)*uSq));	
}


bool Initialise(){

	nVoxels = NLatticeX * NLatticeY * NLatticeZ;	//Total number of voxels

	Solid = new long long[nVoxels];					//Array of element indices or solid = -1

	if(Solid==0){
		cout << "Unable to allocate memory (0)" << endl;
		return false;
	}
	
	bool ReadStat;
	if(GeometryFileBin){
		ReadStat = ReadGeometryBin(&nNonSolid);
	}else{
		ReadStat = ReadGeometry(&nNonSolid);
	}

	if(!ReadStat){									//File read in failed
		return false;
	}


	Nodes = new Node[nNonSolid+1];					//Create node data. Node[0] = Solid

	if(Nodes==0){
		cout << "Unable to allocate memory (1)" << endl;
		return false;
	}

	double dx = 1.0;	//Lattice quantities: dx = grid spacing
	double dt = 1.0;	//dt = time-step. dx/dt = Lattice velocity

	//Quantities related to the lattice velocity constant 'c', precalculated for speed
	cEq = dx/dt;						//c = sqrt(3RT); R=8.31, T=Temp or microscopic lattice velocity
	cEqSq = cEq * cEq;					//c^2
	Inv_cEqSq = 1.0/cEqSq;				//1/(c^2)
	cEq_3 = cEq / 3.0;					//c/3
	//

	//Initialise diagonal of relaxation time matrix S

	double SDiag1 = 1.0 / ( (3.0*Viscosity*Inv_cEqSq/dt) + 0.5 );	//Relaxation times
	double SDiag2 = 8.0*(2.0-SDiag1)/(8.0-SDiag1);

	SDiag[0]  = 0;
	SDiag[1]  = SDiag1;
	SDiag[2]  = SDiag1;
	SDiag[3]  = 0;
	SDiag[4]  = SDiag2;
	SDiag[5]  = 0;
	SDiag[6]  = SDiag2;
	SDiag[7]  = 0;
	SDiag[8]  = SDiag2;
	SDiag[9]  = SDiag1;
	SDiag[10] = SDiag1;
	SDiag[11] = SDiag1;
	SDiag[12] = SDiag1;
	SDiag[13] = SDiag1;
	SDiag[14] = SDiag1;
	SDiag[15] = SDiag1;
	SDiag[16] = SDiag2;
	SDiag[17] = SDiag2;
	SDiag[18] = SDiag2;

	//Initialise dtI_SDiag matrix diagonal = dt * (I - 0.5S)

	for(int i=0; i!=19; i++){		
		dtI_SDiag[i] = dt * (1 - 0.5*SDiag[i]);
	}
	

	//

	//Initialise lattice

	int c = 0;

	for(int z=0; z<NLatticeZ; z++){
	for(int y=0; y<NLatticeY; y++){
	for(int x=0; x<NLatticeX; x++){

		if(Solid[c] == 0){		//Solid
			c++;
			continue;
		}

		long long i = Solid[c];		//Index in nodes array

		Nodes[i].Density = 1.0;
		
		Nodes[i].Velocity[0] = InitialVx;
		Nodes[i].Velocity[1] = InitialVy;
		Nodes[i].Velocity[2] = InitialVz;

		for(int i2=0; i2!=19; i2++){		//Calculate distribution function from initial velocity
			Nodes[i].f[i2] = fEquilibrium(Nodes[i].Density, i2, Nodes[i].Velocity);
		}

		//Find neighbouring nodes
			
		int xp = (x==NLatticeX-1) ? 0 : x+1;	//Loop boundaries
		int xn = (x==0) ? NLatticeX-1 : x-1;
		int yp = (y==NLatticeY-1) ? 0 : y+1;
		int yn = (y==0) ? NLatticeY-1 : y-1;
		int zp = (z==NLatticeZ-1) ? 0 : z+1;
		int zn = (z==0) ? NLatticeZ-1 : z-1;

			//If BoundaryCondition==1 (Solid boundary) and neighbour over loop boundary, set solid

		Nodes[i].Neighbour[1]  = (BoundaryConditionX==1&&xp==0)	? 0 : &Nodes[GetSolid(xp,y,z)];	//+X
		Nodes[i].Neighbour[2]  = (BoundaryConditionX==1&&x==0)	? 0 : &Nodes[GetSolid(xn,y,z)];	//-X
		Nodes[i].Neighbour[3]  = (BoundaryConditionY==1&&yp==0)	? 0 : &Nodes[GetSolid(x,yp,z)];	//+Y
		Nodes[i].Neighbour[4]  = (BoundaryConditionY==1&&y==0)	? 0 : &Nodes[GetSolid(x,yn,z)];	//-Y
		Nodes[i].Neighbour[5]  = (BoundaryConditionZ==1&&zp==0)	? 0 : &Nodes[GetSolid(x,y,zp)];	//+Z
		Nodes[i].Neighbour[6]  = (BoundaryConditionZ==1&&z==0)	? 0 : &Nodes[GetSolid(x,y,zn)];	//-Z

		Nodes[i].Neighbour[7]  = ((BoundaryConditionX==1&&xp==0)||(BoundaryConditionY==1&&yp==0)) ? 0 : &Nodes[GetSolid(xp,yp,z)];	//+X +Y
		Nodes[i].Neighbour[8]  = ((BoundaryConditionX==1&&x==0) ||(BoundaryConditionY==1&&yp==0)) ? 0 : &Nodes[GetSolid(xn,yp,z)];	//-X +Y
		Nodes[i].Neighbour[9]  = ((BoundaryConditionX==1&&xp==0)||(BoundaryConditionY==1&&y==0) ) ? 0 : &Nodes[GetSolid(xp,yn,z)];	//+X -Y
		Nodes[i].Neighbour[10] = ((BoundaryConditionX==1&&x==0) ||(BoundaryConditionY==1&&y==0) ) ? 0 : &Nodes[GetSolid(xn,yn,z)];	//-X -Y

		Nodes[i].Neighbour[11] = ((BoundaryConditionY==1&&yp==0)||(BoundaryConditionZ==1&&zp==0)) ? 0 : &Nodes[GetSolid(x,yp,zp)];	//   +Y +Z
		Nodes[i].Neighbour[12] = ((BoundaryConditionY==1&&y==0) ||(BoundaryConditionZ==1&&zp==0)) ? 0 : &Nodes[GetSolid(x,yn,zp)];	//   -Y +Z
		Nodes[i].Neighbour[13] = ((BoundaryConditionY==1&&yp==0)||(BoundaryConditionZ==1&&z==0) ) ? 0 : &Nodes[GetSolid(x,yp,zn)];	//   +Y -Z
		Nodes[i].Neighbour[14] = ((BoundaryConditionY==1&&y==0) ||(BoundaryConditionZ==1&&z==0) ) ? 0 : &Nodes[GetSolid(x,yn,zn)];	//   -Y -Z

		Nodes[i].Neighbour[15] = ((BoundaryConditionX==1&&xp==0)||(BoundaryConditionZ==1&&zp==0)) ? 0 : &Nodes[GetSolid(xp,y,zp)];	//+X    +Z
		Nodes[i].Neighbour[16] = ((BoundaryConditionX==1&&x==0) ||(BoundaryConditionZ==1&&zp==0)) ? 0 : &Nodes[GetSolid(xn,y,zp)];	//-X    +Z
		Nodes[i].Neighbour[17] = ((BoundaryConditionX==1&&xp==0)||(BoundaryConditionZ==1&&z==0) ) ? 0 : &Nodes[GetSolid(xp,y,zn)];	//+X    -Z
		Nodes[i].Neighbour[18] = ((BoundaryConditionX==1&&x==0) ||(BoundaryConditionZ==1&&z==0) ) ? 0 : &Nodes[GetSolid(xn,y,zn)];	//-X    -Z
		
			//Finally, set any &Nodes[0] pointers to 0 to indicate solid

		for(int i2=0; i2!=19; i2++){
			if(Nodes[i].Neighbour[i2] == &Nodes[0]){
				Nodes[i].Neighbour[i2] = 0;
			}
		}

		//

		c++;
	}
	}
	}

	return true;
}

void Finalise(){

	delete[] Solid;
	delete[] Nodes;

}

void OutputVelocities(int nTimeSteps){

	FILE* OutFile = fopen(OutputVelocityFile, "wb");

	if(!OutFile){
		return;
	}

	if(OutputVelocityFileVTKHeader && !OutputVelocityFileBin && !OutputVelocityFileSparse){

		char Header[2048];
		int Ind = 0;

		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "# vtk DataFile Version 2.0\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "Geometry Output\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ASCII\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DATASET STRUCTURED_POINTS\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DIMENSIONS %i %i %i\n", NLatticeX, NLatticeY, NLatticeZ);
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ORIGIN 0 0 0\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SPACING 1 1 1\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "POINT_DATA %i\n", nVoxels);
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "VECTORS sample_scalars float\n");

		fwrite(Header, sizeof(char), Ind, OutFile);

	}
	
	for(int z=0; z<NLatticeZ; z++){
	for(int y=0; y<NLatticeY; y++){
	for(int x=0; x<NLatticeX; x++){

		long long Index = GetSolid(x, y, z);

		if(OutputVelocityFileSparse && Index == 0){		//Solid
			continue;
		}

		double V[3];

		if(Index == 0){
			V[0] = 0;
			V[1] = 0;
			V[2] = 0;
		}else{
			V[0] = Nodes[Index].Velocity[0];
			V[1] = Nodes[Index].Velocity[1];
			V[2] = Nodes[Index].Velocity[2];
		}

		if(OutputVelocityFileBin){

			fwrite(V, sizeof(double), 3, OutFile);

		}else{

			char Str[128];
			int l;

			if(V[0]==0 && V[1]==0 && V[2]==0){
				l = _snprintf(Str, sizeof(Str), "0 0 0\n");
			}else{
				l = _snprintf(Str, sizeof(Str), "%e %e %e\n", V[0], V[1], V[2]);
			}

			fwrite(Str, sizeof(char), l, OutFile);

		}

	}
	}
	}

	fclose(OutFile);

}

void OutputGeometry(){

	FILE* OutFile = fopen(OutputGeometryFile, "wb");

	if(!OutFile){
		return;
	}

	if(OutputGeometryFileVTKHeader){

		char Header[2048];
		int Ind = 0;

		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "# vtk DataFile Version 2.0\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "Geometry Output\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ASCII\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DATASET STRUCTURED_POINTS\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DIMENSIONS %i %i %i\n", NLatticeX, NLatticeY, NLatticeZ);
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ORIGIN 0 0 0\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SPACING 1 1 1\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "POINT_DATA %i\n", nVoxels);
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SCALARS sample_scalars float\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "LOOKUP_TABLE default\n");

		fwrite(Header, sizeof(char), Ind, OutFile);

	}

	char c0[2] = {'0', ' '};
	char c1[2] = {'1', ' '};

	for(int z=0; z<NLatticeZ; z++){
	for(int y=0; y<NLatticeY; y++){
	for(int x=0; x<NLatticeX; x++){

		long long Index = GetSolid(x, y, z);

		if(OutputGeometryFileBin){

			bool b = (Index==0);
			fwrite(&b, sizeof(bool), 1, OutFile);

		}else{

			if(Index==0){	//Solid
				fwrite(c1, sizeof(char), 2, OutFile);
			}else{
				fwrite(c0, sizeof(char), 2, OutFile);
			}

		}

	}
	}
	}

	fclose(OutFile);
	
}

bool ReadGeometryBin(long long* NumberOfNonSolids){

	FILE* InFile = fopen(GeometryFile, "rb");
	
	if(InFile!=0){
		cout << "Reading in solids from '" << GeometryFile << "'" << endl;
	}else{
		cout << "Could not open solids file" << endl;
		return false;
	}

	bool ret = true;

	long long flength;

#ifdef _WIN32
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif
	
	const long long buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];
	
	int z=0;
	int y=0;
	int x=0;
	
	long long nRead = 0;				//Number of chars read from file in total
	long long readl = 0;				//Number of chars read from file into buffer
	long long readpos = 0;			//Read position in buffer

	int nNonSolids = 0;
	long long Index = 1;

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				cout << "File read error: reached end of geometry file before lattice entirely read" << endl;
				ret = false;
				goto ReadSolidsBinReturnRJMP;
			}
			fread(buf, 1, readl, InFile);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]){		//Solid

			SetSolid(x, y, z, 0);

		}else{					//Non solid

			SetSolid(x, y, z, Index);
			Index++;

			nNonSolids++;
		}

		x++;
		if(x==NLatticeX){
			x=0;
			y++;
			if(y==NLatticeY){
				y=0;
				z++;
				if(z==NLatticeZ){
					goto ReadSolidsBinReturnRJMP;
				}
			}
		}

		readpos++;
	}

ReadSolidsBinReturnRJMP:
	*NumberOfNonSolids = nNonSolids;
	delete[] buf;
	fclose(InFile);
	return ret;

}

bool ReadGeometry(long long* NumberOfNonSolids){

	FILE* InFile = fopen(GeometryFile, "rb");
	
	if(InFile!=0){
		cout << "Reading in solids from '" << GeometryFile << "'" << endl;
	}else{
		cout << "Could not open solids file" << endl;
		return false;
	}

	bool ret = true;
	long long flength;

#ifdef _WIN32
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif
	
	const long long buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	
	long long nRead = 0;				//Number of chars read from file in total
	long long readl = 0;				//Number of chars read from file into buffer
	long long readpos = 0;				//Read position in buffer

	int nNonSolids = 0;

	int NEndLinesRead = 0;
	long long Index = 1;

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				cout << "File read error: reached end of geometry file before lattice entirely read" << endl;
				ret = false;
				goto ReadSolidsReturnRJMP;
			}
			fread(buf, 1, readl, InFile);
			nRead += readl;
			readpos = 0;
		}

		if(GeometryFileVTKHeader && NEndLinesRead!=NLinesGeometryVTKHeader){
			if(buf[readpos]=='\n'){
				NEndLinesRead++;
			}
			readpos++;
			continue;
		}

		if(buf[readpos]=='1'||buf[readpos]=='0'){

			if(buf[readpos]=='1'){	//Solid

				SetSolid(x, y, z, 0);

			}else{					//Non solid

				SetSolid(x, y, z, Index);
				Index++;

				nNonSolids++;
			}

			x++;
			if(x==NLatticeX){
				x=0;
				y++;
				if(y==NLatticeY){
					y=0;
					z++;
					if(z==NLatticeZ){
						goto ReadSolidsReturnRJMP;
					}
				}
			}

		}

		readpos++;
	}

ReadSolidsReturnRJMP:
	*NumberOfNonSolids = nNonSolids;
	delete[] buf;
	fclose(InFile);
	return ret;
}

//Input file reading

int ReadProperty(ifstream& File, char* Property, char* Value){

	int PropertyL = 0;
	int ValueL = 0;

	char c;
	int Stage = 0;					//Reading through property (0), value (1), or comment (2); Invalid line (-1)
	bool InQuotes = false;			//Reading a value in quotes (keep white spaces)
	while(true){

		File >> c;

		if(File.eof()){				//Reached end of file
			if(Stage==0 || Stage==-1){
				return 3;			//Reached end of file with no value
			}
			Property[PropertyL] = '\0';
			Value[ValueL] = '\0';
			return 1;				//Reached end of file 
		}

		if(c=='\r' || c=='\t'){ continue; }
		if(c=='\n'){
			if(Stage==0 || Stage==-1){
				return 2;			//Reached end of invalid line with no value
			}
			Property[PropertyL] = '\0';
			Value[ValueL] = '\0';
			return 0;				//Finished 
		}
		if(c=='/' && (Stage==0||Stage==1)){
			c = ' ';
			File >> c;
			if(c=='/'){				//Reached comment
				if(Stage == 0){		//Invalid line
					Stage = -1;
				}else{				//End of value
					Stage = 2;
				}
			}else{
				c = '/';
				File.unget();
			}
		}

		if(Stage == 0){				//Reading in property name

			if(c==' '){ continue; }	//Ignore white-space

			if(c=='='){
				Stage = 1;			//Move onto value
			}else{
				Property[PropertyL] = tolower(c);
				PropertyL++;
			}

		}else if(Stage == 1){		//Reading in value

			if(!InQuotes && c==' '){ continue; }	//Ignore white-space outside quotes

			if(c=='"'){
				if(!InQuotes){
					InQuotes = true;
				}else{
					InQuotes = false;
					Stage = 2;
				}
			}else{
				Value[ValueL] = c;
				ValueL++;
			}

		}

	}

}

bool CompareStr(char* DataName, char* Property){
	int i=0;
	while(true){
		if(DataName[i]=='\0' && Property[i]=='\0'){
			return true;
		}

		if(tolower(DataName[i])!=Property[i]){
			return false;
		}
		i++;
	}
}

bool ReadBool(char* Value){

	if(Value[1]=='\0'){

		if(Value[0]=='0'){ return false; }
		if(Value[0]=='1'){ return true; }

		if(Value[0]=='n'){ return false; }
		if(Value[0]=='y'){ return true; }

		return false;
	}
	
	if(CompareStr("true",Value)){ return true; }
	if(CompareStr("false",Value)){ return false; }

	if(CompareStr("yes",Value)){ return true; }
	if(CompareStr("no",Value)){ return false; }

	if(CompareStr("on",Value)){ return true; }
	if(CompareStr("off",Value)){ return false; }

	return false;
}

int ReadInputFile(char* InputFileName){

	const int InDataLength = sizeof(InputFileData)/sizeof(InputValueStruct);

	bool InDataSet[InDataLength];
	memset(&InDataSet,false,InDataLength);

	ifstream InFile;
	InFile.open(InputFileName);
	InFile.unsetf(ios::skipws);

	if(InFile.fail()){
		return 1;			//File couldn't be opened
	}

	char Property[256];
	char Value[1024];

	while(true){

		int r = ReadProperty(InFile, Property, Value);

		if(r==0 || r==1){	//Valid input line

			int i=0;
			while(i!=InDataLength){
				if(CompareStr(InputFileData[i].DataName, Property)){

				switch(InputFileData[i].DataType){

					case DataType_Short:			//Interpret as short
						*((short*)(InputFileData[i].Ptr)) = (short)atoi(Value);
						break;
					case DataType_Bool:				//Interpret as bool
						*((bool*)(InputFileData[i].Ptr)) = ReadBool(Value);
						break;
					case DataType_Int:				//Interpret as int
						*((int*)(InputFileData[i].Ptr)) = atoi(Value);
						break;
					case DataType_Double:			//Interpret as double
						*((double*)(InputFileData[i].Ptr)) = atof(Value);
						break;
					case DataType_String:			//Copy over string
						memcpy(InputFileData[i].Ptr, Value, min(InputFileData[i].PtrLength,(int)sizeof(Value)));
						((char*)(InputFileData[i].Ptr))[InputFileData[i].PtrLength-1] = '\0';
						break;
				}

				InDataSet[i] = true;
				}
				i++;
			}

		}

		if(r==1 || r==3){	//End of file
			break;
		}
	
	}

	int nRequired = 0;

	int i=0;
	while(i!=InDataLength){
		if(InputFileData[i].Required && !InDataSet[i]){		//Required value not set
			nRequired++;
		}
		i++;	
	}

	if(nRequired!=0){

		cout << "Input file error: the following value";
		if(nRequired>1){
			cout << "s were ";
		}else{
			cout << " was ";
		}
		cout << "left unspecified." << endl;

		i=0;
		while(i!=InDataLength){
			if(InputFileData[i].Required && !InDataSet[i]){		//Required value not set
				cout << '\t' << InputFileData[i].DataName << endl;
			}
			i++;	
		}

		return 2;
	}

	return 0;		//Complete read-in
}
