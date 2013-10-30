#if _WIN32				//Windows headers
#include <Windows.h>
#else					//Linux headers
#include <sys/time.h>
#define _snprintf snprintf
#endif

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <float.h>
#include <climits>
#include <fstream>

#include <mpi.h>

using namespace std;

//// Input File Variables /////////////////////////////////////////////////////////////////////////////////////////////

char InputFileDefault[1024] = "Input.txt";			//Default input file if not specified in command line

char GeometryFile[1024] = "\0";

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

double BodyForceX = 1.0e-6;							//Body force on fluid
double BodyForceY = 0.0e-6;
double BodyForceZ = 0.0e-6;

int TimeStepsMax		= 100;						//End simulation after this many timesteps (0 = no limit)
int WriteInterval		= 50;						//Write simulation information to console after so many timesteps
int FileOutputInterval	= 1000;						//Output velocity field after so many timesteps

bool UseStopCondition = false;

int nTimeStepsCheck		= 0;						//Timesteps interval between checking stop condition. 0 corresponds to WriteInterval

double RelChangeAvVel  = 0.1;						//Stop when percentage change in quantity is less than value set here.
double RelChangeVarVel = 0.1;						//	AvVel = Average Velocity; VarVel = Variance of velocity field;

char OutputFilesFolder[1024] = "\0";

char OutputVelocityFile[1024] = "\0";

bool OutputVelocityFileBin = false;					//Is velocities file binary format
bool OutputVelocityFileSparse = false;				//Sparse velocity files
bool OutputVelocityFileVTKHeader = true;			//Output file has a VTK header (9 lines)

char OutputGeometryFile[1024] = "\0";				//Output lattice geometry

bool OutputGeometryFileVTKHeader = true;			//Output file has a VTK header (10 lines)
bool OutputGeometryFileBin = false;					//Output file in binary format

char DecompositionLogFile[1024] = "\0";
char DecompositionGeometryFile[1024] = "\0";

bool OutputDecompGeometry = false;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Other settings

int nDPVelOutput = 3;			//Number of decimal places for values in velocity field output file

//

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
		{ "FileOutputInterval", &FileOutputInterval, sizeof(FileOutputInterval), DataType_Int, true },

		{ "UseStopCondition", &UseStopCondition, sizeof(UseStopCondition), DataType_Bool, false },	
		{ "TimeStepsCheck", &nTimeStepsCheck, sizeof(nTimeStepsCheck), DataType_Int, true },
		{ "RelChangeAvVel", &RelChangeAvVel, sizeof(RelChangeAvVel), DataType_Double, false },
		{ "RelChangeVarVel", &RelChangeVarVel, sizeof(RelChangeVarVel), DataType_Double, false },
		
		{ "InputGeometryFile", GeometryFile, sizeof(GeometryFile), DataType_String, true },
			{ "GeometryFileBin", &GeometryFileBin, sizeof(GeometryFileBin), DataType_Bool, true },
			{ "GeometryFileVTKHeader", &GeometryFileVTKHeader, sizeof(GeometryFileVTKHeader), DataType_Bool, true },
			{ "NLinesGeometryVTKHeader", &NLinesGeometryVTKHeader, sizeof(NLinesGeometryVTKHeader), DataType_Int, true },

		{ "OutputFilesFolder", OutputFilesFolder, sizeof(OutputFilesFolder), DataType_String, true },
			
		{ "OutputGeometryFile", OutputGeometryFile, sizeof(OutputGeometryFile), DataType_String, true },
			{ "OutputGeometryFileBin", &OutputGeometryFileBin, sizeof(OutputGeometryFileBin), DataType_Bool, true },
			{ "OutputGeometryFileVTKHeader", &OutputGeometryFileVTKHeader, sizeof(OutputGeometryFileVTKHeader), DataType_Bool, true },
		{ "OutputVelocityFile", OutputVelocityFile, sizeof(OutputVelocityFile), DataType_String, true },
			{ "OutputVelocityFileBin", &OutputVelocityFileBin, sizeof(OutputVelocityFileBin), DataType_Bool, true },
			{ "OutputVelocityFileSparse", &OutputVelocityFileSparse, sizeof(OutputVelocityFileSparse), DataType_Bool, true },
			{ "OutputVelocityFileVTKHeader", &OutputVelocityFileVTKHeader, sizeof(OutputVelocityFileVTKHeader), DataType_Bool, true },
		{ "DecompositionLogFile", DecompositionLogFile, sizeof(DecompositionLogFile), DataType_String, true },
		{ "DecompositionGeometry", DecompositionGeometryFile, sizeof(DecompositionGeometryFile), DataType_String, false },
			{ "OutputDecompGeometry", &OutputDecompGeometry, sizeof(OutputDecompGeometry), DataType_Bool, false },	

	};

char* OutputPath(char* FileName, int Param){

	char* Folder = OutputFilesFolder;

	const int MaxPathLength = 1024;

	static char Path[MaxPathLength];

	int Folderl = (int)strlen(Folder);
	int Filel = (int)strlen(FileName);

	bool fs = true;
	int PathIndex = 0;

	for(int i=0; i<Folderl; i++){

		Path[PathIndex] = Folder[i];
		PathIndex++;

		if(PathIndex == MaxPathLength){ return 0; }

		if(Folder[i]=='\\'){
			fs = false;
		}

	}

	if(!((Folder[Folderl-1]=='\\'||Folder[Folderl-1]=='/') || (FileName[0]=='\\'||FileName[0]=='/'))){
		if(fs){
			Path[PathIndex] = '/';
		}else{
			Path[PathIndex] = '\\';
		}
		PathIndex++;
		if(PathIndex == MaxPathLength){return 0;}
	}

	for(int i=0; i<Filel; i++){

		if(FileName[i]=='%' && i<(Filel-1)){
		if(FileName[i+1]=='N' || FileName[i+1]=='n'){	//Print number of timesteps

				int n = _snprintf(&Path[PathIndex],MaxPathLength-PathIndex,"%d",Param);

				i+=1;
				PathIndex+=n;
				if(PathIndex >= MaxPathLength){return 0;}

				continue;

		}
		}

		Path[PathIndex] = FileName[i];
		PathIndex++;
		if(PathIndex >= MaxPathLength){return 0;}
	}

	Path[PathIndex] = '\0';

	return Path;

}

char* OutputPath(char* FileName){

	return OutputPath(FileName, 0);

}

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

////////////////////

//MPI Variables

struct CPUDomain{		//Defines block

	int ThreadID;		//Thread ID associated with this CPU Domain

	int x0;				//Indices of CPU domain voxel range
	int x1;				//Index of last voxel in x range
	int y0;
	int y1;
	int z0;
	int z1;

	int xn;				//Overlap layer x coordinates
	int xp;
	int yn;				//Overlap layer y coordinates
	int yp;
	int zn;				//Overlap layer z coordinates
	int zp;

	int DomainWidthX;	//Dimensions of voxel array
	int DomainWidthY;
	int DomainWidthZ;
	
};

CPUDomain ThreadDomain;			//Current MPI Thread's block

inline void DomainToLattice(double& x, double& y, double& z){
	x += (double)ThreadDomain.x0;
	y += (double)ThreadDomain.y0;
	z += (double)ThreadDomain.z0;
}

inline void LatticeToDomain(double& x, double& y, double& z){
	x -= (double)ThreadDomain.x0;
	y -= (double)ThreadDomain.y0;
	z -= (double)ThreadDomain.z0;
}

inline void DomainToLattice(int& x, int& y, int& z){
	x += ThreadDomain.x0;
	y += ThreadDomain.y0;
	z += ThreadDomain.z0;
}

inline void LatticeToDomain(int& x, int& y, int& z){
	x -= ThreadDomain.x0;
	y -= ThreadDomain.y0;
	z -= ThreadDomain.z0;
}

char* _DecompositionLattice;		//For decomposition, temporary geometry array

bool _DecompositionGetLattice(int x, int y, int z){

	int Index = NLatticeZ*(NLatticeY*x + y) + z;
	int rShift = Index & 0x00000007;	//Keep last 3 bits
	int ArrIndex = Index>>3;			// divide by 8	

	return (_DecompositionLattice[ArrIndex]>>rShift)&0x01;
}

void _DecompositionSetLattice(int x, int y, int z, bool Value){

	int Index = NLatticeZ*(NLatticeY*x + y) + z;
	int rShift = Index & 0x00000007;	//Keep last 3 bits
	int ArrIndex = Index>>3;			// divide by 8

	if(Value){
		char Op = ((0x01)<<rShift);
		_DecompositionLattice[ArrIndex] |= Op;
	}else{
		char Op = ~((0x01)<<rShift);
		_DecompositionLattice[ArrIndex] &= Op;
	}
}

struct LatticeBlock{
	int c0;					//Lower bound coordinate integer part
	double c0_;				//Lower bound coordinate fractional part
	int c1;					//Upper bound coordinate integer part
	double c1_;				//Upper bound coordinate fractional part

	bool LowerRound;		//Round lower border down or up
	bool UpperRound;		//Round upper border down or up

	int c0Rounded;			//Rounded lower border
	int c1Rounded;			//Rounded upper border

	double n_d;				//Exact number of slices in plane

	int NSquares;			//Number of processor blocks in plane

	LatticeBlock* SubBlocks;
	int SubBlockNum;
};

CPUDomain* Domains;				//Array of domain coordinates

////////////////////

//MPI Transfer Variables

enum Vectors{
	VectorsXp,
	VectorsXn,
	VectorsYp,
	VectorsYn,
	VectorsZp,
	VectorsZn,
	VectorsXpYp,
	VectorsXnYp,
	VectorsXpYn,
	VectorsXnYn,
	VectorsYpZp,
	VectorsYnZp,
	VectorsYpZn,
	VectorsYnZn,
	VectorsXpZp,
	VectorsXnZp,
	VectorsXpZn,
	VectorsXnZn
};

struct TransferVectors{
	int nVectors;
	int Vectors[5];		//Up to 5 vectors transferred
	int XSide;			// 1 for postive, -1 for negative else 0
	int YSide;
	int ZSide;

	int x0;				//Nodes on this side, globally
	int x1;
	int y0;
	int y1;
	int z0;
	int z1;

	int Domainx0;				//Nodes on this side, locally
	int Domainx1;
	int Domainy0;
	int Domainy1;
	int Domainz0;
	int Domainz1;				
};

TransferVectors TransferTypes[18];		//Transfer vectors to transfer on each side initialised in InitialiseMPITransfer

struct TransferRequest{		//Request data from another thread

	int Thread;		
	int RequestID;			//ID of request

	int x0;					//Range of nodes requested (global coordinates)
	int x1;
	int y0;
	int y1;
	int z0;
	int z1;

	int NonSolidCount;		//Number of non-solid elements in range

	int Vectors;			//Which set of vectors requested

	double* Data;			//For requests, allocate array for data

	MPI_Request MPIRequest;
	int RequestStatus;
};

TransferRequest* RecvRequests;	//Requests from this thread to other threads for data
int nRecvRequest;

TransferRequest* SendRequests;	//Requests from other threads to this thread for data
int nSendRequest;

long long* TransferSolid[18];		//Arrays of solid nodes in the extra transfer layer of the domain
int    nTransferNonSolid[18];		//Number of elements in each part of the layer

int nTransferNonSolids;				//Total number of non-solids in transfer layer

////////////////////



//Calculation variables


long long nVoxels;				//Total number of voxels (= NLatticeX * NLatticeY * NLatticeZ)
long long nNonSolid;			//Number of non-solid nodes in domain
long long nBulkNodes;			//Number of non-solid nodes in domain excluding the outer layer

long long nNonSolidLattice;		//Total number of non-solid nodes in lattice

long long* Solid;	//Array, nX*nY*nZ, -1 = Solid, 0 -> NVoxels-1 = index for voxel in data arrays

inline void SetSolid(int x, int y, int z, long long Value){
	Solid[ThreadDomain.DomainWidthX*(ThreadDomain.DomainWidthY*z + y) + x] = Value;
}

inline long long GetSolid(int x, int y, int z){
	return Solid[ThreadDomain.DomainWidthX*(ThreadDomain.DomainWidthY*z + y) + x];
}

long long GetNode(int x, int y, int z){				//May be in transfer layer

	//If domain spans entire lattice
	if(x==ThreadDomain.DomainWidthX && ThreadDomain.x0==0 && ThreadDomain.x1==NLatticeX-1){
		if(BoundaryConditionX == 1){
			return 0;					//Return solid node
		}else{
			x = 0;						//Apply loop boundary
		}
	}
	if(x==-1 && ThreadDomain.x0==0 && ThreadDomain.x1==NLatticeX-1){
		if(BoundaryConditionX == 1){
			return 0;					//Return solid node
		}else{
			x = NLatticeX-1;			//Apply loop boundary
		}
	}
	if(y==ThreadDomain.DomainWidthY && ThreadDomain.y0==0 && ThreadDomain.y1==NLatticeY-1){
		if(BoundaryConditionY == 1){
			return 0;					//Return solid node
		}else{
			y = 0;						//Apply loop boundary
		}
	}
	if(y==-1 && ThreadDomain.y0==0 && ThreadDomain.y1==NLatticeY-1){
		if(BoundaryConditionY == 1){
			return 0;					//Return solid node
		}else{
			y = NLatticeY-1;			//Apply loop boundary
		}
	}
	if(z==ThreadDomain.DomainWidthZ && ThreadDomain.z0==0 && ThreadDomain.z1==NLatticeZ-1){
		if(BoundaryConditionZ == 1){
			return 0;					//Return solid node
		}else{
			z = 0;						//Apply loop boundary
		}
	}
	if(z==-1 && ThreadDomain.z0==0 && ThreadDomain.z1==NLatticeZ-1){
		if(BoundaryConditionZ == 1){
			return 0;					//Return solid node
		}else{
			z = NLatticeZ-1;			//Apply loop boundary
		}
	}

	
	if(x>=0 && x<ThreadDomain.DomainWidthX && y>=0 && y<ThreadDomain.DomainWidthY && z>=0 && z<ThreadDomain.DomainWidthZ){

		return GetSolid(x,y,z);

	}


	int xSide = (x == -1) ? -1 : 0;
	int ySide = (y == -1) ? -1 : 0;
	int zSide = (z == -1) ? -1 : 0;

	xSide = (x == ThreadDomain.DomainWidthX) ? 1 : xSide;
	ySide = (y == ThreadDomain.DomainWidthY) ? 1 : ySide;
	zSide = (z == ThreadDomain.DomainWidthZ) ? 1 : zSide;

	for(int i=0; i<18; i++){

		if(TransferTypes[i].XSide == xSide && TransferTypes[i].YSide == ySide && TransferTypes[i].ZSide == zSide){

			int x_ = (xSide==0) ? x : 0;
			int y_ = (ySide==0) ? y : 0;
			int z_ = (zSide==0) ? z : 0;

			int nX = TransferTypes[i].x1 - TransferTypes[i].x0 + 1;
			int nY = TransferTypes[i].y1 - TransferTypes[i].y0 + 1;

			int Ind = nX*(nY*z_ + y_) + x_;

			return TransferSolid[i][Ind];

		}

	}

	return 0;
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

struct MacroVariables{
	double VAv[3];		//Average velocity
	double VAvMag;		//Magnitude of average velocity

	double VVar[3];		//Variance of the velocity field
	double VVarMag;		//Magnitude of velocity field variance

	double Perm[3];		//Permeability
	double PermMag;		//Magnitude of permeability vector
};

void ZeroMacroVarStruct(MacroVariables* s1){
	s1->VAv[0]  = 0;
	s1->VAv[1]  = 0;
	s1->VAv[2]  = 0;
	s1->VAvMag  = 0;
	s1->VVar[0] = 0;
	s1->VVar[1] = 0;
	s1->VVar[2] = 0;
	s1->VVarMag = 0;
	s1->Perm[0] = 0;
	s1->Perm[1] = 0;
	s1->Perm[2] = 0;
	s1->PermMag = 0;
}

void AddMacroVarStruct(MacroVariables* s1, MacroVariables* s2){	//Add s2 to s1

	s1->VAv[0] += s2->VAv[0];
	s1->VAv[1] += s2->VAv[1];
	s1->VAv[2] += s2->VAv[2];

	s1->VAvMag += s2->VAvMag;

	s1->VVar[0] += s2->VVar[0];
	s1->VVar[1] += s2->VVar[1];
	s1->VVar[2] += s2->VVar[2];

	s1->VVarMag += s2->VVarMag;

	s1->Perm[0] += s2->Perm[0];
	s1->Perm[1] += s2->Perm[1];
	s1->Perm[2] += s2->Perm[2];

	s1->PermMag += s2->PermMag;
}

void SubMacroVarStruct(MacroVariables* s1, MacroVariables* s2){	//Subtract s2 from s1

	s1->VAv[0] -= s2->VAv[0];
	s1->VAv[1] -= s2->VAv[1];
	s1->VAv[2] -= s2->VAv[2];

	s1->VAvMag -= s2->VAvMag;

	s1->VVar[0] -= s2->VVar[0];
	s1->VVar[1] -= s2->VVar[1];
	s1->VVar[2] -= s2->VVar[2];

	s1->VVarMag -= s2->VVarMag;

	s1->Perm[0] -= s2->Perm[0];
	s1->Perm[1] -= s2->Perm[1];
	s1->Perm[2] -= s2->Perm[2];

	s1->PermMag -= s2->PermMag;
}


//Round doubles
double _round(double d){
	double i;
	double f = modf(d, &i);

	return i + (double)((f >= 0.5) ? 1 : 0);
}

//Function declarations
	
	//Input File reading
int ReadInputFile(char* InputFileName);
bool ReadBool(char* Value);
bool CompareStr(char* DataName, char* Property);
int ReadProperty(ifstream& File, char* Property, char* Value);
void DistributeInputFileMPI(int* ret, int ThreadID, int ThreadNum);

	//Program
int Initialise(int ThreadID, int ThreadNum);				//Initialise arrays and read in from files
void Finalise();											//Delete arrays
void OutputVelocities(int ThreadID, int ThreadNum, int nTimeSteps);
void OutputGeometry(int ThreadID, int ThreadNum);			//Output geometry file
class SimulationTimer;										//Timer class

	//Decomposition
int LatticeDecomposition(int NProcessors, int ThreadID);
bool DecompositionReadGeometry(long long* NumberOfNonSolids);
bool DecompositionReadGeometryBin(long long* NumberOfNonSolids);
void DecompositionMinimiseError(LatticeBlock* Plane, int nBlocks, int* nCumulative, int NLattice, double NonSolidsPerCPU, ofstream& DecompositionLog);
void MPIDistributeDecomposition(int* ret, int ThreadID, int ThreadNum);
void OutputDecompositionGeometry(int ThreadNum);
void InitialiseMPITransfer(int ThreadID, int ThreadNum);
void TransferMPISend(int ThreadID, int ThreadNum, int* nWait);
void TransferMPIRecv(int ThreadID, int ThreadNum, int nWait);

	//Calculation
void Collision(long long i0, long long i1);
void StreamNodes();
void StreamTransfer();
void ComputeMacroVariables(MacroVariables* Vars, bool CalcVars, int ThreadID, int ThreadNum);
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

	//MPI

	int ThreadID, ThreadNum;

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ThreadNum);	//Get number of threads running
	MPI_Comm_rank(MPI_COMM_WORLD, &ThreadID);	//Get thread ID (index from 0 -> ThreadNum-1)


	int ret;	//Success flag


	//Input file

	if(ThreadID==0){

		char* InputFile = InputFileDefault;

		if(argc >= 2){		//Assume file specified as second command argument
			InputFile = argv[1];
		}else{
			cout << "No input file specified: defaulting to \"" << InputFileDefault << "\"" << endl;
		}

		ret = ReadInputFile(InputFile);

		if(ret == 1){
			cout << "The input file \"" << InputFile << "\" could not be opened" << endl;
		}else if(ret == 2){
			cout << "The input file did not specify all required values" << endl;
		}

	}

	DistributeInputFileMPI(&ret, ThreadID, ThreadNum);	//Send input variables to all threads

	if(ret!=0){
		goto End;
	}


	//Lattice decomposition

	if(ThreadID==0){

		cout << "Initialising lattice decomposition with " << ThreadNum << " threads" << endl;

		ret = LatticeDecomposition(ThreadNum, ThreadID);	//Decompose lattice

	}

	MPIDistributeDecomposition(&ret, ThreadID, ThreadNum);

	if(ret!=0){
		goto End;
	}


	//Initialise MPI transfer

	InitialiseMPITransfer(ThreadID, ThreadNum);


	//Initialise constants and neighbouring nodes

	ret = Initialise(ThreadID, ThreadNum);

	if(ret!=0){
		goto End;
	}


	//Output geometry

	OutputGeometry(ThreadID, ThreadNum);



	//Simulation

	if(ThreadID==0){
		cout << "Starting Calculation" << endl << endl;
	}
	
	SimulationTimer* Timer = new SimulationTimer();		//Track simulation time

	//

	TimeStepsMax = (TimeStepsMax < 1) ? INT_MAX : TimeStepsMax;		//User setting TimeStepsMax = 0 means no limit
	FileOutputInterval = (FileOutputInterval < 1) ? 0 : FileOutputInterval;
	WriteInterval = (WriteInterval < 1) ? 0 : WriteInterval;
	nTimeStepsCheck = (nTimeStepsCheck > 0) ? min( nTimeStepsCheck, WriteInterval ) : WriteInterval;		//Check stop condition
	
	UseStopCondition = (UseStopCondition && !(RelChangeAvVel<=0.0 && RelChangeVarVel<=0.0));

	//

	MacroVariables Vars;		
	MacroVariables VarsLastOut;
	MacroVariables VarsLastCheck;

	int nWait;
	bool CalcVars = true;

	for(int n=1; n<=TimeStepsMax; n++){

		bool Write = (WriteInterval!=0) && (n%WriteInterval==0);	//Write to console
		bool Check = UseStopCondition && (n%nTimeStepsCheck==0);	//Check stop condition
		CalcVars = (Write || Check);


		//Calculation steps

		Collision(nBulkNodes+1, nNonSolid);		//Collide outer layer

		TransferMPISend(ThreadID, ThreadNum, &nWait);

		Collision(1, nBulkNodes);				//Collide rest of domain
			
		StreamNodes();

		TransferMPIRecv(ThreadID, ThreadNum, nWait);

		StreamTransfer();

		ComputeMacroVariables(&Vars, CalcVars, ThreadID, ThreadNum);


		//Output to console

		if(Write && ThreadID==0){

			printf("After %i time-steps:\n", n);

			if(n == WriteInterval){		//First output

				printf("Average Velocity   = (%.5e, %.5e, %.5e)\n", Vars.VAv[0], Vars.VAv[1], Vars.VAv[2]);
				printf("Velocity magnitude = %.5e\n", Vars.VAvMag);
				printf("Variance magnitude = %.5e\n", Vars.VVarMag);

			}else{

				MacroVariables VarsDiff;
				memcpy(&VarsDiff, &Vars, sizeof(MacroVariables));

				SubMacroVarStruct(&VarsDiff, &VarsLastOut);		//Calculate difference				

				//Percentage change in average velocity and variance since last write interval

				double PercentDiffAvV  = 100 * (VarsDiff.VAvMag / VarsLastOut.VAvMag);
				double PercentDiffVVar = 100 * (VarsDiff.VVarMag / VarsLastOut.VVarMag);

				char Sign1 = (PercentDiffAvV>0)  ? '+' : '-';
				char Sign2 = (PercentDiffVVar>0) ? '+' : '-';

				printf("Average Velocity   = (%.5e, %.5e, %.5e)\n", Vars.VAv[0], Vars.VAv[1], Vars.VAv[2]);
				printf("Velocity magnitude = %.5e (%c%.2f%%)\n", Vars.VAvMag, Sign1, fabs(PercentDiffAvV));
				printf("Variance magnitude = %.5e (%c%.2f%%)\n", Vars.VVarMag, Sign2, fabs(PercentDiffVVar));

				double InvNWrite = 1 / (double)WriteInterval;

				printf("Per time-step percentage changes: ");
				printf("Velocity (%c%.3f%%); ", Sign1, fabs(PercentDiffAvV*InvNWrite));
				printf("Variance (%c%.3f%%)\n", Sign2, fabs(PercentDiffVVar*InvNWrite));
				
			}

			memcpy(&VarsLastOut, &Vars, sizeof(MacroVariables));

			printf("Elapsed time: %.1f seconds\n\n", Timer->GetSimulationTime());

			fflush(stdout);
		}


		//Check stop condition

		if(Check){
			if(n == nTimeStepsCheck){	//First check
				memcpy(&VarsLastCheck, &Vars, sizeof(MacroVariables));
			}else{
				MacroVariables VarsDiff;
				memcpy(&VarsDiff, &Vars, sizeof(MacroVariables));

				SubMacroVarStruct(&VarsDiff, &VarsLastCheck);		//Calculate difference

				double InvNCheck = 1 / (double)nTimeStepsCheck;

				//Percentage change per timestep

				double PercentDiffAvV  = fabs( 100 * (VarsDiff.VAvMag / VarsLastCheck.VAvMag) * InvNCheck );		
				double PercentDiffVVar = fabs( 100 * (VarsDiff.VVarMag / VarsLastCheck.VVarMag) * InvNCheck );

				bool ConditionAvV = (RelChangeAvVel <=0 || PercentDiffAvV < RelChangeAvVel );
				bool ConditionVar = (RelChangeVarVel<=0 || PercentDiffVVar < RelChangeVarVel);

				if(ConditionAvV && ConditionVar){	//Stop condition met

					if(ThreadID==0){
						printf("Stop conditions met after %i time-steps, ending calculation\n\n", n);
						fflush(stdout);
					}

					TimeStepsMax = n;
				}

				memcpy(&VarsLastCheck, &Vars, sizeof(MacroVariables));
			}
		}


		//Output velocity file

		if(FileOutputInterval!=0 && n%FileOutputInterval==0){	

			if(ThreadID==0){
				printf("Outputting velocity field file after %i time-steps\n\n", n);
				fflush(stdout);
			}

			OutputVelocities(ThreadID, ThreadNum, n);

		}

	}

	
	//Final Outputs

	if(FileOutputInterval==0 || TimeStepsMax%FileOutputInterval!=0){

		if(ThreadID==0){
			printf("Outputting final velocity field file after %i time-steps\n", TimeStepsMax);
			fflush(stdout);
		}

		OutputVelocities(ThreadID, ThreadNum, TimeStepsMax);

	}


	//

	delete Timer;

End:
	Finalise();				//Delete variables
	MPI_Finalize();			//End of multithreading

	return 0;
}


void ComputeMacroVariables(MacroVariables* Vars, bool CalcVars, int ThreadID, int ThreadNum){

	MacroVariables VarsSum;
	ZeroMacroVarStruct(&VarsSum);

	const double dt = 1.0;

	double dtF_2[3];
	
	dtF_2[0] = 0.5*dt*BodyForceX;
	dtF_2[1] = 0.5*dt*BodyForceY;
	dtF_2[2] = 0.5*dt*BodyForceZ;

	for(long long i=1; i<=nNonSolid; i++){

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

		if(CalcVars){
		
			VarsSum.VAv[0] += u[0];			//Add velocities for average
			VarsSum.VAv[1] += u[1];
			VarsSum.VAv[2] += u[2];
		
			VarsSum.VVar[0] += (u[0]*u[0]);	//Add velocities square for variance
			VarsSum.VVar[1] += (u[1]*u[1]);
			VarsSum.VVar[2] += (u[2]*u[2]);

		}

	}

	if(!CalcVars){		//No need to calculate simulation variables this time
		return;
	}

	//Obtain sums from all threads

	MacroVariables VarsTot;
	ZeroMacroVarStruct(&VarsTot);
	AddMacroVarStruct(&VarsTot, &VarsSum);

	MacroVariables VarsRecv;

	for(int i=0; i<ThreadNum; i++){
		if(i==ThreadID){
			MPI_Bcast(&VarsSum, sizeof(MacroVariables), MPI_BYTE, ThreadID, MPI_COMM_WORLD);		//Send variables
		}else{
			MPI_Bcast(&VarsRecv, sizeof(MacroVariables), MPI_BYTE, i, MPI_COMM_WORLD);		//Receive
			
			AddMacroVarStruct(&VarsTot, &VarsRecv);
		}
	}

	double InvNonSolidLattice = 1 / (double)nNonSolidLattice;

	Vars->VAv[0] = VarsTot.VAv[0] * InvNonSolidLattice;			//Average velocity
	Vars->VAv[1] = VarsTot.VAv[1] * InvNonSolidLattice;
	Vars->VAv[2] = VarsTot.VAv[2] * InvNonSolidLattice;
	
	Vars->VVar[0] = (VarsTot.VVar[0] * InvNonSolidLattice) - (Vars->VAv[0] * Vars->VAv[0]);		//Variance of velocity field
	Vars->VVar[1] = (VarsTot.VVar[1] * InvNonSolidLattice) - (Vars->VAv[1] * Vars->VAv[1]);
	Vars->VVar[2] = (VarsTot.VVar[2] * InvNonSolidLattice) - (Vars->VAv[2] * Vars->VAv[2]);

	Vars->Perm[0] = Vars->VAv[0] * Viscosity / BodyForceX;
	Vars->Perm[1] = Vars->VAv[1] * Viscosity / BodyForceY;
	Vars->Perm[2] = Vars->VAv[2] * Viscosity / BodyForceZ;

	//Calculate magnitudes of average, variance and permeability

	Vars->VAvMag = sqrt(Vars->VAv[0]*Vars->VAv[0] + Vars->VAv[1]*Vars->VAv[1] + Vars->VAv[2]*Vars->VAv[2]);
	Vars->VVarMag = sqrt(Vars->VVar[0]*Vars->VVar[0] + Vars->VVar[1]*Vars->VVar[1] + Vars->VVar[2]*Vars->VVar[2]);
	Vars->PermMag = sqrt(Vars->Perm[0]*Vars->Perm[0] + Vars->Perm[1]*Vars->Perm[1] + Vars->Perm[2]*Vars->Perm[2]);

}

void StreamNodes(){

	for(long long i=1; i<=nNonSolid; i++){

		double* f = Nodes[i].f;				//Node distribution
		double* F = Nodes[i].fCalc;			//Temporary node distribution

		F[0] = f[0];						//Centre node

		for (int mi=1; mi<19; mi++){

			Node* NeighbourNode = Nodes[i].Neighbour[mi];

			if(NeighbourNode == 0){									//Bounce-back condition	
				F[ReflMap[mi]] = f[mi];
			}else{													//Stream to neighbour node
				NeighbourNode->fCalc[mi] = f[mi];
			}

		}

	}

}

void StreamTransfer(){

	//Stream from transfer layer

	long long NodesIndex = nNonSolid + 1;		//Transfer layer nodes start index

	for(int i=0; i<18; i++){

		int nVectors = TransferTypes[i].nVectors;

		int Vectors[5];

		for(int vi=0; vi<nVectors; vi++){
			Vectors[vi] = TransferTypes[i].Vectors[vi];
		}

		//Stream

		for(int c=0; c<nTransferNonSolid[i]; c++){

			Node* N = &Nodes[NodesIndex];

			for(int vi=0; vi<nVectors; vi++){

				Node* Neighbour = N->Neighbour[Vectors[vi]];

				if(Neighbour!=0){
					Neighbour->fCalc[Vectors[vi]] = N->f[Vectors[vi]];
				}

			}

			NodesIndex++;
		}

	}

}

void Collision(long long i0, long long i1){		//Do collision step on all nodes

	for(long long i=i0; i<=i1; i++){

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

		for(int mi=0; mi<19; mi++){
			f[mi] = MInvC[mi] + f[mi];
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

void TransferMPISend(int ThreadID, int ThreadNum, int* nWait){		//Send vectors and set recv

	//Send data to other threads

	for(int i=0; i<nSendRequest; i++){

		TransferRequest* SRq = &SendRequests[i];

		if(SRq->NonSolidCount==0){
			continue;
		}

		int x0 = SRq->x0;
		int x1 = SRq->x1;
		int y0 = SRq->y0;
		int y1 = SRq->y1;
		int z0 = SRq->z0;
		int z1 = SRq->z1;

		LatticeToDomain(x0,y0,z0);
		LatticeToDomain(x1,y1,z1);

		int nVectors = TransferTypes[SRq->Vectors].nVectors;

		int Vectors[5];

		for(int vi=0; vi<nVectors; vi++){
			Vectors[vi] = TransferTypes[SRq->Vectors].Vectors[vi];
		}

		int c=0;
		for(int z=z0; z<=z1; z++){
		for(int y=y0; y<=y1; y++){
		for(int x=x0; x<=x1; x++){

			long long si = GetSolid(x,y,z);

			if(si==0){
				continue;
			}

			for(int vi=0; vi<nVectors; vi++){
				SRq->Data[c] = Nodes[si].f[Vectors[vi]];
				
				c++;
			}

		}
		}
		}

		

		int DLength = SRq->NonSolidCount*nVectors;

		MPI_Isend(SRq->Data, DLength, MPI_DOUBLE, SRq->Thread, SRq->RequestID, MPI_COMM_WORLD, &(SRq->MPIRequest));

	}


	//Receive data from other threads

	int _nWait = 0;

	for(int i=0; i<nRecvRequest; i++){

		TransferRequest* Rq = &RecvRequests[i];

		if(Rq->NonSolidCount==0){
			continue;
		}

		int DLength = Rq->NonSolidCount * TransferTypes[Rq->Vectors].nVectors;

		MPI_Irecv(Rq->Data, DLength, MPI_DOUBLE, Rq->Thread, Rq->RequestID, MPI_COMM_WORLD, &(Rq->MPIRequest));

		Rq->RequestStatus = 1;		//Awaiting data

		_nWait++;

	}

	*nWait = _nWait;

}

void TransferMPIRecv(int ThreadID, int ThreadNum, int nWait){

	MPI_Status Stat;

	//Wait for all data

	while(nWait>0){					

		for(int i=0; i<nRecvRequest; i++){

			TransferRequest* Rq = &RecvRequests[i];

			if(Rq->NonSolidCount==0 || Rq->RequestStatus==0){	//Data already received or no reception necessary
				continue;
			}

			int Flag;
			MPI_Test(&(Rq->MPIRequest), &Flag, &Stat);			//Test for receipt

			if(!Flag){											//Not received yet
				continue;
			}

			//Put data into transfer layer

			int x0 = Rq->x0;
			int x1 = Rq->x1;
			int y0 = Rq->y0;
			int y1 = Rq->y1;
			int z0 = Rq->z0;
			int z1 = Rq->z1;

			TransferVectors* Transfer = &TransferTypes[Rq->Vectors];

			int nVectors = Transfer->nVectors;

			int Vectors[5];

			for(int vi=0; vi<nVectors; vi++){
				Vectors[vi] = Transfer->Vectors[vi];
			}

			int c=0;
			for(int z=z0; z<=z1; z++){
			for(int y=y0; y<=y1; y++){
			for(int x=x0; x<=x1; x++){

				int x_ = x - Transfer->x0;
				int y_ = y - Transfer->y0;
				int z_ = z - Transfer->z0;

				int nX = Transfer->x1 - Transfer->x0 + 1;
				int nY = Transfer->y1 - Transfer->y0 + 1;

				int Ind = nX*(nY*z_ + y_) + x_;

				long long si = TransferSolid[Rq->Vectors][Ind];

				if(si==0){
					continue;
				}

				for(int vi=0; vi<nVectors; vi++){
					Nodes[si].f[Vectors[vi]] = Rq->Data[c];
					c++;
				}

			}
			}
			}

			//

			Rq->RequestStatus = 0;
			nWait--;

		}

	}

}


int Initialise(int ThreadID, int ThreadNum){

	nVoxels = (long long)NLatticeX * (long long)NLatticeY * (long long)NLatticeZ;

	long long nNodes = nNonSolid + nTransferNonSolids;

	Nodes = new Node[nNodes+1];					//Create node data. Node[0] = Solid
	
	if(Nodes==0){
		cout << "Unable to allocate memory (1)" << endl;
		return 1;
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

	//Initialise transfer layer nodes

	long long Ind = nNonSolid + 1;				//Start index for transfer layer nodes
	for(int i=0; i<18; i++){

		int x0 = TransferTypes[i].x0;	//Range of transfer layer section
		int x1 = TransferTypes[i].x1;
		int y0 = TransferTypes[i].y0;
		int y1 = TransferTypes[i].y1;
		int z0 = TransferTypes[i].z0;
		int z1 = TransferTypes[i].z1;

		int c=0;
		for(int z=z0; z<=z1; z++){
		for(int y=y0; y<=y1; y++){
		for(int x=x0; x<=x1; x++){

			if(TransferSolid[i][c] != 0){		//Non-solid
				TransferSolid[i][c] = Ind;
				Ind++;
			}

			c++;
		}
		}
		}

	}

	//Initialise lattice

	long long NodeIndex = 1;

	for(int z=1; z<ThreadDomain.DomainWidthZ-1; z++){		//Bulk nodes
	for(int y=1; y<ThreadDomain.DomainWidthY-1; y++){
	for(int x=1; x<ThreadDomain.DomainWidthX-1; x++){

		long long SInd = GetSolid(x,y,z);

		if(SInd == 0){		//Solid
			continue;
		}

		SetSolid(x,y,z,NodeIndex);					//Index in nodes array
		NodeIndex++;
		
	}
	}
	}

	nBulkNodes = NodeIndex - 1;

	for(int z=0; z<ThreadDomain.DomainWidthZ; z++){		//Outer layer nodes
	for(int y=0; y<ThreadDomain.DomainWidthY; y++){
	for(int x=0; x<ThreadDomain.DomainWidthX; x++){

		bool CondX = (x==0 || x==ThreadDomain.DomainWidthX-1);
		bool CondY = (y==0 || y==ThreadDomain.DomainWidthY-1);
		bool CondZ = (z==0 || z==ThreadDomain.DomainWidthZ-1);

		if(!CondX && !CondY && !CondZ){
			continue;
		}

		long long SInd = GetSolid(x,y,z);

		if(SInd == 0){		//Solid
			continue;
		}
		
		SetSolid(x,y,z,NodeIndex);					//Index in nodes array
		NodeIndex++;

	}
	}
	}


	for(int z=0; z<ThreadDomain.DomainWidthZ; z++){
	for(int y=0; y<ThreadDomain.DomainWidthY; y++){
	for(int x=0; x<ThreadDomain.DomainWidthX; x++){

		long long i = GetSolid(x,y,z);

		if(i == 0){		//Solid
			continue;
		}

		Nodes[i].Density = 1.0;
		
		Nodes[i].Velocity[0] = InitialVx;
		Nodes[i].Velocity[1] = InitialVy;
		Nodes[i].Velocity[2] = InitialVz;

		for(int i2=0; i2!=19; i2++){		//Calculate distribution function from initial velocity
			Nodes[i].f[i2] = fEquilibrium(Nodes[i].Density, i2, Nodes[i].Velocity);
		}

		//Find neighbouring nodes
			
		int xp = x+1;
		int xn = x-1;
		int yp = y+1;
		int yn = y-1;
		int zp = z+1;
		int zn = z-1;

		Nodes[i].Neighbour[1]  = &Nodes[GetNode(xp,y,z)];	//+X
		Nodes[i].Neighbour[2]  = &Nodes[GetNode(xn,y,z)];	//-X
		Nodes[i].Neighbour[3]  = &Nodes[GetNode(x,yp,z)];	//+Y
		Nodes[i].Neighbour[4]  = &Nodes[GetNode(x,yn,z)];	//-Y
		Nodes[i].Neighbour[5]  = &Nodes[GetNode(x,y,zp)];	//+Z
		Nodes[i].Neighbour[6]  = &Nodes[GetNode(x,y,zn)];	//-Z

		Nodes[i].Neighbour[7]  = &Nodes[GetNode(xp,yp,z)];	//+X +Y
		Nodes[i].Neighbour[8]  = &Nodes[GetNode(xn,yp,z)];	//-X +Y
		Nodes[i].Neighbour[9]  = &Nodes[GetNode(xp,yn,z)];	//+X -Y
		Nodes[i].Neighbour[10] = &Nodes[GetNode(xn,yn,z)];	//-X -Y

		Nodes[i].Neighbour[11] = &Nodes[GetNode(x,yp,zp)];	//   +Y +Z
		Nodes[i].Neighbour[12] = &Nodes[GetNode(x,yn,zp)];	//   -Y +Z
		Nodes[i].Neighbour[13] = &Nodes[GetNode(x,yp,zn)];	//   +Y -Z
		Nodes[i].Neighbour[14] = &Nodes[GetNode(x,yn,zn)];	//   -Y -Z

		Nodes[i].Neighbour[15] = &Nodes[GetNode(xp,y,zp)];	//+X    +Z
		Nodes[i].Neighbour[16] = &Nodes[GetNode(xn,y,zp)];	//-X    +Z
		Nodes[i].Neighbour[17] = &Nodes[GetNode(xp,y,zn)];	//+X    -Z
		Nodes[i].Neighbour[18] = &Nodes[GetNode(xn,y,zn)];	//-X    -Z
		
			//Finally, set any &Nodes[0] pointers to 0 to indicate solid
		
		for(int i2=1; i2!=19; i2++){
			if(Nodes[i].Neighbour[i2] == &Nodes[0]){
				Nodes[i].Neighbour[i2] = 0;
			}
		}

	}
	}
	}


	//Initialise transfer layer neigbours

	for(int i=0; i<18; i++){

		int x0 = TransferTypes[i].Domainx0;	//Range of transfer layer section
		int x1 = TransferTypes[i].Domainx1;
		int y0 = TransferTypes[i].Domainy0;
		int y1 = TransferTypes[i].Domainy1;
		int z0 = TransferTypes[i].Domainz0;
		int z1 = TransferTypes[i].Domainz1;

		int c=0;
		for(int z=z0; z<=z1; z++){
		for(int y=y0; y<=y1; y++){
		for(int x=x0; x<=x1; x++){

			if(TransferSolid[i][c] == 0){		//Solid
				c++;
				continue;
			}

			long long sInd = TransferSolid[i][c];

			for(int vi=0; vi<TransferTypes[i].nVectors; vi++){

				int x_ = x + BasisVector[TransferTypes[i].Vectors[vi]][0];
				int y_ = y + BasisVector[TransferTypes[i].Vectors[vi]][1];
				int z_ = z + BasisVector[TransferTypes[i].Vectors[vi]][2];

				if(x_>=0&&x_<ThreadDomain.DomainWidthX && y_>=0&&y_<ThreadDomain.DomainWidthY && z_>=0&&z_<ThreadDomain.DomainWidthZ){

					long long Ind = GetSolid(x_,y_,z_);

					if(Ind!=0){
						Nodes[sInd].Neighbour[TransferTypes[i].Vectors[vi]] = &Nodes[Ind];
					}else{
						Nodes[sInd].Neighbour[TransferTypes[i].Vectors[vi]] = 0;
					}

				}else{

						Nodes[sInd].Neighbour[TransferTypes[i].Vectors[vi]] = 0;

				}


			}

			c++;
		}
		}
		}

	}

	return 0;
}



void Finalise(){

	for(int i=0; i<nSendRequest; i++){
		delete[] SendRequests[i].Data;
	}

	for(int i=0; i<nRecvRequest; i++){
		delete[] RecvRequests[i].Data;
	}

	delete[] SendRequests;
	delete[] RecvRequests;

	for(int i=0; i<18; i++){
		delete[] TransferSolid[i];
	}

	delete[] Nodes;
	delete[] Solid;
	
	delete[] Domains;

}

void OutputDomainGeometry(int ThreadID, int ThreadNum){

	char Str[256];
	_snprintf(Str, 256, "D:/LB/Geo%d.vtk", ThreadID);

	FILE* OutFile = fopen(Str, "wb");

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
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DIMENSIONS %i %i %i\n", ThreadDomain.DomainWidthX, ThreadDomain.DomainWidthY, ThreadDomain.DomainWidthZ);
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ORIGIN 0 0 0\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SPACING 1 1 1\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "POINT_DATA %i\n", (ThreadDomain.DomainWidthX*ThreadDomain.DomainWidthY*ThreadDomain.DomainWidthZ));
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SCALARS sample_scalars float\n");
		Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "LOOKUP_TABLE default\n");

		fwrite(Header, sizeof(char), Ind, OutFile);

	}

	char c0[2] = {'0', ' '};
	char c1[2] = {'1', ' '};

	for(int z=0; z<ThreadDomain.DomainWidthZ; z++){
	for(int y=0; y<ThreadDomain.DomainWidthY; y++){
	for(int x=0; x<ThreadDomain.DomainWidthX; x++){

		long long Index = GetSolid(x, y, z);

		if(Index==0){	//Solid
			fwrite(c1, sizeof(char), 2, OutFile);
		}else{
			fwrite(c0, sizeof(char), 2, OutFile);
		}

	}
	}
	}

	fclose(OutFile);

}

void OutputDomainVelocity(int ThreadID, int ThreadNum){

	char Str[256];
	_snprintf(Str, 256, "D:/LB/Vel%d.vtk", ThreadID);

	FILE* OutFile = fopen(Str, "wb");

	if(!OutFile){
		return;
	}

	char Header[2048];
	int Ind = 0;

	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "# vtk DataFile Version 2.0\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "Geometry Output\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ASCII\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DATASET STRUCTURED_POINTS\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DIMENSIONS %i %i %i\n", ThreadDomain.DomainWidthX, ThreadDomain.DomainWidthY, ThreadDomain.DomainWidthZ);
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ORIGIN 0 0 0\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SPACING 1 1 1\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "POINT_DATA %i\n", (ThreadDomain.DomainWidthX*ThreadDomain.DomainWidthY*ThreadDomain.DomainWidthZ));
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "VECTORS sample_scalars float\n");

	fwrite(Header, sizeof(char), Ind, OutFile);


	for(int z=0; z<ThreadDomain.DomainWidthZ; z++){
	for(int y=0; y<ThreadDomain.DomainWidthY; y++){
	for(int x=0; x<ThreadDomain.DomainWidthX; x++){

		long long Index = GetSolid(x, y, z);

		int l;

		if(Index==0){	//Solid
			l = _snprintf(Str, sizeof(Str), "0 0 0\n");
		}else{

			double V[3];

			V[0] = Nodes[Index].Velocity[0];
			V[1] = Nodes[Index].Velocity[1];
			V[2] = Nodes[Index].Velocity[2];

			l = _snprintf(Str, sizeof(Str), "%e %e %e\n", V[0], V[1], V[2]);
		}
				
		fwrite(Str, sizeof(char), l, OutFile);

	}
	}
	}

	fclose(OutFile);

}

void OutputVelocities(int ThreadID, int ThreadNum, int nTimeSteps){

	MPI_Status Stat;

	int Coords[4];	//[0] = x0, [1] = x1, [2] = y, [3] = z

	if(ThreadID==0){
		
		FILE* OutFile = fopen(OutputPath(OutputVelocityFile, nTimeSteps), "wb");

		if(!OutFile){
			Coords[0] = -1;
			for(int i=1; i<ThreadNum; i++){
				MPI_Send(&Coords, 4, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
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
			Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "POINT_DATA %lli\n", nVoxels);
			Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "VECTORS sample_scalars float\n");

			fwrite(Header, sizeof(char), Ind, OutFile);

		}

		for(int z=0; z<NLatticeZ; z++){
		for(int y=0; y<NLatticeY; y++){
		for(int x=0; x<NLatticeX;    ){

			int x1;
			int tId;

			for(int i=0; i<ThreadNum; i++){
				if(x>=Domains[i].x0&&x<=Domains[i].x1 && y>=Domains[i].y0&&y<=Domains[i].y1 && z>=Domains[i].z0&&z<=Domains[i].z1){
					x1 = Domains[i].x1;
					tId = i;
					break;
				}
			}

			Coords[0] = x;
			Coords[1] = x1;
			Coords[2] = y;
			Coords[3] = z;

			int nx = Coords[1] - Coords[0] + 1;

			double* V = new double[nx*3];

			if(tId!=ThreadID){

				MPI_Send(&Coords, 4, MPI_INT, tId, 0, MPI_COMM_WORLD);

				MPI_Recv(V, nx*3, MPI_DOUBLE, tId, 0, MPI_COMM_WORLD, &Stat);

			}else{

				LatticeToDomain(Coords[0], Coords[2], Coords[3]);

				for(int i=0; i<nx; i++){

					long long si = GetSolid(Coords[0]+i, Coords[2], Coords[3]);

					if(si!=0){
						Node* N = &Nodes[si];

						V[i*3    ] = N->Velocity[0];
						V[i*3 + 1] = N->Velocity[1];
						V[i*3 + 2] = N->Velocity[2];
					}else{
						V[i*3    ] = 0;
						V[i*3 + 1] = 0;
						V[i*3 + 2] = 0;
					}

				}

			}

			char Str[128];
			int l;

			for(int i=0; i<nx; i++){

				if(V[i*3] == 0 && V[i*3 + 1] == 0 && V[i*3 + 2] == 0){
					l = _snprintf(Str, sizeof(Str), "0 0 0\n");
				}else{
					l = _snprintf(Str, sizeof(Str), "%.*e %.*e %.*e\n", nDPVelOutput, V[i*3], nDPVelOutput, V[i*3 + 1], nDPVelOutput, V[i*3 + 2]);
				}
				
				fwrite(Str, sizeof(char), l, OutFile);

			}

			delete[] V;

			x = x1 + 1;

		}
		}
		}

		Coords[0] = -1;
		for(int i=1; i<ThreadNum; i++){
			MPI_Send(&Coords, 4, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		fclose(OutFile);

	}else{

		while(true){

			int Coords[4];	//[0] = x0, [1] = x1, [2] = y, [3] = z
			MPI_Recv(&Coords, 4, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);

			if(Coords[0] == -1){
				return;
			}

			int nx = Coords[1] - Coords[0] + 1;

			LatticeToDomain(Coords[0], Coords[2], Coords[3]);

			double* V = new double[nx*3];

			for(int i=0; i<nx; i++){

				long long si = GetSolid(Coords[0]+i, Coords[2], Coords[3]);

				if(si!=0){
					Node* N = &Nodes[si];

					V[i*3  ] = N->Velocity[0];
					V[i*3+1] = N->Velocity[1];
					V[i*3+2] = N->Velocity[2];

				}else{
					V[i*3  ] = 0;
					V[i*3+1] = 0;
					V[i*3+2] = 0;
				}

			}

			MPI_Send(V, nx*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

			delete[] V;

		}
	}

}


void OutputGeometry(int ThreadID, int ThreadNum){

	MPI_Status Stat;

	int Coords[4];	//[0] = x0, [1] = x1, [2] = y, [3] = z

	if(ThreadID==0){
		
		FILE* OutFile = fopen(OutputPath(OutputGeometryFile), "wb");

		if(!OutFile){
			Coords[0] = -1;
			for(int i=1; i<ThreadNum; i++){
				MPI_Send(&Coords, 4, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
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
			Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "POINT_DATA %lli\n", nVoxels);
			Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SCALARS sample_scalars float\n");
			Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "LOOKUP_TABLE default\n");

			fwrite(Header, sizeof(char), Ind, OutFile);

		}

		for(int z=0; z<NLatticeZ; z++){
		for(int y=0; y<NLatticeY; y++){
		for(int x=0; x<NLatticeX;    ){

			int x1;
			int tId;

			for(int i=0; i<ThreadNum; i++){
				if(x>=Domains[i].x0&&x<=Domains[i].x1 && y>=Domains[i].y0&&y<=Domains[i].y1 && z>=Domains[i].z0&&z<=Domains[i].z1){
					x1 = Domains[i].x1;
					tId = i;
					break;
				}
			}

			Coords[0] = x;
			Coords[1] = x1;
			Coords[2] = y;
			Coords[3] = z;

			int nx = Coords[1] - Coords[0] + 1;

			bool* Solid = new bool[nx];

			if(tId!=ThreadID){

				MPI_Send(&Coords, 4, MPI_INT, tId, 0, MPI_COMM_WORLD);

				MPI_Recv(Solid, nx, MPI_BYTE, tId, 0, MPI_COMM_WORLD, &Stat);

			}else{

				LatticeToDomain(Coords[0], Coords[2], Coords[3]);

				for(int i=0; i<nx; i++){
					Solid[i] = (GetSolid(Coords[0]+i, Coords[2], Coords[3])==0);
				}

			}
			
			char c0[2] = {'0', ' '};
			char c1[2] = {'1', ' '};

			for(int i=0; i<nx; i++){

				if(Solid[i]){	//Solid
					fwrite(c1, sizeof(char), 2, OutFile);
				}else{
					fwrite(c0, sizeof(char), 2, OutFile);
				}

			}

			delete[] Solid;

			x = x1 + 1;

		}
		}
		}

		Coords[0] = -1;
		for(int i=1; i<ThreadNum; i++){
			MPI_Send(&Coords, 4, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		fclose(OutFile);

	}else{

		while(true){

			int Coords[4];	//[0] = x0, [1] = x1, [2] = y, [3] = z
			MPI_Recv(&Coords, 4, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);

			if(Coords[0] == -1){
				return;
			}

			int nx = Coords[1] - Coords[0] + 1;

			LatticeToDomain(Coords[0], Coords[2], Coords[3]);

			bool* Solid = new bool[nx];

			for(int i=0; i<nx; i++){
				Solid[i] = (GetSolid(Coords[0]+i, Coords[2], Coords[3])==0);
			}

			MPI_Send(Solid, nx, MPI_BYTE, 0, 0, MPI_COMM_WORLD);

			delete[] Solid;

		}
	}

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

void DistributeInputFileMPI(int* ret, int ThreadID, int ThreadNum){

	MPI_Bcast(ret, 1, MPI_INT, 0, MPI_COMM_WORLD);		//Share input file success flag

	const int InDataLength = sizeof(InputFileData)/sizeof(InputValueStruct);

	for(int i=0; i<InDataLength; i++){
		MPI_Bcast(InputFileData[i].Ptr, InputFileData[i].PtrLength, MPI_BYTE, 0, MPI_COMM_WORLD);
	}

}


void InitialiseMPITransfer(int ThreadID, int ThreadNum){

	MPI_Status Stat;


	//Initialise transfer information

	memset(TransferTypes, 0, sizeof(TransferTypes));

	TransferTypes[0].XSide = 1;			//+X
	TransferTypes[1].XSide = -1;		//-X
	TransferTypes[2].YSide = 1;			//   +Y
	TransferTypes[3].YSide = -1;		//   -Y
	TransferTypes[4].ZSide = 1;			//      +Z
	TransferTypes[5].ZSide = -1;		//      -Z
	TransferTypes[6].XSide = 1;			//+X +Y
	TransferTypes[6].YSide = 1;
	TransferTypes[7].XSide = -1;		//-X +Y
	TransferTypes[7].YSide = 1;
	TransferTypes[8].XSide = 1;			//+X -Y
	TransferTypes[8].YSide = -1;
	TransferTypes[9].XSide = -1;		//-X -Y
	TransferTypes[9].YSide = -1;
	TransferTypes[10].YSide = 1;		//   +Y +Z
	TransferTypes[10].ZSide = 1;
	TransferTypes[11].YSide = -1;		//   -Y +Z
	TransferTypes[11].ZSide = 1;
	TransferTypes[12].YSide = 1;		//   +Y -Z
	TransferTypes[12].ZSide = -1;
	TransferTypes[13].YSide = -1;		//   -Y -Z
	TransferTypes[13].ZSide = -1;
	TransferTypes[14].XSide = 1;		//+X    +Z
	TransferTypes[14].ZSide = 1;
	TransferTypes[15].XSide = -1;		//-X    +Z
	TransferTypes[15].ZSide = 1;
	TransferTypes[16].XSide = 1;		//+X    -Z
	TransferTypes[16].ZSide = -1;
	TransferTypes[17].XSide = -1;		//-X    -Z
	TransferTypes[17].ZSide = -1;

	for(int i=0; i<18; i++){

		for(int i2=1; i2<19; i2++){

			bool CondX = (TransferTypes[i].XSide == 0 || BasisVector[i2][0] == TransferTypes[i].XSide);
			bool CondY = (TransferTypes[i].YSide == 0 || BasisVector[i2][1] == TransferTypes[i].YSide);
			bool CondZ = (TransferTypes[i].ZSide == 0 || BasisVector[i2][2] == TransferTypes[i].ZSide);

			if(CondX && CondY && CondZ){
				TransferTypes[i].Vectors[TransferTypes[i].nVectors] = ReflMap[i2];
				TransferTypes[i].nVectors++;
			}

		}

		nTransferNonSolid[i] = 0;		//Number of non-solids in each part of the transfer layer

	}

	nTransferNonSolids = 0;				//Total number of non-solids in transfer layer


	//Find adjoining domains

	CPUDomain* Dm = &ThreadDomain;

	int nAdjoining = 0;
	int c=0;

	for(int Loop=0; Loop<2; Loop++){		//Run through twice - once to find number of adjoining domains, secondly to fill array

		for(int i=0; i<18; i++){	

			int x0 = Dm->x0;	//Ranges sought
			int x1 = Dm->x1;
			int y0 = Dm->y0;
			int y1 = Dm->y1;
			int z0 = Dm->z0;
			int z1 = Dm->z1;

			int x0_Local = Dm->x0;
			int x1_Local = Dm->x1;
			int y0_Local = Dm->y0;
			int y1_Local = Dm->y1;
			int z0_Local = Dm->z0;
			int z1_Local = Dm->z1;

			LatticeToDomain(x0_Local, y0_Local, z0_Local);
			LatticeToDomain(x1_Local, y1_Local, z1_Local);

			if(TransferTypes[i].XSide == 1){
				x0 = Dm->xp;
				x1 = Dm->xp;
				x0_Local = Dm->DomainWidthX;
				x1_Local = Dm->DomainWidthX;
			}
			if(TransferTypes[i].XSide == -1){
				x0 = Dm->xn;
				x1 = Dm->xn;
				x0_Local = -1;
				x1_Local = -1;
			}
			if(TransferTypes[i].YSide == 1){
				y0 = Dm->yp;
				y1 = Dm->yp;
				y0_Local = Dm->DomainWidthY;
				y1_Local = Dm->DomainWidthY;
			}
			if(TransferTypes[i].YSide == -1){
				y0 = Dm->yn;
				y1 = Dm->yn;
				y0_Local = -1;
				y1_Local = -1;
			}
			if(TransferTypes[i].ZSide == 1){
				z0 = Dm->zp;
				z1 = Dm->zp;
				z0_Local = Dm->DomainWidthZ;
				z1_Local = Dm->DomainWidthZ;
			}
			if(TransferTypes[i].ZSide == -1){
				z0 = Dm->zn;
				z1 = Dm->zn;
				z0_Local = -1;
				z1_Local = -1;
			}

			TransferTypes[i].x0 = x0;
			TransferTypes[i].x1 = x1;
			TransferTypes[i].y0 = y0;
			TransferTypes[i].y1 = y1;
			TransferTypes[i].z0 = z0;
			TransferTypes[i].z1 = z1;

			TransferTypes[i].Domainx0 = x0_Local;
			TransferTypes[i].Domainx1 = x1_Local;
			TransferTypes[i].Domainy0 = y0_Local;
			TransferTypes[i].Domainy1 = y1_Local;
			TransferTypes[i].Domainz0 = z0_Local;
			TransferTypes[i].Domainz1 = z1_Local;


			//Create array of solids for each part of the overlap layer

			TransferSolid[i] = new long long[ (x1-x0+1)*(y1-y0+1)*(z1-z0+1) ];
			memset(TransferSolid[i], 0, sizeof(long long)*(x1-x0+1)*(y1-y0+1)*(z1-z0+1));


			//Find adjoining domains
	
			for(int i2=0; i2<ThreadNum; i2++){

				CPUDomain* Dm2 = &Domains[i2];

				bool CondX = (Dm2->x0 <= x1 && Dm2->x1 >= x0) && !(x0 == x1 && (x0 == 0 || x0 == NLatticeX-1) && BoundaryConditionX == 1);
				bool CondY = (Dm2->y0 <= y1 && Dm2->y1 >= y0) && !(y0 == y1 && (y0 == 0 || y0 == NLatticeY-1) && BoundaryConditionY == 1);
				bool CondZ = (Dm2->z0 <= z1 && Dm2->z1 >= z0) && !(z0 == z1 && (z0 == 0 || z0 == NLatticeZ-1) && BoundaryConditionZ == 1);

				if(CondX && CondY && CondZ){	//Overlap of domain on range

					if(i2==ThreadID){			//Deal with loop boundaries in Initialise()
						continue;
					}

					if(Loop==0){
						nAdjoining++;
					}else{

						TransferRequest* Rq = &RecvRequests[c];

						Rq->RequestID = c;
						Rq->Thread = i2;
						Rq->Vectors = i;
						Rq->x0 = max( x0, Dm2->x0 );
						Rq->x1 = min( x1, Dm2->x1 );
						Rq->y0 = max( y0, Dm2->y0 );
						Rq->y1 = min( y1, Dm2->y1 );
						Rq->z0 = max( z0, Dm2->z0 );
						Rq->z1 = min( z1, Dm2->z1 );

						c++;
					}

				}

			}

		}

		if(Loop==0){
			RecvRequests = new TransferRequest[nAdjoining];
			nRecvRequest = nAdjoining;
		}

	}

	//Inform each thread of how many sends are required from it

	nSendRequest = 0;

	for(int tId=0; tId<ThreadNum; tId++){

		if(tId == ThreadID){

			for(int i=0; i<ThreadNum; i++){

				if(i==ThreadID){
					continue;
				}

				int nRecvThread = 0;

				for(int ci=0; ci<nRecvRequest; ci++){
					if(RecvRequests[ci].Thread == i){
						nRecvThread++;
					}
				}

				MPI_Send(&nRecvThread, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

			}

		}else{

			int nSendThread;
			MPI_Recv(&nSendThread, 1, MPI_INT, tId, 0, MPI_COMM_WORLD, &Stat);

			nSendRequest += nSendThread;

		}

	}

	SendRequests = new TransferRequest[nSendRequest];


	//Obtain non-solids from neighbours for each part of the overlap layer

	int SReqIndex = 0;

	for(int tId=0; tId<ThreadNum; tId++){

		if(tId == ThreadID){	//Requester
			
			int Rect[8];

			for(int i=0; i<nRecvRequest; i++){

				Rect[0] = RecvRequests[i].x0;
				Rect[1] = RecvRequests[i].x1;
				Rect[2] = RecvRequests[i].y0;
				Rect[3] = RecvRequests[i].y1;
				Rect[4] = RecvRequests[i].z0;
				Rect[5] = RecvRequests[i].z1;

				Rect[6] = RecvRequests[i].RequestID;
				Rect[7] = RecvRequests[i].Vectors;

				int x0 = Rect[0];		//Index of first x in range
				int x1 = Rect[1];		//Index of last (inclusive) x in range
				int y0 = Rect[2];
				int y1 = Rect[3];
				int z0 = Rect[4];
				int z1 = Rect[5];

				int DatLength = (x1-x0+1)*(y1-y0+1)*(z1-z0+1);
				bool* SDat = new bool[DatLength];

				MPI_Send(Rect, 8, MPI_INT, RecvRequests[i].Thread, 0, MPI_COMM_WORLD);

				MPI_Recv(SDat, DatLength, MPI_BYTE, RecvRequests[i].Thread, 0, MPI_COMM_WORLD, &Stat);

				int nRqNonSolids = 0;

				int c=0;
				for(int z=z0; z<=z1; z++){
				for(int y=y0; y<=y1; y++){
				for(int x=x0; x<=x1; x++){

					int x_ = x - TransferTypes[RecvRequests[i].Vectors].x0;
					int y_ = y - TransferTypes[RecvRequests[i].Vectors].y0;
					int z_ = z - TransferTypes[RecvRequests[i].Vectors].z0;

					int nX = TransferTypes[RecvRequests[i].Vectors].x1 - TransferTypes[RecvRequests[i].Vectors].x0 + 1;
					int nY = TransferTypes[RecvRequests[i].Vectors].y1 - TransferTypes[RecvRequests[i].Vectors].y0 + 1;

					int Ind = nX*(nY*z_ + y_) + x_;

					TransferSolid[RecvRequests[i].Vectors][Ind] = SDat[c] ? 0 : 1;	//0 if solid, 1 if non-solid

					if(!SDat[c]){	//Non-solid
						nTransferNonSolid[RecvRequests[i].Vectors]++;
						nTransferNonSolids++;
						nRqNonSolids++;
					}

					c++;
				}
				}
				}

				delete[] SDat;

				RecvRequests[i].NonSolidCount = nRqNonSolids;

				RecvRequests[i].Data = new double[TransferTypes[RecvRequests[i].Vectors].nVectors * nRqNonSolids];

			}

			//Inform all threads to continue

			Rect[0] = -1;
			for(int i=0; i<ThreadNum; i++){
				if(i!=tId){
					MPI_Send(Rect, 8, MPI_INT, i, 0, MPI_COMM_WORLD);
				}
			}

		}else{					//Servant

			while(true){

				int Rect[8];		//Range on lattice thread wants solids from {x0,x1,y0,y1,z0,z1}
				MPI_Recv(Rect, 8, MPI_INT, tId, 0, MPI_COMM_WORLD, &Stat);

				if(Rect[0]==-1){	//Finished or no data needed from current thread
					break;
				}

				int x0 = Rect[0];		//Index of first x in range
				int x1 = Rect[1];		//Index of last (inclusive) x in range
				int y0 = Rect[2];
				int y1 = Rect[3];
				int z0 = Rect[4];
				int z1 = Rect[5];

				int DatLength = (x1-x0+1)*(y1-y0+1)*(z1-z0+1);
				bool* SDat = new bool[DatLength];				//Data to send back

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				int nRqNonSolids = 0;

				int c=0;
				for(int z=z0; z<=z1; z++){
				for(int y=y0; y<=y1; y++){
				for(int x=x0; x<=x1; x++){

					SDat[c] = (GetSolid(x,y,z)==0);

					if(!SDat[c]){			//Non-solid
						nRqNonSolids++;
					}

					c++;

				}
				}
				}

				MPI_Send(SDat, DatLength, MPI_BYTE, tId, 0, MPI_COMM_WORLD);

				delete[] SDat;

				//Add entry to send requests

				SendRequests[SReqIndex].Thread = tId;	//Send to current master thread

				SendRequests[SReqIndex].x0 = Rect[0];
				SendRequests[SReqIndex].x1 = Rect[1];
				SendRequests[SReqIndex].y0 = Rect[2];
				SendRequests[SReqIndex].y1 = Rect[3];
				SendRequests[SReqIndex].z0 = Rect[4];
				SendRequests[SReqIndex].z1 = Rect[5];

				SendRequests[SReqIndex].RequestID = Rect[6];
				SendRequests[SReqIndex].Vectors = Rect[7];

				SendRequests[SReqIndex].NonSolidCount = nRqNonSolids;

				SendRequests[SReqIndex].Data = new double[TransferTypes[Rect[7]].nVectors * nRqNonSolids];

				SReqIndex++;

			}

		}

	}

}


int LatticeDecomposition(int ThreadNum, int ThreadID){

	ofstream DecompositionLog;
	DecompositionLog.open(OutputPath(DecompositionLogFile));

	long long TotalNonSolids;										//Total number of fluid nodes

	long long nVoxels = ((long long)NLatticeX)*((long long)NLatticeY)*((long long)NLatticeZ);

	int NArray = (int)(nVoxels>>3) + 1;
	_DecompositionLattice = new char[NArray];						//Bit array of solids

	memset(_DecompositionLattice,0,NArray);

	bool ReadStat;

	if(GeometryFileBin){											//Read in geometry to solids array
		ReadStat = DecompositionReadGeometryBin(&TotalNonSolids);
	}else{
		ReadStat = DecompositionReadGeometry(&TotalNonSolids);
	}

	if(!ReadStat){
		return 1;
	}

	double MemRq = (double)(sizeof(long long)*nVoxels + sizeof(Node)*TotalNonSolids) / 1073741824.0;	//Gigabytes

	printf("Memory required: %.2fGB\n", MemRq);
	
	double ycoefficient = (double)NLatticeY / (double)NLatticeX;
	double zcoefficient = (double)NLatticeZ / (double)NLatticeX;
	double N = pow((double)ThreadNum / (ycoefficient*zcoefficient) , 1.0/3.0);

	int nZ = (int)_round( zcoefficient*N );							//Number of partitions in Z

	int* NonSolidsZSlices = new int[NLatticeZ];						//Count number of non solids in each slice in z
	int* NonSolidsZCumulative = new int[NLatticeZ];					//Cumulative number of non solids in z

	for(int z=0; z!=NLatticeZ; z++){
		NonSolidsZSlices[z]=0;

		for(int x=0; x!=NLatticeX; x++){
		for(int y=0; y!=NLatticeY; y++){

			if(!_DecompositionGetLattice(x,y,z)){
				NonSolidsZSlices[z]++;
			}

		}
		}

		if(z!=0){
			NonSolidsZCumulative[z] = NonSolidsZCumulative[z-1] + NonSolidsZSlices[z];
		}else{
			NonSolidsZCumulative[0] = NonSolidsZSlices[0];
		}

	}

	double Porosity = ((double)TotalNonSolids)/((double)(NLatticeX*NLatticeY*NLatticeZ));
	double NonSolidsPerCPU = ((double)TotalNonSolids)/((double)ThreadNum);

	DecompositionLog << "Porosity: " << Porosity << endl;
	DecompositionLog << "Non-solids per CPU: " << NonSolidsPerCPU << endl << endl;

	//Divide in Z

	LatticeBlock* Plane = new LatticeBlock[nZ];			//Divided into blocks in Z according to nZ

	for(int i=0; i<nZ; i++){
		Plane[i].NSquares = ThreadNum/nZ + ( (i < ThreadNum%nZ) ? 1 : 0 );
	}

	DecompositionLog << "Dividing z into " << nZ << " partitions of average " << ((double)TotalNonSolids)/((double)nZ) << " non-solids" << endl;

		//Find precise range
	
	int z0 = 0;
	double z0_ = 0;

	int z1 = 0;
	double z1_ = 0;

	double NSolidsPreceding = 0;

	for(int i=0; i<nZ; i++){

		double NSolidsInSlice = ((double)(Plane[i].NSquares) * NonSolidsPerCPU);

		while(true){

			double N0 = 0;
			if(z1!=0){
				N0 = NonSolidsZCumulative[z1-1] - NSolidsPreceding;
			}
			double N1 = NonSolidsZCumulative[z1] - NSolidsPreceding;

			if(z1==NLatticeZ){

				z1_ = 0;
				break;

			}else if(N0<NSolidsInSlice && N1>=NSolidsInSlice){

				z1_ = (NSolidsInSlice - N0)/NonSolidsZSlices[z1];
				break;

			}

			z1++;
		}

		NSolidsPreceding = NonSolidsZCumulative[z1-1] + z1_*NonSolidsZSlices[z1];

		Plane[i].c0 = z0;
		Plane[i].c0_ = z0_;
		Plane[i].c1 = z1;
		Plane[i].c1_ = z1_;
		Plane[i].n_d = ((double)z1 + z1_) - ((double)z0 + z0_);

		DecompositionLog << "Block range: " << ((double)z0 + z0_) << " -> " << ((double)z1 + z1_) << " with " << NSolidsInSlice << " non-solids (" << Plane[i].NSquares << " CPUs)" << endl;

		z0 = z1;
		z0_ = z1_;

	}

		//Set rounded range

	if(nZ>1){

		DecompositionMinimiseError(Plane, nZ, NonSolidsZCumulative, NLatticeZ, NonSolidsPerCPU, DecompositionLog);

	}else{

		Plane[0].c0Rounded = 0;
		Plane[0].c1Rounded = NLatticeZ;

	}

	//For each z block, divide into x,y cuboids

	for(int zBlock=0; zBlock<nZ; zBlock++){

		double ycoefficient = (double)NLatticeY / (double)NLatticeX;
		double N = sqrt( (double)Plane[zBlock].NSquares/ycoefficient );

		int nY = (int)_round( ycoefficient*N );	//Number of blocks in y


		int* NonSolidsYPlane = new int[NLatticeY];			//Count number of non solids in each slice in y of this z block
		int* NonSolidsYCumulative = new int[NLatticeY];		//Cumulative number of non solids in z
		int TotalNonSolidsY = 0;

		Plane[zBlock].SubBlocks = new LatticeBlock[nY];
		Plane[zBlock].SubBlockNum = nY;

		for(int i=0; i<nY; i++){
			Plane[zBlock].SubBlocks[i].NSquares = Plane[zBlock].NSquares/nY + ( (i < Plane[zBlock].NSquares%nY) ? 1 : 0 );
		}

		for(int y=0; y<NLatticeY; y++){

			NonSolidsYPlane[y] = 0;

			for(int z=Plane[zBlock].c0Rounded; z<Plane[zBlock].c1Rounded; z++){
			for(int x=0; x<NLatticeX; x++){

				if(!_DecompositionGetLattice(x,y,z)){
					NonSolidsYPlane[y]++;
					TotalNonSolidsY++;
				}

			}
			}

			NonSolidsYCumulative[y] = TotalNonSolidsY;

		}

		DecompositionLog << endl;
		DecompositionLog << "z partition " << zBlock << " containing " << TotalNonSolidsY << " non-solids" << endl;

		double NonSolidsPerCPU = ((double)TotalNonSolidsY)/((double)Plane[zBlock].NSquares);
		double NSolidsPerSlice = ((double)TotalNonSolidsY)/((double)nY);

		DecompositionLog << "Non-solids per CPU: " << NonSolidsPerCPU << endl;
		DecompositionLog << '\t' << "Dividing y into " << nY << " partitions of average " << NSolidsPerSlice << " non-solids" << endl;

		int y0 = 0;
		double y0_ = 0;

		int y1 = 0;
		double y1_ = 0;

		double NSolidsPreceding = 0;

		for(int i=0; i<nY; i++){

			double NSolidsInSlice = (((double)(Plane[zBlock].SubBlocks[i].NSquares)) * NonSolidsPerCPU);

			while(true){

				double N0 = (y1!=0) ? (NonSolidsYCumulative[y1-1] - NSolidsPreceding) : 0;
				double N1 = NonSolidsYCumulative[y1] - NSolidsPreceding;

				if(y1==NLatticeY){

					y1_ = 0;
					break;

				}else if(N0<NSolidsInSlice && N1>=NSolidsInSlice){

					y1_ = (NSolidsInSlice - N0)/NonSolidsYPlane[y1];
					break;

				}

				y1++;
			}

			NSolidsPreceding = NonSolidsYCumulative[y1-1] + y1_*NonSolidsYPlane[y1];

			Plane[zBlock].SubBlocks[i].c0 = y0;
			Plane[zBlock].SubBlocks[i].c0_ = y0_;
			Plane[zBlock].SubBlocks[i].c1 = y1;
			Plane[zBlock].SubBlocks[i].c1_ = y1_;
			Plane[zBlock].SubBlocks[i].n_d = ((double)y1 + y1_) - ((double)y0 + y0_);

			DecompositionLog << '\t' << "Block range: " << ((double)y0 + y0_) << " -> " << ((double)y1 + y1_) << " with " << NSolidsInSlice << " non-solids (" << Plane[zBlock].SubBlocks[i].NSquares << " CPUs)" << endl;
			
			y0 = y1;
			y0_ = y1_;

		}
		
		if(nY>1){

			DecompositionMinimiseError(Plane[zBlock].SubBlocks, nY, NonSolidsYCumulative, NLatticeY, NonSolidsPerCPU, DecompositionLog);

		}else{		

			Plane[zBlock].SubBlocks[0].c0Rounded = 0;
			Plane[zBlock].SubBlocks[0].c1Rounded = NLatticeY;

		}

		//For each y block, divide in x to obtain decomposed domain

		for(int yBlock = 0; yBlock<nY; yBlock++){

			LatticeBlock* Block = &(Plane[zBlock].SubBlocks[yBlock]);

			int nX = Block->NSquares;	//Number of blocks in x

			int* NonSolidsXPlane = new int[NLatticeX];			//Count number of non solids in each slice in x of this y block
			int* NonSolidsXCumulative = new int[NLatticeX];		//Cumulative number of non solids in x
			int TotalNonSolidsX = 0;

			Block->SubBlocks = new LatticeBlock[nX];
			Block->SubBlockNum = nX;

			for(int i=0; i<nX; i++){
				Block->SubBlocks[i].NSquares = 1;
			}

			for(int x=0; x<NLatticeX; x++){

				NonSolidsXPlane[x] = 0;

				for(int z=Plane[zBlock].c0Rounded; z<Plane[zBlock].c1Rounded; z++){
				for(int y=Block->c0Rounded; y<Block->c1Rounded; y++){

					if(!_DecompositionGetLattice(x,y,z)){
						NonSolidsXPlane[x]++;
						TotalNonSolidsX++;
					}

				}
				}

				NonSolidsXCumulative[x] = TotalNonSolidsX;

			}

			DecompositionLog << endl;
			DecompositionLog << '\t' << "y partition " << yBlock << " containing " << TotalNonSolidsX << " non-solids" << endl;

			double NonSolidsPerCPU = ((double)TotalNonSolidsX)/((double)Block->NSquares);
			double NSolidsPerSlice = ((double)TotalNonSolidsX)/((double)nX);

			DecompositionLog << "\tNon-solids per CPU: " << NonSolidsPerCPU << endl;
			DecompositionLog << '\t' << "Dividing x into " << nX << " partitions of average " << NSolidsPerSlice << " non-solids" << endl;

			int x0 = 0;
			double x0_ = 0;

			int x1 = 0;
			double x1_ = 0;

			double NSolidsPreceding = 0;

			for(int i=0; i<nX; i++){

				while(true){

					double N0 = 0;
					if(x1!=0){
						N0 = NonSolidsXCumulative[x1-1] - NSolidsPreceding;
					}
					double N1 = NonSolidsXCumulative[x1] - NSolidsPreceding;

					if(x1==NLatticeX){

						x1_ = 0;
						break;

					}else if(N0<NonSolidsPerCPU && N1>=NonSolidsPerCPU){

						x1_ = (NonSolidsPerCPU - N0)/NonSolidsXPlane[x1];
						break;

					}

					x1++;
				}

				NSolidsPreceding = NonSolidsXCumulative[x1-1] + x1_*NonSolidsXPlane[x1];

				Block->SubBlocks[i].c0 = x0;
				Block->SubBlocks[i].c0_ = x0_;
				Block->SubBlocks[i].c1 = x1;
				Block->SubBlocks[i].c1_ = x1_;
				Block->SubBlocks[i].n_d = ((double)x1 + x1_) - ((double)x0 + x0_);

				DecompositionLog << "\t\t" << "Block range: " << ((double)x0 + x0_) << " -> " << ((double)x1 + x1_) << " with " << NonSolidsPerCPU << " non-solids" << endl;

				x0 = x1;
				x0_ = x1_;

			}

			if(nX>1){

				DecompositionMinimiseError(Plane[zBlock].SubBlocks[yBlock].SubBlocks, nX, NonSolidsXCumulative, NLatticeX, NonSolidsPerCPU, DecompositionLog);

			}else{

				Block->SubBlocks[0].c0Rounded = 0;
				Block->SubBlocks[0].c1Rounded = NLatticeX;

			}

			delete[] NonSolidsXPlane;
			delete[] NonSolidsXCumulative;

		}

		delete[] NonSolidsYPlane;
		delete[] NonSolidsYCumulative;

	}

	delete[] NonSolidsZSlices;
	delete[] NonSolidsZCumulative;


	//Organise CPU domains

	Domains = new CPUDomain[ThreadNum];		//Thread 0 domains

	int clength[6] = {0,0,0,0,0,0};

	int i=0;
	for(int z=0; z<nZ; z++){

		LatticeBlock* ZBlock = &(Plane[z]);

		for(int y=0; y<ZBlock->SubBlockNum; y++){

			LatticeBlock* YBlock = &(ZBlock->SubBlocks[y]);

			for(int x=0; x<YBlock->SubBlockNum; x++){

				LatticeBlock* XBlock = &(YBlock->SubBlocks[x]);

				Domains[i].x0 = XBlock->c0Rounded;
				Domains[i].x1 = XBlock->c1Rounded-1;
				Domains[i].y0 = YBlock->c0Rounded;
				Domains[i].y1 = YBlock->c1Rounded-1;
				Domains[i].z0 = ZBlock->c0Rounded;
				Domains[i].z1 = ZBlock->c1Rounded-1;

				Domains[i].xp = (Domains[i].x1==NLatticeX-1) ? 0 : Domains[i].x1+1;
				Domains[i].xn = (Domains[i].x0==0) ? NLatticeX-1 : Domains[i].x0-1;
				Domains[i].yp = (Domains[i].y1==NLatticeY-1) ? 0 : Domains[i].y1+1;
				Domains[i].yn = (Domains[i].y0==0) ? NLatticeY-1 : Domains[i].y0-1;
				Domains[i].zp = (Domains[i].z1==NLatticeZ-1) ? 0 : Domains[i].z1+1;
				Domains[i].zn = (Domains[i].z0==0) ? NLatticeZ-1 : Domains[i].z0-1;

				Domains[i].DomainWidthX = Domains[i].x1 - Domains[i].x0 + 1;
				Domains[i].DomainWidthY = Domains[i].y1 - Domains[i].y0 + 1;
				Domains[i].DomainWidthZ = Domains[i].z1 - Domains[i].z0 + 1;

				//

				int tlength[6] = {1,1,1,1,1,1};
				for(int c=10; Domains[i].x0 > c; c*=10){ tlength[0]++; }
				for(int c=10; Domains[i].x1 > c; c*=10){ tlength[1]++; }
				for(int c=10; Domains[i].y0 > c; c*=10){ tlength[2]++; }
				for(int c=10; Domains[i].y1 > c; c*=10){ tlength[3]++; }
				for(int c=10; Domains[i].z0 > c; c*=10){ tlength[4]++; }
				for(int c=10; Domains[i].z1 > c; c*=10){ tlength[5]++; }

				clength[0] = max( clength[0], tlength[0] );
				clength[1] = max( clength[1], tlength[1] );
				clength[2] = max( clength[2], tlength[2] );
				clength[3] = max( clength[3], tlength[3] );
				clength[4] = max( clength[4], tlength[4] );
				clength[5] = max( clength[5], tlength[5] );

				//

				i++;

			}

			delete[] YBlock->SubBlocks;

		}

		delete[] ZBlock->SubBlocks;

	}

	DecompositionLog << endl << "Domain ranges:" << endl;

	int il = 1;
	for(int c=10; (ThreadNum-1) > c; c*=10){ il++; }

	char HChar[96];
	memset(HChar, ' ', sizeof(HChar));

	HChar[4 + il + clength[0]] = 'x';
	HChar[10 + il + clength[0] + clength[1] + clength[2]] = 'y';
	HChar[16 + il + clength[0] + clength[1] + clength[2] + clength[3] + clength[4]] = 'z';
	HChar[16 + il + clength[0] + clength[1] + clength[2] + clength[3] + clength[4] + 1] = '\0';

	DecompositionLog << HChar << endl;

	for(int i=0; i<ThreadNum; i++){

		char CoordChar[96];
		memset(CoordChar, ' ', sizeof(CoordChar));

		int Ind  = 0;

		_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "(%d)         ", i);
		Ind += il + 3;

		_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "%d           ", Domains[i].x0);
		Ind += clength[0];
		_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, " -> %d       ", Domains[i].x1);
		Ind += clength[1] + 6;
		_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "%d           ", Domains[i].y0);
		Ind += clength[2];
		_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, " -> %d       ", Domains[i].y1);
		Ind += clength[3] + 6;
		_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, "%d           ", Domains[i].z0);
		Ind += clength[4];
		_snprintf(&CoordChar[Ind], sizeof(CoordChar)-Ind, " -> %d\0     ", Domains[i].z1);
		Ind += clength[5] + 6;

		DecompositionLog << CoordChar << endl;
		

	}

	DecompositionLog.close();

	OutputDecompositionGeometry(ThreadNum);

	return 0;

}

void MPIDistributeDecomposition(int* ret, int ThreadID, int ThreadNum){

	MPI_Status Stat;


	//Distribute return flag for file read-in

	MPI_Bcast(ret, 1, MPI_INT, 0, MPI_COMM_WORLD);	

	if((*ret)!=0){
		goto MPIDist_End;
	}


	//Distribute domain decomposition data

	if(ThreadID!=0){
		Domains = new CPUDomain[ThreadNum];
	}

	MPI_Bcast(Domains, sizeof(CPUDomain)*ThreadNum, MPI_BYTE, 0, MPI_COMM_WORLD);

	memcpy(&ThreadDomain, &Domains[ThreadID], sizeof(CPUDomain));	//Thread's domain


	//Distribute solids

	long long nNonSolidTot = 0;

	if(ThreadID==0){

		for(int tId=1; tId<ThreadNum; tId++){

			CPUDomain* Dm = &Domains[tId];

			int nNodes = Dm->DomainWidthX * Dm->DomainWidthY * Dm->DomainWidthZ;

			bool* Geometry = new bool[nNodes];

			if(Geometry==0){
				cout << "Decomposition memory allocation error (0)" << endl;
				return;
			}

			//Obtain geometry

			int c=0;
			for(int z=0; z<Dm->DomainWidthZ; z++){
			for(int y=0; y<Dm->DomainWidthY; y++){
			for(int x=0; x<Dm->DomainWidthX; x++){

				int Vx = Dm->x0 + x;
				int Vy = Dm->y0 + y;
				int Vz = Dm->z0 + z;

				Geometry[c] = _DecompositionGetLattice(Vx, Vy, Vz);

				if(!Geometry[c]){
					nNonSolidTot++;
				}

				c++;
			}
			}
			}

			//Send thread geometry

			int dCount = 0;
			int dLength = 1048576;		//1MB

			while(dCount < nNodes){
				int l = min( dLength, nNodes-dCount );

				MPI_Send(&Geometry[dCount], l, MPI_BYTE, tId, 0, MPI_COMM_WORLD);			//Send geometry to thread

				dCount += l;
			}

			delete[] Geometry;

		}

		//Thread 0 geometry

		CPUDomain* Dm = &ThreadDomain;

		int nNodes = Dm->DomainWidthX * Dm->DomainWidthY * Dm->DomainWidthZ;

		//Create solids array

		Solid = new long long[nNodes];

		if(Solid==0){
			cout << "Decomposition memory allocation error (2)" << endl;
			return;
		}

		long long nInd = 1;

		int c=0;
		for(int z=0; z<Dm->DomainWidthZ; z++){
		for(int y=0; y<Dm->DomainWidthY; y++){
		for(int x=0; x<Dm->DomainWidthX; x++){

			int Vx = Dm->x0 + x;
			int Vy = Dm->y0 + y;
			int Vz = Dm->z0 + z;

			if(_DecompositionGetLattice(Vx, Vy, Vz)){	//Solid
				Solid[c] = 0;
			}else{
				Solid[c] = nInd;
				nInd++;
				nNonSolidTot++;
			}

			c++;
		}
		}
		}

		nNonSolid = nInd-1;

	}else{

		CPUDomain* Dm = &ThreadDomain;

		int nNodes = Dm->DomainWidthX * Dm->DomainWidthY * Dm->DomainWidthZ;

		bool* Geometry = new bool[nNodes];

		if(Geometry==0){
			cout << "Decomposition memory allocation error (1)" << endl;
			return;
		}

		//Receive geometry

		int dCount = 0;
		int dLength = 1048576;		//1MB

		while(dCount < nNodes){
			int l = min( dLength, nNodes-dCount );

			MPI_Recv(&Geometry[dCount], l, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);			//Send geometry to thread

			dCount += l;
		}


		//Create solids array

		Solid = new long long[nNodes];

		if(Solid==0){
			cout << "Decomposition memory allocation error (2)" << endl;
			return;
		}

		long long nInd = 1;

		for(int i=0; i<nNodes; i++){

			if(Geometry[i]){	//Solid

				Solid[i] = 0;

			}else{

				Solid[i] = nInd;
				nInd++;

			}

		}

		nNonSolid = nInd-1;

		delete[] Geometry;

	}

	//Distribute total number of non-solids

	MPI_Bcast(&nNonSolidTot, sizeof(long long), MPI_BYTE, 0, MPI_COMM_WORLD);

	nNonSolidLattice = nNonSolidTot;

	//

MPIDist_End:

	delete[] _DecompositionLattice;
	_DecompositionLattice = 0;

	return;
}

void OutputDecompositionGeometry(int ThreadNum){

	if(!OutputDecompGeometry || _DecompositionLattice==0){
		return;
	}

	FILE* OutFile = fopen(OutputPath(DecompositionGeometryFile), "wb");

	if(!OutFile){
		cout << "Unable to output decomposed geometry to \"" << OutputPath(DecompositionGeometryFile) << "\"" << endl;
		return;
	}

	cout << "Outputting decomposed geometry VTK to \"" << DecompositionGeometryFile << "\"" << endl;

	//Output header

	char Header[2048];
	int Ind = 0;

	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "# vtk DataFile Version 2.0\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "Decomposed Geometry\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ASCII\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DATASET STRUCTURED_POINTS\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "DIMENSIONS %i %i %i\n", NLatticeX, NLatticeY, NLatticeZ);
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "ORIGIN 0 0 0\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SPACING 1 1 1\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "POINT_DATA %i\n", NLatticeX*NLatticeY*NLatticeZ);
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "SCALARS sample_scalars float\n");
	Ind += _snprintf(&Header[Ind], sizeof(Header)-Ind, "LOOKUP_TABLE default\n");

	fwrite(Header, sizeof(char), Ind, OutFile);


	//Output geometry

	char c0[2] = {'0', ' '};
	char c1[16];

	int c1L = 0;

	for(int z=0; z<NLatticeZ; z++){
	for(int y=0; y<NLatticeY; y++){
	for(int x=0; x<NLatticeX;    ){

		int x1;

		for(int i=0; i<ThreadNum; i++){
			if(x>=Domains[i].x0&&x<=Domains[i].x1 && y>=Domains[i].y0&&y<=Domains[i].y1 && z>=Domains[i].z0&&z<=Domains[i].z1){
				x1 = Domains[i].x1;
				c1L = _snprintf(c1, 16, "%d ", i+1);
				break;
			}
		}

		while(x <= x1){

			if(_DecompositionGetLattice(x,y,z)){	//Solid
				fwrite(c1, sizeof(char), c1L, OutFile);
			}else{
				fwrite(c0, sizeof(char), 2, OutFile);
			}

			x++;
		}

	}
	}
	}

	fclose(OutFile);

}

bool DecompositionReadGeometryBin(long long* NumberOfNonSolids){

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

			_DecompositionSetLattice(x, y, z, true);

		}else{					//Non solid

			_DecompositionSetLattice(x, y, z, false);

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

bool DecompositionReadGeometry(long long* NumberOfNonSolids){

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

				_DecompositionSetLattice(x, y, z, true);

			}else{					//Non solid

				_DecompositionSetLattice(x, y, z, false);
				
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

void DecompositionMinimiseError(LatticeBlock* Plane, int nBlocks, int* nCumulative, int NLattice, double NonSolidsPerCPU, ofstream& DecompositionLog){

	//Obtain minimum error combination of rounded borders

	int nBorders = nBlocks - 1;

	unsigned long long FlagInt = (0xFFFFFFFFFFFFFFFF << nBorders);	//1111111111....1110000 (nZ rightmost 0s)

	unsigned long long MinErrorCombination;
	double MinError = DBL_MAX;

	while(FlagInt!=0){	//Wait until it overflows

		unsigned long long Bits = FlagInt;

		double Error = 0;
		for(int i=0; i!=nBorders; i++){

			double NSolidsInSlice = ((double)(Plane[i].NSquares) * NonSolidsPerCPU);

			bool Round = (bool)(Bits&1);	//1 or 0

			Plane[i].UpperRound = Round;
			Plane[i+1].LowerRound = Round;

			Bits >>= 1;

			if(i==0){

				int rIndex;
				
				if(Round){
					rIndex = (int)ceil(Plane[0].c1_);
				}else{
					rIndex = (int)floor(Plane[0].c1_);
				}

				int zIndex = Plane[0].c1 + rIndex;

				double delta = nCumulative[zIndex-1] - NSolidsInSlice;

				Error += (delta*delta);


				if(nBorders==1){	//2 Partitions

					NSolidsInSlice = ((double)(Plane[1].NSquares) * NonSolidsPerCPU);

					if(Plane[0].UpperRound){
						rIndex = (int)ceil(Plane[1].c0_);
					}else{
						rIndex = (int)floor(Plane[1].c0_);
					}

					zIndex = Plane[1].c0 + rIndex;

					delta = (nCumulative[NLattice-1] - nCumulative[zIndex-1]) - NSolidsInSlice;

					Error += (delta*delta);

				}

			}else{

				int rIndex0;
				if(Plane[i].LowerRound){
					rIndex0 = (int)ceil(Plane[i].c0_);
				}else{
					rIndex0 = (int)floor(Plane[i].c0_);
				}

				int rIndex1;
				if(Plane[i].UpperRound){
					rIndex1 = (int)ceil(Plane[i].c1_);
				}else{
					rIndex1 = (int)floor(Plane[i].c1_);
				}

				int zIndex0 = Plane[i].c0 + rIndex0;
				int zIndex1 = Plane[i].c1 + rIndex1;

				double delta = (nCumulative[zIndex1-1] - nCumulative[zIndex0-1]) - NSolidsInSlice;

				Error += (delta*delta);


				if(i==(nBorders-1)){

					NSolidsInSlice = ((double)(Plane[nBorders].NSquares) * NonSolidsPerCPU);

					if(Plane[i].UpperRound){
						rIndex0 = (int)ceil(Plane[nBorders].c0_);
					}else{
						rIndex0 = (int)floor(Plane[nBorders].c0_);
					}

					zIndex0 = Plane[nBorders].c0 + rIndex0;

					delta = (nCumulative[NLattice-1] - nCumulative[zIndex0-1]) - NSolidsInSlice;

					Error += (delta*delta);

				}

			}

		}

		if(Error < MinError){

			MinErrorCombination = FlagInt;
			MinError = Error;

		}

		FlagInt++;
	}

	//Set lowest error combination

//	DecompositionLog << endl;
//	DecompositionLog << "Rounding borders to the minimum error combination: " << endl;

	unsigned long long Bits = MinErrorCombination;

	for(int i=0; i!=nBorders; i++){

		double NSolidsInSlice = ((double)(Plane[i].NSquares) * NonSolidsPerCPU);

		bool Round = (bool)(Bits&1);	//1 or 0

		Plane[i].UpperRound = Round;
		Plane[i+1].LowerRound = Round;

		Bits >>= 1;

		if(i==0){

			int rIndex;
				
			if(Round){
				rIndex = (int)ceil(Plane[0].c1_);
			}else{
				rIndex = (int)floor(Plane[0].c1_);
			}

			int zIndex = Plane[0].c1 + rIndex;

			double delta = nCumulative[zIndex-1] - NSolidsInSlice;
			double Err = fabs(delta)/NSolidsInSlice;

			Plane[0].c0Rounded = 0;
			Plane[0].c1Rounded = zIndex;


			if(nBorders==1){		//2 Partitions

				NSolidsInSlice = ((double)(Plane[1].NSquares) * NonSolidsPerCPU);

				if(Plane[0].UpperRound){
					rIndex = (int)ceil(Plane[1].c0_);
				}else{
					rIndex = (int)floor(Plane[1].c0_);
				}

				zIndex = Plane[1].c0 + rIndex;

				delta = (nCumulative[NLattice-1] - nCumulative[zIndex-1]) - NSolidsInSlice;
				Err = fabs(delta)/NSolidsInSlice;

				Plane[nBorders].c0Rounded = zIndex;
				Plane[nBorders].c1Rounded = NLattice;

			}

		}else{

			int rIndex0;
			if(Plane[i].LowerRound){
				rIndex0 = (int)ceil(Plane[i].c0_);
			}else{
				rIndex0 = (int)floor(Plane[i].c0_);
			}

			int rIndex1;
			if(Plane[i].UpperRound){
				rIndex1 = (int)ceil(Plane[i].c1_);
			}else{
				rIndex1 = (int)floor(Plane[i].c1_);
			}

			int zIndex0 = Plane[i].c0 + rIndex0;
			int zIndex1 = Plane[i].c1 + rIndex1;

			double delta = (nCumulative[zIndex1-1] - nCumulative[zIndex0-1]) - NSolidsInSlice;
			double Err = fabs(delta)/NSolidsInSlice;

			Plane[i].c0Rounded = zIndex0;
			Plane[i].c1Rounded = zIndex1;

			if(i==(nBorders-1)){

				NSolidsInSlice = ((double)(Plane[nBorders].NSquares) * NonSolidsPerCPU);

				if(Plane[i].UpperRound){
					rIndex0 = (int)ceil(Plane[nBorders].c0_);
				}else{
					rIndex0 = (int)floor(Plane[nBorders].c0_);
				}

				zIndex0 = Plane[nBorders].c0 + rIndex0;

				delta = (nCumulative[NLattice-1] - nCumulative[zIndex0-1]) - NSolidsInSlice;
				Err = fabs(delta)/NSolidsInSlice;

				Plane[nBorders].c0Rounded = zIndex0;
				Plane[nBorders].c1Rounded = NLattice;

			}

		}

	}

}