//Input Parameters:///////////////////////////////////////////////////////////////////////////////////////

#define Using_WinMultiThreading	0		//Windows multithreading (shared memory)
#define WinMultiThreading_TNum	7		//Number of threads to run in this mode

#define Using_MPI 0						//Use parallel computing (1 or 0, not true or false)
#define MPI_Mode 0						//Mode 0 runs full geometry on all CPUs, dividing particles between them, Mode 1 divides up lattice among CPUs
#define MPIMode1_DecomposeZ 0			//Whether or not to decompose the lattice in the z (flow) direction for MPI mode 1

#define NLattice_x 320					//Lattice of NLattice_x x NLattice_y x NLattice_z
#define NLattice_y 320
#define NLattice_z 640

#define BoundaryConditionX 0			//Boundary conditions in X,Y,Z. Loop boundary = 0, Solid boundary = 1
#define BoundaryConditionY 0
#define BoundaryConditionZ 0

#define ParticleNum 0					//Number of tracer particles. Set to 0 to use particles per voxel...
#define ParticlesPerVoxel 1				//For uniform initialisation, number of particles per voxel. Use a cube number eg 1, 8, 27, 64...

#define DiffusionCoefficient 0.00001	//Diffusion coefficient for random walk length

#define SimulationTimestep 0			//Timestep for the simulation. Set to 0 for automatic (optimum) value
#define OutputTimeInterval 1000000		//Time interval for writing output data
#define SimulationTimeMax  1000000000		//End time of simulation. Set to 0 for no limit
#define SimulationTimeScale 10000000		//Value of dt corresponding physically to 1 second

#define LatticeResolutionMicron 4.9			//Resolution of lattice in micrometres

#define RelativePropagatorBinWidth 0.01		//Bin width relative to average flow distance for propagator


//Input folder
#define InputFilesFolder "H:\\2012 ChemEng PhD\\VC++ Projects\\Input Files\\Portland Carbonate\\"
//#define InputFilesFolder "H:\\2012 ChemEng PhD\\VC++ Projects\\Input Files\\Bentheimer Sandstone\\"
//#define InputFilesFolder "H:\\2012 ChemEng PhD\\VC++ Projects\\Input Files\\Sandpack\\"
//#define InputFilesFolder "H:\\2012 ChemEng PhD\\VC++ Projects\\Input Files\\Circular Tunnel\\"
//#define InputFilesFolder "C:\\Users\\Farrel\\Desktop\\Input Files\\Sandpack\\"
//#define InputFilesFolder "C:\\Users\\Farrel\\Desktop\\Input Files\\Bentheimer Sandstone\\"
//#define InputFilesFolder "/home/ps/ce-fs1/fg709/Dispersion/CircularTunnel/Input/"


//Input files
/*
#define SolidsFileName "CircularTunnel a=20.txt"			//File specifying solid voxels
#define VelocitiesFileName "VectorField a=20.vtk"		//File containing vector field
#define SolidsFileName "Bentheimer_Geometry_250x250x500.vtk"			//File specifying solid voxels
#define VelocitiesFileName "Bentheimer_Vectorfield_250x250x500.vtk"		//File containing vector field
#define SolidsFileName "Bentheimer_Geometry_400x400x800.vtk"			//File specifying solid voxels
#define VelocitiesFileName "Bentheimer_Vectorfield_400x400x800.vtk"		//File containing vector field
#define SolidsFileName "Sandpack_Geometry_400x200x200.vtk"			//File specifying solid voxels
#define VelocitiesFileName "Sandpack_VectorField_400x200x200.vtk"		//File containing vector field
*/
#define SolidsFileName "Portland_Geometry_320x320x640.vtk"			//File specifying solid voxels
#define VelocitiesFileName "Portland_Vectorfield_320x320x640.vtk"		//File containing vector field


//Output folder
#define OutputFilesFolder "H:\\2012 ChemEng PhD\\VC++ Projects\\Output Files\\"
//#define OutputFilesFolder "C:\\Users\\Farrel\\Desktop\\Output Files\\"
//#define OutputFilesFolder "/home/ps/ce-fs1/fg709/Dispersion/CircularTunnel/Output/"

//Output files
#define DispersionOutputFile "Dispersion Coefficient Portland.txt"		//File for dispersion coefficient
#define DistributionOutputFile "Particle Distribution.txt"					//File for particle position distributions on lattice. Use %ST for simulation time, %RT for real scaled time
#define PropagatorOutputFile "Propagator Distribution Portland t=%RT.txt"	//File for particle propagator distribution. Use %ST for simulation time, %RT for real scaled time
#define DecompositionLogFile "Decomposition Log.txt"						//File for lattice decomposition log (MPI Mode 1)

/////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////

//Debug mode:	0 = off
//				1 = minimal (catches infinite loops only - fastest)
//				2 = detailed (captures infinite loops and invalidate exit locations and generates debug data - medium)
//				3 = full (catches infinite loops, infinities, invalid exit locations and generates full debug data - slow)

#define DebugAdvectionMode 1	

bool DebugParticle = false;

ofstream AdvectDebug;
#define AdvectionDebugFile "H:\\2012 ChemEng PhD\\VC++ Projects\\Output Files\\AdvectionDebug.txt"
//#define AdvectionDebugFile "C:\\Users\\Farrel\\Desktop\\Output Files\\AdvectionDebug.txt"

//////////////////////////////

#if Using_WinMultiThreading
#include <Windows.h>
#include <process.h>
#endif

#if Using_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <cmath>
#include <limits.h>
#include <math.h>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <float.h>
#include <ctime>
#include <string.h>

using namespace std;

#define PI 3.14159265358979323846264338328
#define Min_dV (1.0E-9)							//Minimum difference between opposing face velocities for which equations hold (otherwise treat as special cases)
#define MIN_VECTOR_RESOLUTION Min_dV			//Vectors read in smaller than this are set to 0
						
int NParticles;							//Number of particles (per thread in MPI)

struct Voxel{				//sizeof(Voxel) = 80
double Vx;					//Voxel vector x
double Vy;					//Voxel vector y
double Vz;					//Voxel vector z
double u1;					//Voxel face velocity x1
double u2;					//Voxel face velocity x2
double v1;					//Voxel face velocity y1
double v2;					//Voxel face velocity y2
double w1;					//Voxel face velocity z1
double w2;					//Voxel face velocity z2
unsigned char Solid;		//bits correspond to neighbouring solid voxels |0|NegativeZSolid|PositiveZSolid|NegativeYSolid|PositiveYSolid|NegativeXSolid|PositiveXSolid|VoxelSolid|
};

//Bits of Solid byte eg if(Solid&VoxelSolid).. or if(Solid&PositiveYSolid).. etc
#define VoxelSolid 0x01
#define PositiveXSolid 0x02
#define NegativeXSolid 0x04
#define PositiveYSolid 0x08
#define NegativeYSolid 0x10
#define PositiveZSolid 0x20
#define NegativeZSolid 0x40

struct Particle{			//sizeof(Particle) = 72 (72 unpadded)
double X;					//X,Y,Z position in whole lattice
double Y;
double Z;
double InitialX;			//Particle starting position
double InitialY;
double InitialZ;
short VoxelX;				//Voxel the particle is in
short VoxelY;
short VoxelZ;
short LastVoxelX;
short LastVoxelY;
short LastVoxelZ;
int XCycles;				//Number of times the particle has run through boundary (+ve in positive direction etc)
int YCycles;
int ZCycles;
};

struct Vector{
double X;
double Y;
double Z;
};

struct Coords{
short x;
short y;
short z;
short xp;
short yp;
short zp;
short xn;
short yn;
short zn;
};

struct LatticeInformation{
	Vector AverageFlowVelocity;
	double MaxFlowVelocityComponent;
	int NNonSolids;
	int NSolids;
	int NVoxels;
	double Porosity;
};

struct CPUDomain{
	short x0;				//Indeces of CPU domain voxel range
	short x1;				//Index of last voxel in x range
	short y0;
	short y1;
	short z0;
	short z1;

	/*
	short xn;				//Coordinate of 1 Lattice padding in negative x
	short xp;				//Coordinate of 1 Lattice padding in positive x
	short yn;				//Coordinate of 1 Lattice padding in negative y
	short yp;				//Coordinate of 1 Lattice padding in positive y
	short zn;				//Coordinate of 1 Lattice padding in negative z
	short zp;				//Coordinate of 1 Lattice padding in positive z
	*/

	short DomainWidthX;		//Dimensions of voxel array
	short DomainWidthY;
	short DomainWidthZ;
		
	int* NeighbourDomainsXn;	//Array of neighbouring domains in xn side
	int nNeighbourDomainsXn;

	int* NeighbourDomainsXp;	//Array of neighbouring domains in xp side
	int nNeighbourDomainsXp;

	int* NeighbourDomainsYn;	//Array of neighbouring domains in yn side
	int nNeighbourDomainsYn;

	int* NeighbourDomainsYp;	//Array of neighbouring domains in yp side
	int nNeighbourDomainsYp;

	int* NeighbourDomainsZn;	//Array of neighbouring domains in zn side
	int nNeighbourDomainsZn;

	int* NeighbourDomainsZp;	//Array of neighbouring domains in zp side
	int nNeighbourDomainsZp;

	LatticeInformation LatticeInfo;	//CPU Domain lattice information

	int NParticlesDomain;		//Number of particles in domain
};

class SimulationTimer{

private:

	clock_t SimulationStart;
	clock_t LastStepTime;

public:

	SimulationTimer(){
		SimulationStart = clock();
	}

	double GetTimeSinceLastStep(){

		clock_t StepTime = clock();

		double dt = ((double)(StepTime - LastStepTime))/CLOCKS_PER_SEC;

		LastStepTime = StepTime;

		return dt;		//dt in seconds
	}

	double GetSimulationTime(){

		clock_t Time = clock();

		return ((double)(Time - SimulationStart))/CLOCKS_PER_SEC;	//time in seconds

	}
	
};

//Class for generating uniform random numbers modified and extended from
//J. D. Cook, “Simple random number generation"
//Accessed May 2012. http://www.codeproject.com/Articles/25172/Simple-Random-Number-Generation.

class RandomNumber{

private:

	unsigned int m_w;
	unsigned int m_z;

	unsigned int GetUint(){

		m_z = 36969 * (m_z & 65535) + (m_z >> 16);
		m_w = 18000 * (m_w & 65535) + (m_w >> 16);

		return (m_z << 16) + m_w;
	}

public:

	RandomNumber(){

		m_w = 521288629;	//Default seeds
		m_z = 362436069;

	}

	void SetSeed(unsigned int mw, unsigned int mz){

		m_w = mw;
		m_z = mz;

		if(m_w==0x464FFFFF||m_w==0x8C9FFFFE||m_w==0xD2EFFFFD){	//w = 0x464FFFFF, 0x8C9FFFFE, or 0xD2EFFFFD fixed points
			m_w = 521288629;	//Default seed
		}

		if(m_z==0x9068FFFF){		//z = 0x9068FFFF fixed point
			m_z = 362436069;
		}

	}

	double GetUniform() {
            // 0 <= u < 2^32
		unsigned int u = GetUint();
            // The magic number below is 1/(2^32 + 2).
            // The result is strictly between 0 and 1.0
		return (u + 1.0) * 2.328306435454494e-10;
	}

	void TimeSeed(){	//Randomise according to local time

		time_t timet;
		tm* TimeInfo;

		time(&timet);
		TimeInfo = localtime(&timet);

		double TFraction = (TimeInfo->tm_sec + TimeInfo->tm_hour + TimeInfo->tm_wday + 1.0)/(61.0 + 23.0 + 6.0 + 2.0);

		long double Tfr = ((long double)UINT_MAX)*TFraction;
		
		m_w = (unsigned int)floor(Tfr);

		if(m_w==0x464FFFFF||m_w==0x8C9FFFFE||m_w==0xD2EFFFFD){	//w = 0x464FFFFF, 0x8C9FFFFE, or 0xD2EFFFFD fixed points
			m_w = 521288629;	//Default seed
		}

		long double Tfr2 = ((long double)UINT_MAX)*(TFraction*46.4 + TimeInfo->tm_sec + 1.0)/(1.0 + 62.0 + 46.4);

		m_z = (unsigned int)floor(Tfr2);

		if(m_z==0x9068FFFF){		//z = 0x9068FFFF fixed point
			m_z = 362436069;
		}

		GetUniform();

	}

	void ThreadSeed(int ThreadID, int ThreadNum){	//Randomise according to local time and thread id

		time_t timet;
		tm* TimeInfo;

		time(&timet);
		TimeInfo = localtime(&timet);

		double TFraction = (TimeInfo->tm_sec + TimeInfo->tm_hour + TimeInfo->tm_wday + 1.0)/(61.0 + 23.0 + 6.0 + 2.0);

		long double Tfr = ((long double)UINT_MAX)*TFraction;
		
		m_w = (unsigned int)floor(Tfr);

		if(m_w==0x464FFFFF||m_w==0x8C9FFFFE||m_w==0xD2EFFFFD){	//w = 0x464FFFFF, 0x8C9FFFFE, or 0xD2EFFFFD fixed points
			m_w = 521288629;	//Default seed
		}
		
		long double Tfr2 = ((long double)UINT_MAX)*(ThreadID + 1.0 + TFraction)/(ThreadNum + 1.1);

		m_z = (unsigned int)floor(Tfr2);

		if(m_z==0x9068FFFF){		//z = 0x9068FFFF fixed point
			m_z = 362436069;
		}

		GetUniform();

	}

};

RandomNumber Random;								//Initialise class called Random. Get next number with: Random.GetUniform();
SimulationTimer SimulationTime;						//Simulation timer gives time between calculation steps

void CreateLattice();													//Sets lattice entirely non solid, all flow vectors to 0
void ReadVelocities();													//Sets lattice velocity field
void ReadSolids(unsigned int* NumberOfNonSolids);						//Reads in the solids file for non-decomposed lattice
void InitialiseLattice();												//Calculates and sets voxel face velocities
void InitialiseParticles();												//Called from non-mpi main to set initial particle positions
void InitialiseParticlesAtPoint(double x, double y, double z);			//Can be called from InitialiseParticles() to set all particles at a single point (x,y,z)
void InitialiseParticlesRandomXY(int z, int npreceding);				//Distributes particles randomly in the x-y plane at a given z
void InitialiseParticlesUniform();										//Disitrbutes particles evenly throughout the lattice
void SetLatticeInfo();													//Calculates average velocities, porosity etc in LatticeInfo struct
void WriteOutParameters();												//Writes out average velocities, porosity, Peclet number to the console
void RunTimeStep(double dt);											//Carries out an advection and diffusion step on all particles
void RunTimeStep(double dt, int ParticleIndex0, int ParticleIndex1);	//Carries out an advection and diffusion step on particles between and including given indices
void GetCoordinates(short x, short y, short z, Coords* C);
void RandomWalkParticle(int i, double dt, double RandomWalkLength);
void VoxelAfterDisplacement(Particle* P, Vector* V, Coords* C);
void AdvectParticle(int i, double dt);
void CaseNoSolids(int i,Coords* C, Voxel* LatticeElement, double* dt);
void CaseOneSolid(int i,Coords* C, Voxel* LatticeElement, double* dt);
void CaseTwoSolids(int i,Coords* C, Voxel* LatticeElement, double* dt);
void CaseThreeSolids(int i,Coords* C, Voxel* LatticeElement, double* dt);
void CaseFourSolids(int i,Coords* C, Voxel* LatticeElement, double* dt);
void UpdateParticleHistory(int i);
void OutputParticleDistribution(int XBinsPerVoxel, int YBinsPerVoxel, int ZBinsPerVoxel);	//Non-MPI outputs particle distribution on the lattice (not propagator)

void CreateLattice(unsigned int** LatticeBase_, Voxel** Lattice_);
void ReadVelocities(unsigned int* LatticeBase_, Voxel* Lattice_);					//Sets lattice velocity field
void ReadSolids(unsigned int* NumberOfNonSolids, unsigned int* LatticeBase_);		//Reads in the solids file for non-decomposed lattice
void InitialiseLattice(unsigned int* LatticeBase_, Voxel* Lattice_);	

#if Using_MPI	//MPI Functions

void InitialiseParticlesMPI(int n, int ThreadID, int ThreadNum, int* ThreadParticleCount, int ThreadParticleIndex);		//Initialises particles for each processor
void SendDispersionStatsToThread0();																					//Calculates particle statistics for CPU's particles and sends the data to thread 0 to calculate the coefficient
double Thread0CalculateDispersionCoefficient(int ThreadNum, double* LastVariance, double t, double dt, int LogInterval, ofstream& outfile);		//Receives other threads' particle statistics and calculates the dispersion coefficient and writes it to file
void SendPropagatorStatsToThread0(double RelativeBinWidth, double t);					//Calculates propagator distribution for thread's particles and sends them to thread 0
void Thread0OutputPropagatorGraph(double RelativeBinWidth, double t, int ThreadNum);	//Receives threads' propagator data and writes to file

int LatticeDecomposition(int NProcessors, int ThreadID);				//Thread 0 decomposes lattice and sends data to other threads - returns 0 on success
int ReceiveLatticeDecompositionData(int NProcessors, int ThreadID);		//Threads receive lattice decomposition data
void _DecompositionGetVectorField(ifstream& InFile, __int64 FileLength, __int64* StreamPtr, double* VectorField, CPUDomain* Domain);	//Gets domain's vector field including overlap
void _DecompositionGetVectorField(ifstream& InFile, __int64 FileLength, __int64* StreamPtr, Voxel* VectorField, CPUDomain* Domain);		//Gets domain's vector field including overlap
unsigned char _DecompositionGetVoxelSolid(int x, int y, int z);										//Obtains voxel solid flags
void _DecompositionGetStreamPtrArray(ifstream& InFile, __int64* FileLength, __int64* StreamPtr);	//Create array of filestream pointers to z,y coordinates in file
void DecompositionInitialiseLattice(int NProcessors, int ThreadID);		//Initialise face velocities amongst threads
void DomainInitialiseParticlesUniform(int NProcessors, int ThreadID, int* NParticlesThread);	//MPI Lattice division uniform particle initialisation

struct ParticleStatusData;

void DomainRunTimeStep(double dt);																//Carries out an advection and diffusion step on all particles in domain
void DomainRandomWalkParticle(int i, double dt, double RandomWalkLength, ParticleStatusData* ParticleStatus);	//Random walk particle in domain
void DomainAdvectParticle(int i, double dt, ParticleStatusData* ParticleStatus);								//Advects particle in domain
void DomainUpdateParticleHistory(int i, ParticleStatusData* ParticleStatus);									//Update particle struct information and domain transfer
void DomainGetCoordinates(short x, short y, short z, Coords* C);
void DomainVoxelAfterDisplacement(Particle* P, Vector* V, Coords* C);

#endif

unsigned int* LatticeBase;				//Array for sparse storage. Element = 0 for solid node; Element = [index] for non solid, where [index] is the index of the voxel struct in Lattice array
Voxel* Lattice;							//Array of lattice elements. Element 0 is a solid voxel

#if Using_MPI && (MPI_Mode == 1)

struct CPUDomain;

CPUDomain ProcessorDomain;				//Current thread's domain
CPUDomain* Domains;						//Array of all threads' domains

int _CPUDomainWidthX;					//Number of X elements in array
int _CPUDomainWidthY;					//Number of Y elements in array

int ParticleArrayLength;				//Number of elements in particle array

/*	With overlap
inline void DomainToLattice(short& x, short& y, short& z){
	x += ProcessorDomain.xn;
	y += ProcessorDomain.yn;
	z += ProcessorDomain.zn;
	while(x>=NLattice_x){ x-=NLattice_x; }
	while(y>=NLattice_y){ y-=NLattice_y; }
	while(z>=NLattice_z){ z-=NLattice_z; }
}

inline void LatticeToDomain(short& x, short& y, short& z){
	x -= ProcessorDomain.xn;
	y -= ProcessorDomain.yn;
	z -= ProcessorDomain.zn;
	while(x<0){ x+=NLattice_x; }
	while(y<0){ y+=NLattice_y; }
	while(z<0){ z+=NLattice_z; }
}
*/

inline void DomainToLattice(double& x, double& y, double& z){
	x += (double)ProcessorDomain.x0;
	y += (double)ProcessorDomain.y0;
	z += (double)ProcessorDomain.z0;
}

inline void LatticeToDomain(double& x, double& y, double& z){
	x -= (double)ProcessorDomain.x0;
	y -= (double)ProcessorDomain.y0;
	z -= (double)ProcessorDomain.z0;
}

inline void DomainToLattice(short& x, short& y, short& z){
	x += ProcessorDomain.x0;
	y += ProcessorDomain.y0;
	z += ProcessorDomain.z0;
}

inline void LatticeToDomain(short& x, short& y, short& z){
	x -= ProcessorDomain.x0;
	y -= ProcessorDomain.y0;
	z -= ProcessorDomain.z0;
}

inline bool IsElementSolid(short x, short y, short z){
	return (LatticeBase[_CPUDomainWidthX*(_CPUDomainWidthY*z + y) + x]==0);
}

inline Voxel* GetLattice(int x, int y, int z){
	return &(Lattice[LatticeBase[_CPUDomainWidthX*(_CPUDomainWidthY*z + y) + x]]);
}

#else

inline bool IsElementSolid(short x, short y, short z){
	return (LatticeBase[NLattice_z*(NLattice_y*x + y) + z]==0);
}

Voxel* GetLattice(int x, int y, int z){
	return &(Lattice[LatticeBase[NLattice_z*(NLattice_y*x + y) + z]]);
}

#endif

Particle* Particles;				//Dynamic array to allocate number of particles assigned to each thread

LatticeInformation LatticeInfo;		//Struct to be filled with (full) lattice information



// MPI Mode 1 Lattice Decomposition declarations

char* _DecompositionLattice;		//For MPI Mode 1, temporary geometry array

bool _DecompositionGetLattice(int x, int y, int z){

	int Index = NLattice_z*(NLattice_y*x + y) + z;
	int rShift = Index & 0x00000007;	//Keep last 3 bits
	int ArrIndex = Index>>3;			// divide by 8	

	return (_DecompositionLattice[ArrIndex]>>rShift)&0x01;
}

void _DecompositionSetLattice(int x, int y, int z, bool Value){

	int Index = NLattice_z*(NLattice_y*x + y) + z;
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

struct LatticeZBlock;
struct LatticeYBlock;
struct LatticeXBlock;

struct LatticeZBlock{
	int z0;				//Lower bound z integer part
	double z0_;			//Lower bound z fractional part
	int z1;				//Upper bound z integer part
	double z1_;			//Upper bound z fractional part

	bool LowerRound;	//Round lower border down or up
	bool UpperRound;	//Round upper border down or up

	int z0Rounded;		//Rounded lower border
	int z1Rounded;		//Rounded upper border

	double nZ_d;		//Exact number of z slices in plane

	int NSquares;		//Number of processor blocks in plane

	LatticeYBlock* YBlocks;
	int YBlockNum;
};

struct LatticeYBlock{
	int y0;				//Lower bound y integer part
	double y0_;			//Lower bound y fractional part
	int y1;				//Upper bound y integer part
	double y1_;			//Upper bound y fractional part

	bool LowerRound;	//Round lower border down or up
	bool UpperRound;	//Round upper border down or up

	int y0Rounded;		//Rounded lower border
	int y1Rounded;		//Rounded upper border

	double nY_d;		//Exact number of y slices in plane

	int NSquares;		//Number of processor blocks in plane

	LatticeXBlock* XBlocks;
	int XBlockNum;
};

struct LatticeXBlock{
	int x0;				//Lower bound y integer part
	double x0_;			//Lower bound y fractional part
	int x1;				//Upper bound y integer part
	double x1_;			//Upper bound y fractional part

	bool LowerRound;	//Round lower border down or up
	bool UpperRound;	//Round upper border down or up

	int x0Rounded;		//Rounded lower border
	int x1Rounded;		//Rounded upper border

	double nX_d;		//Exact number of y slices in plane
};

//

char* OutputFilePath(char* FileName, double t){

	char* Folder = OutputFilesFolder;

	const int MaxPathLength = 1024;

	static char Path[MaxPathLength];

	int Folderl = (int)strlen(Folder);
	int Filel = (int)strlen(FileName);

	if(Folderl==0){
		return FileName;
	}

	bool fs = true;

	int PathIndex = 0;
	int i = 0;
	while(i!=Folderl){

		Path[PathIndex] = Folder[i];
		PathIndex++;
		if(PathIndex == MaxPathLength){return 0;}

		if(Folder[i]=='\\'){
			fs = false;
		}

		i++;
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

	i=0;
	while(i!=Filel){
		if(FileName[i]=='%' && i<(Filel-2)){
			if(FileName[i+1]=='S' && FileName[i+2]=='T'){			//Print simulation time

				int n = _snprintf(&Path[PathIndex],MaxPathLength-PathIndex,"%.0f",t);

				i+=3;
				PathIndex+=n;
				if(PathIndex >= MaxPathLength){return 0;}
				continue;

			}else if(FileName[i+1]=='R' && FileName[i+2]=='T'){		//Print scaled time

				double tscaled = t / (double)SimulationTimeScale;
				double tscinterval = (double)OutputTimeInterval / (double)SimulationTimeScale;

				int ndp = 0;
				if(tscinterval < 1.0){
					ndp = 1;
				}
				if(tscinterval < 0.1){
					ndp = 2;
				}
				if(tscinterval < 0.01){
					ndp = 3;
				}
				if(tscinterval < 0.001){
					ndp = 4;
				}

				int n = _snprintf(&Path[PathIndex],MaxPathLength-PathIndex,"%.*f",ndp,tscaled);

				i+=3;
				PathIndex+=n;
				if(PathIndex >= MaxPathLength){return 0;}
				continue;

			}
		}

		Path[PathIndex] = FileName[i];
		PathIndex++;
		if(PathIndex >= MaxPathLength){return 0;}

		i++;
	}

	Path[PathIndex] = '\0';

	return Path;

}

char* OutputFilePath(char* FileName){

	return OutputFilePath(FileName, 0);

}

char* InputFilePath(char* FileName){

	char* Folder = InputFilesFolder;

	const int MaxPathLength = 1024;

	static char Path[MaxPathLength];

	int Folderl = (int)strlen(Folder);
	int Filel = (int)strlen(FileName);

	if(Folderl==0){
		return FileName;
	}

	bool fs = true;

	int PathIndex = 0;
	int i = 0;
	while(i!=Folderl){

		Path[PathIndex] = Folder[i];
		PathIndex++;
		if(PathIndex == MaxPathLength){return 0;}

		if(Folder[i]=='\\'){
			fs = false;
		}

		i++;
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

	i=0;
	while(i!=Filel){

		Path[PathIndex] = FileName[i];
		PathIndex++;
		if(PathIndex >= MaxPathLength){return 0;}

		i++;
	}

	Path[PathIndex] = '\0';

	return Path;

}

void ReadSolids(unsigned int* NumberOfNonSolids, unsigned int* LatticeBase_){

	ifstream InFile;
	InFile.open(InputFilePath(SolidsFileName),ios::binary);
	
	if(InFile!=0){
		cout << "Reading in solids from '" << SolidsFileName << "'" << endl;
	}else{
		cout << "Could not open solids file" << endl;
		return;
	}

	InFile.seekg(0,ios::end);
	__int64 flength = InFile.tellg();	//File length
	InFile.seekg(0,ios::beg);

	const __int64 buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	
	__int64 nRead = 0;				//Number of chars read from file in total
	__int64 readl = 0;				//Number of chars read from file into buffer
	__int64 readpos = 0;			//Read position in buffer

	int nNonSolids = 0;

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				cout << "Error: reached end of geometry file before lattice entirely read" << endl;
				goto ReadSolidsReturnRJMP;
			}
			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]=='1'||buf[readpos]=='0'){

			if(buf[readpos]=='1'){	//Solid

				LatticeBase_[NLattice_z*(NLattice_y*x + y) + z] = 0;

			}else{	//Non solid

				LatticeBase_[NLattice_z*(NLattice_y*x + y) + z] = nNonSolids+1;
				nNonSolids++;
			}

			x++;
			if(x==NLattice_x){
				x=0;
				y++;
				if(y==NLattice_y){
					y=0;
					z++;
					if(z==NLattice_z){
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
	InFile.close();
	return;
}

void ReadSolids(unsigned int* NumberOfNonSolids){

	ifstream InFile;
	InFile.open(InputFilePath(SolidsFileName),ios::binary);
	
	if(InFile!=0){
		cout << "Reading in solids from '" << SolidsFileName << "'" << endl;
	}else{
		cout << "Could not open solids file" << endl;
		return;
	}

	InFile.seekg(0,ios::end);
	__int64 flength = InFile.tellg();	//File length
	InFile.seekg(0,ios::beg);

	const __int64 buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	
	__int64 nRead = 0;				//Number of chars read from file in total
	__int64 readl = 0;				//Number of chars read from file into buffer
	__int64 readpos = 0;			//Read position in buffer

	int nNonSolids = 0;

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				cout << "Error: reached end of geometry file before lattice entirely read" << endl;
				goto ReadSolidsReturnRJMP;
			}
			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]=='1'||buf[readpos]=='0'){

			//Voxel* LatticeElement = GetLattice(x,y,z);

			if(buf[readpos]=='1'){	//Solid

				LatticeBase[NLattice_z*(NLattice_y*x + y) + z] = 0;

			}else{	//Non solid

				LatticeBase[NLattice_z*(NLattice_y*x + y) + z] = nNonSolids+1;
				nNonSolids++;
			}

			x++;
			if(x==NLattice_x){
				x=0;
				y++;
				if(y==NLattice_y){
					y=0;
					z++;
					if(z==NLattice_z){
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
	InFile.close();
	return;
}

void ReadVelocities(unsigned int* LatticeBase_, Voxel* Lattice_){

	ifstream InFile(InputFilePath(VelocitiesFileName),ios::binary);

	if(InFile!=0){
		cout << "Reading in velocity field from '" << VelocitiesFileName << "'" << endl;
	}else{
		cout << "Could not open velocity field file" << endl;
		return;
	}

	InFile.seekg(0,ios::end);
	__int64 flength = InFile.tellg();	//File length
	InFile.seekg(0,ios::beg);

	char initChar;
	InFile.read(&initChar,1);

	InFile.seekg(0,ios::beg);

	const __int64 buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	int component=0;

	bool state = 1;											//Read char (true) or whitespace/return (false)
	if(initChar==' '||initChar=='\n'||initChar=='\r'){		//Between values
		state = 0;
	}

	char VBuffer[24];				//Buffer value
	int VBufl = 0;

	double v[3];					//Vector
	
	__int64 nRead = 0;				//Number of chars read from file in total
	__int64 readl = 0;				//Number of chars read from file into buffer
	__int64 readpos = 0;			//Read position in buffer

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				delete[] buf;
				return;
			}
			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]==' '||buf[readpos]=='\n'||buf[readpos]=='\r'){	//Between values

			if(VBufl!=0){

				VBuffer[VBufl]='\0';
				v[component] = atof(VBuffer);
						
				if(fabs(v[component]) < MIN_VECTOR_RESOLUTION){
					v[component] = 0;
				}

				VBufl = 0;
			}

			state = 0;

		}else if(state==0){							//Beginning of next value

			state = 1;

			VBuffer[VBufl] = buf[readpos];
			VBufl++;

			component++;
			if(component==3){

				if(LatticeBase_[NLattice_z*(NLattice_y*x + y) + z]!=0){

					Voxel* LatticeElement = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*x + y) + z]]);

					LatticeElement->Vx = v[0];
					LatticeElement->Vy = v[1];
					LatticeElement->Vz = v[2];

				}
				
			component=0;
			x++;
			if(x==NLattice_x){
				x=0;
				y++;
				if(y==NLattice_y){
					y=0;
					z++;
					if(z==NLattice_z){
						delete[] buf;
						return;
					}
				}
			}

			}

		}else{		//Number 

			VBuffer[VBufl] = buf[readpos];
			VBufl++;

		}

		readpos++;
	}

}

void ReadVelocities(){

	ifstream InFile(InputFilePath(VelocitiesFileName),ios::binary);

	if(InFile!=0){
		cout << "Reading in velocity field from '" << VelocitiesFileName << "'" << endl;
	}else{
		cout << "Could not open velocity field file" << endl;
		return;
	}

	InFile.seekg(0,ios::end);
	__int64 flength = InFile.tellg();	//File length
	InFile.seekg(0,ios::beg);

	char initChar;
	InFile.read(&initChar,1);

	InFile.seekg(0,ios::beg);

	const __int64 buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	int component=0;

	bool state = 1;											//Read char (true) or whitespace/return (false)
	if(initChar==' '||initChar=='\n'||initChar=='\r'){		//Between values
		state = 0;
	}

	char VBuffer[24];				//Buffer value
	int VBufl = 0;

	double v[3];					//Vector
	
	__int64 nRead = 0;				//Number of chars read from file in total
	__int64 readl = 0;				//Number of chars read from file into buffer
	__int64 readpos = 0;			//Read position in buffer

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				delete[] buf;
				return;
			}
			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]==' '||buf[readpos]=='\n'||buf[readpos]=='\r'){	//Between values

			if(VBufl!=0){

				VBuffer[VBufl]='\0';
				v[component] = atof(VBuffer);
						
				if(fabs(v[component]) < MIN_VECTOR_RESOLUTION){
					v[component] = 0;
				}

				VBufl = 0;
			}

			state = 0;

		}else if(state==0){							//Beginning of next value

			state = 1;

			VBuffer[VBufl] = buf[readpos];
			VBufl++;

			component++;
			if(component==3){

				if(!IsElementSolid(x,y,z)){

					Voxel* LatticeElement = GetLattice(x,y,z);

					LatticeElement->Vx = v[0];
					LatticeElement->Vy = v[1];
					LatticeElement->Vz = v[2];

				}
				
			component=0;
			x++;
			if(x==NLattice_x){
				x=0;
				y++;
				if(y==NLattice_y){
					y=0;
					z++;
					if(z==NLattice_z){
						delete[] buf;
						return;
					}
				}
			}

			}

		}else{		//Number 

			VBuffer[VBufl] = buf[readpos];
			VBufl++;

		}

		readpos++;
	}

}

void CreateLattice(unsigned int** LatticeBase__, Voxel** Lattice__){		//Allocate lattice and zero

	*LatticeBase__ = new unsigned int[NLattice_x*NLattice_y*NLattice_z];	//0 for a non solid, array index for solid

	unsigned int* LatticeBase_ = *LatticeBase__;

	unsigned int nNonSolids;						//Obtain number of non solids

	ReadSolids(&nNonSolids, LatticeBase_);			//Read in solids file

	*Lattice__ = new Voxel[nNonSolids+1];	//Array of non solid voxels + element 0 (solid)

	Voxel* Lattice_ = *Lattice__;

	if(Lattice_==0){
		cout << "Lattice memory allocation failed." << endl;
		return;
	}

	Lattice_[0].Solid = VoxelSolid;
	Lattice_[0].Vx = 0;
	Lattice_[0].Vy = 0;
	Lattice_[0].Vz = 0;

	//Zero Lattice

	int x=0;
	while(x!=NLattice_x){
	int y=0;
	while(y!=NLattice_y){
	int z=0;
	while(z!=NLattice_z){

	if(LatticeBase_[NLattice_z*(NLattice_y*x + y) + z]!=0){

		Voxel* LatticeElement = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*x + y) + z]]);

		LatticeElement->Solid = 0;

		LatticeElement->Vx = 0.0;
		LatticeElement->Vy = 0.0;
		LatticeElement->Vz = 0.0;	//0.00035625;

	}

	z++;
	}
	y++;
	}
	x++;
	}

	ReadVelocities(LatticeBase_, Lattice_);					//Read in velocities

	InitialiseLattice(LatticeBase_, Lattice_);				//Calculate face velocities

}

void CreateLattice(){		//Allocate lattice and zero

	LatticeBase = new unsigned int[NLattice_x*NLattice_y*NLattice_z];	//0 for a non solid, array index for solid

	unsigned int nNonSolids;						//Obtain number of non solids

	ReadSolids(&nNonSolids);			//Read in solids file

	Lattice = new Voxel[nNonSolids+1];	//Array of non solid voxels + element 0 (solid)

	if(Lattice==0){
		cout << "Lattice memory allocation failed." << endl;
		return;
	}

	Lattice[0].Solid = VoxelSolid;
	Lattice[0].Vx = 0;
	Lattice[0].Vy = 0;
	Lattice[0].Vz = 0;

	//Zero Lattice

	int x=0;
	while(x!=NLattice_x){
	int y=0;
	while(y!=NLattice_y){
	int z=0;
	while(z!=NLattice_z){

	if(!IsElementSolid(x,y,z)){

		Voxel* LatticeElement = GetLattice(x,y,z);

		LatticeElement->Solid = 0;

		LatticeElement->Vx = 0.0;
		LatticeElement->Vy = 0.0;
		LatticeElement->Vz = 0.0;	//0.00035625;

	}

	z++;
	}
	y++;
	}
	x++;
	}

	ReadVelocities();					//Read in velocities

	InitialiseLattice();				//Calculate face velocities

}

void InitialiseLattice(unsigned int* LatticeBase_, Voxel* Lattice_){

	int x=0;
	while(x!=NLattice_x){
	int y=0;
	while(y!=NLattice_y){
	int z=0;
	while(z!=NLattice_z){

		Coords C;
		GetCoordinates(x,y,z,&C);

		int xp = C.xp;				//Voxel to the right
		int xn = C.xn;				//Voxel to the left
		int yp = C.yp;				
		int yn = C.yn;				
		int zp = C.zp;				
		int zn = C.zn;

		if(LatticeBase_[NLattice_z*(NLattice_y*x + y) + z]!=0){

		Voxel* LatticeElement = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*x + y) + z]]);

		//Calculate velocities at the faces
		Voxel* LatticeElement_xn = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*xn + y) + z]]);
		LatticeElement->u1 = 0.5*(LatticeElement->Vx + LatticeElement_xn->Vx);	//Velocity of negative x face

			if(LatticeElement_xn->Solid&VoxelSolid  || (BoundaryConditionX==1 && x==0)){
				LatticeElement->u1=0;
				LatticeElement->Solid |= NegativeXSolid;
			}

		Voxel* LatticeElement_xp = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*xp + y) + z]]);
		LatticeElement->u2 = 0.5*(LatticeElement->Vx + LatticeElement_xp->Vx);	//Velocity of positive x face

			if(LatticeElement_xp->Solid&VoxelSolid || (BoundaryConditionX==1 && x==(NLattice_x-1))){
				LatticeElement->u2=0;
				LatticeElement->Solid |= PositiveXSolid;
			}


		Voxel* LatticeElement_yn = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*x + yn) + z]]);
		LatticeElement->v1 = 0.5*(LatticeElement->Vy + LatticeElement_yn->Vy);	//Velocity of negative y face

			if(LatticeElement_yn->Solid&VoxelSolid || (BoundaryConditionY==1 && y==0)){
				LatticeElement->v1=0;
				LatticeElement->Solid |= NegativeYSolid;
			}

		Voxel* LatticeElement_yp = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*x + yp) + z]]);
		LatticeElement->v2 = 0.5*(LatticeElement->Vy + LatticeElement_yp->Vy);	//Velocity of positive y face

			if(LatticeElement_yp->Solid&VoxelSolid || (BoundaryConditionY==1 && y==(NLattice_y-1))){
				LatticeElement->v2=0;
				LatticeElement->Solid |= PositiveYSolid;
			}

		Voxel* LatticeElement_zn = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*x + y) + zn]]);
		LatticeElement->w1 = 0.5*(LatticeElement->Vz + LatticeElement_zn->Vz);	//Velocity of negative z face

			if(LatticeElement_zn->Solid&VoxelSolid || (BoundaryConditionZ==1 && z==0)){
				LatticeElement->w1=0;
				LatticeElement->Solid |= NegativeZSolid;
			}

		Voxel* LatticeElement_zp = &(Lattice_[LatticeBase_[NLattice_z*(NLattice_y*x + y) + zp]]);
		LatticeElement->w2 = 0.5*(LatticeElement->Vz + LatticeElement_zp->Vz);	//Velocity of positive z face

			if(LatticeElement_zp->Solid&VoxelSolid || (BoundaryConditionZ==1 && z==(NLattice_z-1))){
				LatticeElement->w2=0;
				LatticeElement->Solid |= PositiveZSolid;
			}

		}

	z++;
	}
	y++;
	}
	x++;
	}

}

void InitialiseLattice(){

	int x=0;
	while(x!=NLattice_x){
	int y=0;
	while(y!=NLattice_y){
	int z=0;
	while(z!=NLattice_z){

		Coords C;
		GetCoordinates(x,y,z,&C);

		int xp = C.xp;				//Voxel to the right
		int xn = C.xn;				//Voxel to the left
		int yp = C.yp;				
		int yn = C.yn;				
		int zp = C.zp;				
		int zn = C.zn;

		if(!IsElementSolid(x,y,z)){

		Voxel* LatticeElement = GetLattice(x,y,z);

		//Calculate velocities at the faces
		Voxel* LatticeElement_xn = GetLattice(xn,y,z);
		LatticeElement->u1 = 0.5*(LatticeElement->Vx + LatticeElement_xn->Vx);	//Velocity of negative x face

			if(LatticeElement_xn->Solid&VoxelSolid  || (BoundaryConditionX==1 && x==0)){
				LatticeElement->u1=0;
				LatticeElement->Solid |= NegativeXSolid;
			}

		Voxel* LatticeElement_xp = GetLattice(xp,y,z);
		LatticeElement->u2 = 0.5*(LatticeElement->Vx + LatticeElement_xp->Vx);	//Velocity of positive x face

			if(LatticeElement_xp->Solid&VoxelSolid || (BoundaryConditionX==1 && x==(NLattice_x-1))){
				LatticeElement->u2=0;
				LatticeElement->Solid |= PositiveXSolid;
			}


		Voxel* LatticeElement_yn = GetLattice(x,yn,z);
		LatticeElement->v1 = 0.5*(LatticeElement->Vy + LatticeElement_yn->Vy);	//Velocity of negative y face

			if(LatticeElement_yn->Solid&VoxelSolid || (BoundaryConditionY==1 && y==0)){
				LatticeElement->v1=0;
				LatticeElement->Solid |= NegativeYSolid;
			}

		Voxel* LatticeElement_yp = GetLattice(x,yp,z);
		LatticeElement->v2 = 0.5*(LatticeElement->Vy + LatticeElement_yp->Vy);	//Velocity of positive y face

			if(LatticeElement_yp->Solid&VoxelSolid || (BoundaryConditionY==1 && y==(NLattice_y-1))){
				LatticeElement->v2=0;
				LatticeElement->Solid |= PositiveYSolid;
			}

		Voxel* LatticeElement_zn = GetLattice(x,y,zn);
		LatticeElement->w1 = 0.5*(LatticeElement->Vz + LatticeElement_zn->Vz);	//Velocity of negative z face

			if(LatticeElement_zn->Solid&VoxelSolid || (BoundaryConditionZ==1 && z==0)){
				LatticeElement->w1=0;
				LatticeElement->Solid |= NegativeZSolid;
			}

		Voxel* LatticeElement_zp = GetLattice(x,y,zp);
		LatticeElement->w2 = 0.5*(LatticeElement->Vz + LatticeElement_zp->Vz);	//Velocity of positive z face

			if(LatticeElement_zp->Solid&VoxelSolid || (BoundaryConditionZ==1 && z==(NLattice_z-1))){
				LatticeElement->w2=0;
				LatticeElement->Solid |= PositiveZSolid;
			}

		}

	z++;
	}
	y++;
	}
	x++;
	}

}

void InitialiseParticlesAtPoint(double x, double y, double z){

	int i=0;
	while(i!=NParticles){
		
		Particles[i].X = x;
		Particles[i].Y = y;
		Particles[i].Z = z;

		Particles[i].VoxelX = (short)floor(Particles[i].X);
		Particles[i].VoxelY = (short)floor(Particles[i].Y);
		Particles[i].VoxelZ = (short)floor(Particles[i].Z);

		Particles[i].InitialX = Particles[i].X;
		Particles[i].InitialY = Particles[i].Y;
		Particles[i].InitialZ = Particles[i].Z;

		Particles[i].LastVoxelX = Particles[i].VoxelX;
		Particles[i].LastVoxelY = Particles[i].VoxelY;
		Particles[i].LastVoxelZ = Particles[i].VoxelZ;

		Particles[i].XCycles = 0;
		Particles[i].YCycles = 0;
		Particles[i].ZCycles = 0;

		i++;
	}

}

void InitialiseParticles(){			//Initialise particles (not using MPI)
	
	if(ParticleNum == 0){	//Number of particles set per voxel

		NParticles = LatticeInfo.NNonSolids * ParticlesPerVoxel;
		//cout << "With " << LatticeInfo.NNonSolids << " non-solids, and " << ParticlesPerVoxel << " particles per voxel, NParticles = " << NParticles << endl;
	}else{

		NParticles = ParticleNum;

	}

	Particles = new Particle[NParticles];

	if(Particles == 0){
		cout << "Memory allocation failure for particle array";
	}

	cout << "Number of particles: " << NParticles << endl;

	InitialiseParticlesUniform();
	//InitialiseParticlesRandomXY(0,0);
	//InitialiseParticlesAtPoint(1.5,1.5,0);

}

void InitialiseParticlesUniform(){		//Distribute particles uniformly throughout the lattice

	double pSpacing = pow( (double)LatticeInfo.NNonSolids , 1.0/3.0 ) / pow( (double)NParticles , 1.0/3.0 );	//Distance between particles

	cout << "Particle spacing: " << pSpacing << "lu" << endl;

	double Px = pSpacing/2;
	double Py = pSpacing/2;
	double Pz = pSpacing/2;

	int i=0;
	while(i!=NParticles){

		Voxel* LatticeElement = GetLattice((int)floor(Px),(int)floor(Py),(int)floor(Pz));
				
		if(!(LatticeElement->Solid&VoxelSolid)){

			Particles[i].X = Px;
			Particles[i].Y = Py;
			Particles[i].Z = Pz;

			Particles[i].VoxelX = (short)floor(Particles[i].X);
			Particles[i].VoxelY = (short)floor(Particles[i].Y);
			Particles[i].VoxelZ = (short)floor(Particles[i].Z);

			Particles[i].InitialX = Particles[i].X;
			Particles[i].InitialY = Particles[i].Y;
			Particles[i].InitialZ = Particles[i].Z;

			Particles[i].LastVoxelX = Particles[i].VoxelX;
			Particles[i].LastVoxelY = Particles[i].VoxelY;
			Particles[i].LastVoxelZ = Particles[i].VoxelZ;

			Particles[i].XCycles = 0;
			Particles[i].YCycles = 0;
			Particles[i].ZCycles = 0;

			i++;
		}

		Px += pSpacing;

		if(Px >= (double)NLattice_x){

			Px = pSpacing/2;
			Py += pSpacing;

			if(Py >= (double)NLattice_y){

				Py = pSpacing/2;
				Pz += pSpacing;

				if(Pz >= (double)NLattice_z){

					Pz = pSpacing/2;

				}
			}
		}

	}

}

void InitialiseParticlesRandomXY(int z, int npreceding){	//Distributes particles randomly about x-y plane at given z starting from voxel (initx, inity)

	int i=0;
	int count=0;

		int x=0;
		while(true){

			int y=0;
			while(y!=NLattice_y){

				Voxel* LatticeElement = GetLattice(x,y,z);
				
				if(!(LatticeElement->Solid&VoxelSolid)){

					if(count<npreceding){
						count++;
					}else{

					double px = ((double)x + Random.GetUniform());
					double py = ((double)y + Random.GetUniform());

					Particles[i].X = px;
					Particles[i].Y = py;
					Particles[i].Z = (double)z;

					Particles[i].VoxelX = (short)floor(Particles[i].X);
					Particles[i].VoxelY = (short)floor(Particles[i].Y);
					Particles[i].VoxelZ = (short)floor(Particles[i].Z);

					Particles[i].InitialX = Particles[i].X;
					Particles[i].InitialY = Particles[i].Y;
					Particles[i].InitialZ = Particles[i].Z;

					Particles[i].LastVoxelX = Particles[i].VoxelX;
					Particles[i].LastVoxelY = Particles[i].VoxelY;
					Particles[i].LastVoxelZ = Particles[i].VoxelZ;

					Particles[i].XCycles = 0;
					Particles[i].YCycles = 0;
					Particles[i].ZCycles = 0;

					i++;

					if(i>=NParticles){
						return;
					}

					}

				}

				y++;
			}

			x++;

			if(x==NLattice_x){x=0;}
		}

}

void SetLatticeInfo(){

	double VMax = 0;

	Vector V = {0,0,0};
	int ncount=0;

	int x=0;
	while(x!=NLattice_x){
	int y=0;
	while(y!=NLattice_y){
	int z=0;
	while(z!=NLattice_z){

		Voxel* LatticeElement = GetLattice(x,y,z);

		if(!(LatticeElement->Solid&VoxelSolid)){
			V.X += LatticeElement->Vx;
			V.Y += LatticeElement->Vy;
			V.Z += LatticeElement->Vz;
			ncount++;

			if(fabs(LatticeElement->Vx) > VMax){
				VMax = fabs(LatticeElement->Vx);
			}
			if(fabs(LatticeElement->Vy) > VMax){
				VMax = fabs(LatticeElement->Vy);
			}
			if(fabs(LatticeElement->Vz) > VMax){
				VMax = fabs(LatticeElement->Vz);
			}
		}

	z++;
	}
	y++;
	}
	x++;
	}

	V.X /= ((double)ncount);
	V.Y /= ((double)ncount);
	V.Z /= ((double)ncount);

	LatticeInfo.AverageFlowVelocity.X = V.X;
	LatticeInfo.AverageFlowVelocity.Y = V.Y;
	LatticeInfo.AverageFlowVelocity.Z = V.Z;

	LatticeInfo.MaxFlowVelocityComponent = VMax;

	LatticeInfo.NSolids = (NLattice_x*NLattice_y*NLattice_z) - ncount;
	LatticeInfo.NNonSolids = ncount;

	//cout << "From " << (NLattice_x*NLattice_y*NLattice_z) << ", " << ncount << " non solids" << endl;

	LatticeInfo.Porosity = ((double)ncount)/((double)(NLattice_x*NLattice_y*NLattice_z));
	
}

#if Using_MPI	//MPI Functions

struct ParticleStatusData{
	int ParticleStatus;		//Direction particle has moved into neighbouring domain xn=1,xp=2,yn=3,yp=4,zn=5,zp=6
	double dt;				//Remaining dt to be carried out in new domain
	Vector RandomWalk;		//Remaining random walk vector to be carried out in new domain
};

void DomainRunTimeStep(double dt){ 

	double RandomWalkLength = sqrt(6*DiffusionCoefficient*dt);	//Computationally intensive, precalculate where possible

	ParticleStatusData PStat;

	int i=0;
	while(i!=NParticles){

		DomainAdvectParticle(i, dt, &PStat);						//Moves particle along vector field

		if(PStat.ParticleStatus==0){
		DomainRandomWalkParticle(i, dt, RandomWalkLength, &PStat);	//Moves particle along a random path
		}

		DomainUpdateParticleHistory(i, &PStat);						//Checks if has cycled around in any direction and updates

	i++;
	}

}

void DomainUpdateParticleHistory(int i, ParticleStatusData* PStat){

	Particle* P = &(Particles[i]);

	if(PStat->ParticleStatus == 0){	//Still in same domain
		P->LastVoxelX = P->VoxelX;
		P->LastVoxelY = P->VoxelY;
		P->LastVoxelZ = P->VoxelZ;
		return;
	}

	if(PStat->ParticleStatus == 1){	//Transfer into xn domain

	}

}

#endif

void RunTimeStep(double dt, int ParticleIndex0, int ParticleIndex1){

	double RandomWalkLength = sqrt(6*DiffusionCoefficient*dt);	//Computationally intensive, precalculate where possible

	int i = ParticleIndex0;
	while(i <= ParticleIndex1){

		AdvectParticle(i,dt);						//Moves particle along vector field
		RandomWalkParticle(i,dt,RandomWalkLength);	//Moves particle along a random path
		UpdateParticleHistory(i);					//Checks if has cycled around in any direction and updates

	i++;
	}

}

void RunTimeStep(double dt){

	double RandomWalkLength = sqrt(6*DiffusionCoefficient*dt);	//Computationally intensive, precalculate where possible

	int count = 0;

	int i=0;
	while(i!=NParticles){

#if DebugAdvectionMode==3
		if(count == 50000){
			cout << "Doing particle " << i << endl;
			count = 0;
		}
		count++;
#endif

		AdvectParticle(i,dt);						//Moves particle along vector field
		RandomWalkParticle(i,dt,RandomWalkLength);	//Moves particle along a random path
		UpdateParticleHistory(i);					//Checks if has cycled around in any direction and updates

	i++;
	}

}

void WriteOutParameters(){
	//cout << sizeof(Particle) << ", " << sizeof(Voxel) << ", " << sizeof(unsigned short) << endl << endl;
	
	double AverageSpeed = sqrt((LatticeInfo.AverageFlowVelocity.X*LatticeInfo.AverageFlowVelocity.X) + (LatticeInfo.AverageFlowVelocity.Y*LatticeInfo.AverageFlowVelocity.Y) + (LatticeInfo.AverageFlowVelocity.Z*LatticeInfo.AverageFlowVelocity.Z));
	double Peclet = AverageSpeed / DiffusionCoefficient;

	cout << "Porosity: " << LatticeInfo.Porosity*100 << "%" << endl;
	cout << "Average velocity = (" << LatticeInfo.AverageFlowVelocity.X << ", " << LatticeInfo.AverageFlowVelocity.Y << ", " << LatticeInfo.AverageFlowVelocity.Z << ")" << endl;
	cout << "Diffusion Coefficient Dm = " << DiffusionCoefficient << endl;
	cout << "Peclet number = " << Peclet << " x Characteristic Length" << endl;

}

//MPI Specific functions

#if Using_MPI

void InitialiseParticlesMPI(int n, int ThreadID, int ThreadNum, int* ThreadParticleCount, int ThreadParticleIndex){

	NParticles = n;							//Number of particles dealt with by thread
	Particles = new Particle[NParticles];

	InitialiseParticlesRandomXY(0,ThreadParticleIndex);

}

struct MPIDispersionStats{
	Vector SumR;
	double SumRSq;
	Vector SumRSqComponent;
};

void SendDispersionStatsToThread0(){

	MPIDispersionStats Data;

	int i=0;

	Data.SumR.X = 0;
	Data.SumR.Y = 0;
	Data.SumR.Z = 0;

	Data.SumRSq = 0;

	Data.SumRSqComponent.X = 0;
	Data.SumRSqComponent.Y = 0;
	Data.SumRSqComponent.Z = 0;

	while(i!=NParticles){

		double dX = (double)(Particles[i].XCycles*NLattice_x) + Particles[i].X;
		double dY = (double)(Particles[i].YCycles*NLattice_y) + Particles[i].Y;
		double dZ = (double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z;
		
		Data.SumR.X += dX;
		Data.SumR.Y += dY;
		Data.SumR.Z += dZ;

		Data.SumRSqComponent.X += (dX*dX);
		Data.SumRSqComponent.Y += (dY*dY);
		Data.SumRSqComponent.Z += (dZ*dZ);

		Data.SumRSq += (dX*dX + dY*dY + dZ*dZ);

		i++;
	}

	MPI_Send(&Data, sizeof(MPIDispersionStats), MPI_BYTE, 0, 0, MPI_COMM_WORLD);

}

double Thread0CalculateDispersionCoefficient(int ThreadNum, double* LastVariance, double t, double dt, int LogInterval, ofstream& outfile){

	Vector SumR;
	double SumRSq;
	Vector SumRSqComponent;

	SumR.X = 0;
	SumR.Y = 0;
	SumR.Z = 0;

	SumRSq = 0;

	SumRSqComponent.X = 0;
	SumRSqComponent.Y = 0;
	SumRSqComponent.Z = 0;

	int i=0;
	while(i!=NParticles){

		double dX = (double)(Particles[i].XCycles*NLattice_x) + Particles[i].X;
		double dY = (double)(Particles[i].YCycles*NLattice_y) + Particles[i].Y;
		double dZ = (double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z;
		
		SumR.X += dX;
		SumR.Y += dY;
		SumR.Z += dZ;

		SumRSqComponent.X += (dX*dX);
		SumRSqComponent.Y += (dY*dY);
		SumRSqComponent.Z += (dZ*dZ);

		SumRSq += (dX*dX + dY*dY + dZ*dZ);

		i++;
	}



	MPIDispersionStats Data;
	MPI_Status Stat;

	i=1;
	while(i!=ThreadNum){	//Receive and append statistics data from each thread

		MPI_Recv(&Data, sizeof(MPIDispersionStats), MPI_BYTE, i, 0, MPI_COMM_WORLD, &Stat);

		SumR.X += Data.SumR.X;
		SumR.Y += Data.SumR.Y;
		SumR.Z += Data.SumR.Z;

		SumRSq += Data.SumRSq;

		SumRSqComponent.X += Data.SumRSqComponent.X;
		SumRSqComponent.Y += Data.SumRSqComponent.Y;
		SumRSqComponent.Z += Data.SumRSqComponent.Z;

		i++;
	}

	SumR.X /= ((double)ParticleNum);	// /= total number of particles in simulation
	SumR.Y /= ((double)ParticleNum);
	SumR.Z /= ((double)ParticleNum);

	SumRSq /= ((double)ParticleNum);

	SumRSqComponent.X /= ((double)ParticleNum);
	SumRSqComponent.Y /= ((double)ParticleNum);
	SumRSqComponent.Z /= ((double)ParticleNum);

	double ExpR_Sq = ((SumR.X*SumR.X) + (SumR.Y*SumR.Y) + (SumR.Z*SumR.Z));

	double Variance = SumRSq - ExpR_Sq;

	Vector VarianceComponent = {0,0,0};
	VarianceComponent.X = SumRSqComponent.X - (SumR.X*SumR.X);
	VarianceComponent.Y = SumRSqComponent.Y - (SumR.Y*SumR.Y);
	VarianceComponent.Z = SumRSqComponent.Z - (SumR.Z*SumR.Z);

	double DispersionCoefficient = (Variance - *LastVariance)/(2*LogInterval*dt);

	printf("Calculation time: %.3lfs\n", SimulationTime.GetTimeSinceLastStep());

	cout << "Sigma^2: " << Variance << endl;
	cout << "Dispersion Coefficient: " << DispersionCoefficient << endl;

	outfile << t << '\t' << SumR.X << '\t' << VarianceComponent.X << '\t' << SumR.Y << '\t' << VarianceComponent.Y << '\t'  << SumR.Z << '\t' << VarianceComponent.Z << '\t' << Variance << '\t' << DispersionCoefficient << endl;

	*LastVariance = Variance;

	return DispersionCoefficient;
}

void SendPropagatorStatsToThread0(double RelativeBinWidth, double t){
	if(t==0){
		return;
	}

	double AverageDisplacement = t*LatticeInfo.AverageFlowVelocity.Z;
	double binwidth = RelativeBinWidth*AverageDisplacement;

	double dzmin = DBL_MAX;		//limits of graph
	double dzmax = DBL_MIN;

	int i=0;
	while(i!=NParticles){

		double dZ = ((double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z) - Particles[i].InitialZ;	//Distance particle has moved in the z direction from its starting position

		if(dZ<dzmin){
			dzmin = dZ;
		}
		if(dZ>dzmax){
			dzmax = dZ;
		}

		i++;
	}

	int lowerbinindex = (int)floor(dzmin/binwidth);		//Number of bins relative to dz=0 to the lower and upper bins
	int upperbinindex = (int)floor(dzmax/binwidth)+1;

	int nbins = upperbinindex - lowerbinindex;

	int* binvalue = new int[nbins];	//Array of ints

	i=0;
	while(i!=nbins){
		binvalue[i]=0;
		i++;
	}

	i=0;
	while(i!=NParticles){

		double dZ = ((double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z) - Particles[i].InitialZ;	//Distance particle has moved in the z direction from its starting position

		int particlebinindex = (int)floor(dZ/binwidth) - lowerbinindex;

		binvalue[particlebinindex]++;

		i++;
	}

	//Send to thread0 - lowerindex, nbins, array

	MPI_Send(&lowerbinindex, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	MPI_Send(&nbins, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);

	MPI_Send(binvalue, nbins, MPI_INT, 0, 1, MPI_COMM_WORLD);

	delete[] binvalue;

}

void Thread0OutputPropagatorGraph(double RelativeBinWidth, double t, int ThreadNum){
	if(t==0){
		return;
	}

	double AverageDisplacement = t*LatticeInfo.AverageFlowVelocity.Z;
	double binwidth = RelativeBinWidth*AverageDisplacement;

	double dzmin = DBL_MAX;		//limits of graph
	double dzmax = DBL_MIN;

	int i=0;
	while(i!=NParticles){

		double dZ = ((double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z) - Particles[i].InitialZ;	//Distance particle has moved in the z direction from its starting position

		if(dZ<dzmin){
			dzmin = dZ;
		}
		if(dZ>dzmax){
			dzmax = dZ;
		}

		i++;
	}

	int lowerbinindex = (int)floor(dzmin/binwidth);		//Number of bins relative to dz=0 to the lower and upper bins
	int upperbinindex = (int)floor(dzmax/binwidth)+1;

	int nbins = upperbinindex - lowerbinindex;

	int* binvalue = new int[nbins];	//Array of ints

	i=0;
	while(i!=nbins){
		binvalue[i]=0;
		i++;
	}

	i=0;
	while(i!=NParticles){

		double dZ = ((double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z) - Particles[i].InitialZ;	//Distance particle has moved in the z direction from its starting position

		int particlebinindex = (int)floor(dZ/binwidth) - lowerbinindex;

		binvalue[particlebinindex]++;

		i++;
	}

	//Receive and put together data from other threads

	int* ThreadLowerBinIndex = new int[ThreadNum];
	int* ThreadBinNum = new int[ThreadNum];

	ThreadLowerBinIndex[0] = lowerbinindex;
	ThreadBinNum[0] = nbins;

	MPI_Status Stat;

	i = 1;
	while(i!=ThreadNum){

		MPI_Recv(&ThreadLowerBinIndex[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, &Stat);
		MPI_Recv(&ThreadBinNum[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, &Stat);

		if(ThreadLowerBinIndex[i]<lowerbinindex){
			lowerbinindex = ThreadLowerBinIndex[i];
		}

		int ThreadUpperBinIndex = ThreadLowerBinIndex[i] + ThreadBinNum[i];

		if(ThreadUpperBinIndex>upperbinindex){
			upperbinindex = ThreadUpperBinIndex;
		}

		i++;
	}

	nbins = upperbinindex - lowerbinindex;

	int* BinValue = new int[nbins];	//Total from all threads
	i=0;
	while(i!=nbins){
		BinValue[i] = 0;
		i++;
	}

	//Thread 0 data
	int c=0;
	int threadoffset = ThreadLowerBinIndex[0] - lowerbinindex;
	while(c!=ThreadBinNum[0]){

		BinValue[threadoffset + c] = binvalue[c];

		c++;
	}

	//Other threads' data
	i=1;
	while(i!=ThreadNum){

		int* DataArray = new int[ThreadBinNum[i]];
		MPI_Recv(DataArray, ThreadBinNum[i], MPI_INT, i, 1, MPI_COMM_WORLD, &Stat);

		c=0;
		int threadoffset = ThreadLowerBinIndex[i] - lowerbinindex;
		while(c!=ThreadBinNum[i]){

			BinValue[threadoffset + c] += DataArray[c];
			c++;
		}

		delete[] DataArray;

		i++;
	}

	//Write out to file
	ofstream outfile;
	outfile.open(OutputFilePath(PropagatorOutputFile,t),ios_base::trunc);
	outfile << "Relative Displacement\tDisplacement\tNumber of Particles\tProbability" << endl;

	i=0;
	while(i!=nbins){

		double displacement = (lowerbinindex + i)*binwidth;
		double reldisplacement = displacement/(t*LatticeInfo.AverageFlowVelocity.Z);

		outfile << reldisplacement << '\t' << displacement << '\t' << BinValue[i] << '\t' << (BinValue[i]/((double)ParticleNum)) << endl;

		i++;
	}

	outfile.close();

	delete[] binvalue;
	delete[] BinValue;
	delete[] ThreadLowerBinIndex;
	delete[] ThreadBinNum;

}

#endif

//Non MPI function returns dispersion coefficient
double OutputDispersionCoefficient(double t, double dt, int LogInterval, std::ofstream& outfile, double* LastVariance){

	int i=0;
	double avgdr=0;
	double avg=0;
	Vector AverageR = {0,0,0};				//Average (vector) position
	Vector AverageComponentSq = {0,0,0};	//Average component^2 vector
	double AverageRSq = 0;					//Average position^2

	while(i!=NParticles){
		double dX = (double)(Particles[i].XCycles*NLattice_x) + Particles[i].X;
		double dY = (double)(Particles[i].YCycles*NLattice_y) + Particles[i].Y;
		double dZ = (double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z;
		
		AverageR.X += dX;
		AverageR.Y += dY;
		AverageR.Z += dZ;

		AverageComponentSq.X += (dX*dX);
		AverageComponentSq.Y += (dY*dY);
		AverageComponentSq.Z += (dZ*dZ);

		AverageRSq += (dX*dX + dY*dY + dZ*dZ);

		i++;
	}

	AverageR.X /= ((double)NParticles);
	AverageR.Y /= ((double)NParticles);
	AverageR.Z /= ((double)NParticles);

	AverageComponentSq.X /= ((double)NParticles);
	AverageComponentSq.Y /= ((double)NParticles);
	AverageComponentSq.Z /= ((double)NParticles);

	AverageRSq /= ((double)NParticles);

	double ExpR_Sq = ((AverageR.X*AverageR.X) + (AverageR.Y*AverageR.Y) + (AverageR.Z*AverageR.Z));

	double Variance = AverageRSq - ExpR_Sq;

	Vector VarianceComponent = {0,0,0};
	VarianceComponent.X = AverageComponentSq.X - (AverageR.X*AverageR.X);
	VarianceComponent.Y = AverageComponentSq.Y - (AverageR.Y*AverageR.Y);
	VarianceComponent.Z = AverageComponentSq.Z - (AverageR.Z*AverageR.Z);

	double DispersionCoefficient = (Variance - *LastVariance)/(2*LogInterval*dt);

	printf("Calculation time: %.3lfs\n", SimulationTime.GetTimeSinceLastStep());

	cout << "Sigma^2: " << Variance << endl;

	if(t!=0){
		cout << "Dispersion Coefficient: " << DispersionCoefficient << endl;
		outfile << t << '\t' << AverageR.X << '\t' << VarianceComponent.X << '\t' << AverageR.Y << '\t' << VarianceComponent.Y << '\t'  << AverageR.Z << '\t' << VarianceComponent.Z << '\t' << Variance << '\t' << DispersionCoefficient << endl;
	}else{
		outfile << t << '\t' << AverageR.X << '\t' << VarianceComponent.X << '\t' << AverageR.Y << '\t' << VarianceComponent.Y << '\t'  << AverageR.Z << '\t' << VarianceComponent.Z << '\t' << Variance << '\t' << endl;
	}

	*LastVariance = Variance;

	return DispersionCoefficient;
}

//Non-MPI function outputs propagator distribution
void OutputPropagatorGraph(double RelativeBinWidth, double t){
	if(t==0){
		return;
	}

	double AverageDisplacement = t*LatticeInfo.AverageFlowVelocity.Z;
	double binwidth = RelativeBinWidth*AverageDisplacement;

	double dzmin = DBL_MAX;		//limits of graph
	double dzmax = DBL_MIN;

	int i=0;
	while(i!=NParticles){

		double dZ = ((double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z) - Particles[i].InitialZ;	//Distance particle has moved in the z direction from its starting position

		if(dZ<dzmin){
			dzmin = dZ;
		}
		if(dZ>dzmax){
			dzmax = dZ;
		}

		i++;
	}

	int lowerbinindex = (int)floor(dzmin/binwidth);		//Number of bins relative to dz=0 to the lower and upper bins
	int upperbinindex = (int)floor(dzmax/binwidth)+1;

	int nbins = upperbinindex - lowerbinindex;

	int* binvalue = new int[nbins];	//Array of ints

	i=0;
	while(i!=nbins){
		binvalue[i]=0;
		i++;
	}

	i=0;
	while(i!=NParticles){

		double dZ = ((double)(Particles[i].ZCycles*NLattice_z) + Particles[i].Z) - Particles[i].InitialZ;	//Distance particle has moved in the z direction from its starting position

		int particlebinindex = (int)floor(dZ/binwidth) - lowerbinindex;

		binvalue[particlebinindex]++;

		i++;
	}

	ofstream outfile;
	outfile.open(OutputFilePath(PropagatorOutputFile,t),ios_base::trunc);

	outfile << "Relative Displacement\tDisplacement\tNumber of Particles\tProbability" << endl;

	i=0;
	while(i!=nbins){

		double displacement = (lowerbinindex + i)*binwidth;
		double reldisplacement = displacement/(t*LatticeInfo.AverageFlowVelocity.Z);

		outfile << reldisplacement << '\t' << displacement << '\t' << binvalue[i] << '\t' << (binvalue[i]/((double)NParticles)) << endl;

		i++;
	}

	outfile.close();
	delete[] binvalue;

}

#if Using_WinMultiThreading

struct WinThread{
	int ID;						//Thread index
	HANDLE THandle;				//Handle to the thread
	int ParticlesIndex0;		//Index of first particle under thread's auspices
	int ParticlesIndex1;		//Index of last particle under thread's auspices
	int Status;					//Running or waiting

	double t;					//Simulation time
	double tmax;				//Max simulation time
	double dt;					//Step time
	int TimeSteps;				//Number of time steps dt carried out
	int LogInterval;			//Number of timesteps between recalculating dispersion coefficient or propagator graph
	double RelativeBinWidth;	//Proportion of average flow distance for z bins in propagator distribution

};

unsigned int __stdcall ThreadMain(void* data){

	WinThread* Thread = (WinThread*)data;

	int ThreadID = Thread->ID;

	WinThread* Threads;					//Array of all threads
	if(ThreadID == 0){
		Threads = (WinThread*)data;
	}

	int ThreadNum = WinMultiThreading_TNum;
	int NParticlesThread = Thread->ParticlesIndex1 - Thread->ParticlesIndex0 + 1;

	printf("Greetings from thread %i, dealing with %i particles from index %i\n" , Thread->ID, NParticlesThread, Thread->ParticlesIndex0);

	double t = Thread->t;									//Simulation time
	double tmax = Thread->tmax;								//Max simulation time
	double dt = Thread->dt;									//Step time
	int TimeSteps = Thread->TimeSteps;						//Number of time steps dt carried out
	int LogInterval = Thread->LogInterval;					//Number of timesteps between recalculating dispersion coefficient or propagator graph
	double RelativeBinWidth = Thread->RelativeBinWidth;		//Proportion of average flow distance for z bins in propagator distribution

	int ParticlesIndex0 = Thread->ParticlesIndex0;
	int ParticlesIndex1 = Thread->ParticlesIndex1;

	double Variance = 0;			//Keeps track of last variance for dVar/dt (thread 0)

	ofstream outfile;

	if(ThreadID == 0){		//Write initial output

		outfile.open(OutputFilePath(DispersionOutputFile), ios::trunc);		//Write output file data header

		outfile << "Time\tAv X\tSigma^2 X\tAv Y\tSigma^2 Y\tAv Z\tSigma^2 Z\tSigma^2\tDispersion Coefficient" << endl;
		OutputDispersionCoefficient(t, dt, LogInterval, outfile, &Variance);

		outfile.close();

		int i=1;
		while(i!=ThreadNum){
			Threads[i].Status = 0;		//Allow threads to start
			i++;
		}	

	}else{

		while(true){		//Wait for thread 0 to output before continuing
			if(Thread->Status == 0){
				break;
			}
			Sleep(1);
		}

	}

	//Simulation
	
	while(t <= tmax){

		RunTimeStep(dt,ParticlesIndex0,ParticlesIndex1);		//Simulate thread's particles

		if(DebugAdvectionMode>=2){
			cout << "Completed timestep run [Thread " << ThreadID << "]" << endl;
		}

		TimeSteps++;
		t+=dt;

		if(TimeSteps == LogInterval){

			double dt_remainder = (double)OutputTimeInterval - (double)TimeSteps*dt;		//Carry out any remaining time
			if(dt_remainder > (10E-16)){
				RunTimeStep(dt_remainder,ParticlesIndex0,ParticlesIndex1);
				t += dt_remainder;
			}

			if(ThreadID == 0){

				while(true){	//Wait until all threads have carried out requisite number of time-steps
					bool Proceed = true;

					int i=1;
					while(i!=ThreadNum){
						if(Threads[i].Status != 1){
							Proceed = false;
						}
						i++;
					}

					if(Proceed){break;}
					Sleep(1);
				}

				outfile.open(OutputFilePath(DispersionOutputFile), ios::app);		//Write output file data header

				OutputDispersionCoefficient(t, dt, LogInterval, outfile, &Variance);

				outfile.close();

				OutputPropagatorGraph(RelativePropagatorBinWidth,t);

				int i=1;
				while(i!=ThreadNum){
					Threads[i].Status = 0;		//Allow threads to continue
					i++;
				}

			}else{

				Thread->t = t;
				Thread->TimeSteps = TimeSteps;

				Thread->Status = 1;

				while(true){		//Wait for thread 0 to output before continuing
					if(Thread->Status == 0){
						break;
					}
					Sleep(1);
				}

			}
			 
			TimeSteps = 0;

		}

	}

	return 0;
}

#endif

//Initialises lattice and particles for non MPI (used in graphics)
void Initialise(){

	CreateLattice();			//Allocate lattice and read in all data

	SetLatticeInfo();			//Fills in global LatticeInfo struct

	InitialiseParticles();		//Set up particle positions
}

double GetTimestep(){			//Calculates ideal value of dt for simulation

	double x = ( ( -sqrt(6*DiffusionCoefficient) + sqrt(6*DiffusionCoefficient + 4*LatticeInfo.MaxFlowVelocityComponent) )/(2*LatticeInfo.MaxFlowVelocityComponent) );
	double dt = x*x;

	return dt;
}

int main(int argc, char *argv[]){

	if(DebugAdvectionMode>=2){
		cout << "[Running in debug mode " << DebugAdvectionMode << "]" << endl;
		AdvectDebug.open(AdvectionDebugFile,ios::trunc);
	}

	double t = 0;						//Simulation time
	double dt = SimulationTimestep;		//Timestep
	double tmax = SimulationTimeMax;	//Simulation end time

	if(SimulationTimeMax == 0){
		tmax = DBL_MAX;					//No limit
	}
	
	int TimeSteps = 0;					//Number of time steps dt carried out between log intervals

	double Variance = 0;				//Keeps track of last variance for dVar/dt

#if !Using_MPI && !Using_WinMultiThreading	/////////////////// Non-MPI 1 CPU /////////////////////////////////////////////////

	CreateLattice();			//Allocate lattice and read in all data

	SetLatticeInfo();			//Fills in global LatticeInfo struct

	if(SimulationTimestep == 0){
		dt = GetTimestep();				//Optimum timestep
	}

	int LogInterval = (int)floor((double)OutputTimeInterval / dt);	//Number of timesteps between recalculating dispersion coefficient or propagator graph

	cout << "Simulation timestep dt = " << dt << endl;

	Random.TimeSeed();			//Seed random number generator based on local time

	InitialiseParticles();		//Initialise particles positions

	WriteOutParameters();		//Write out diffusion coefficient and average velocities to console

	//OutputParticleDistribution(1,1,1);
	
	ofstream outfile;
	outfile.open(OutputFilePath(DispersionOutputFile), ios::trunc);		//Write output file data header
	outfile << "Time\tAv X\tSigma^2 X\tAv Y\tSigma^2 Y\tAv Z\tSigma^2 Z\tSigma^2\tDispersion Coefficient" << endl;

	OutputDispersionCoefficient(t, dt, LogInterval, outfile, &Variance);

	while(t < tmax){

		RunTimeStep(dt);

		if(DebugAdvectionMode>=2){
			cout << "Completed timestep run" << endl;
			//return 0;
		}

		TimeSteps++;
		t+=dt;

		if(TimeSteps == LogInterval){

			double dt_remainder = (double)OutputTimeInterval - (double)TimeSteps*dt;		//Carry out any remaining time
			if(dt_remainder > (10E-16)){
				RunTimeStep(dt_remainder);
				t += dt_remainder;
			}

			OutputDispersionCoefficient(t, dt, LogInterval, outfile, &Variance);
			//OutputParticleDistribution(1,1,1);
			OutputPropagatorGraph(RelativePropagatorBinWidth,t);

			TimeSteps = 0;
		}

	}

	outfile.close();

#elif !Using_MPI && Using_WinMultiThreading	//////////////// Non-MPI Windows Multithreading ///////////////////////////////////

	CreateLattice();			//Allocate lattice and read in all data

	SetLatticeInfo();			//Fills in global LatticeInfo struct

	if(SimulationTimestep == 0){
		dt = GetTimestep();				//Optimum timestep
	}

	int LogInterval = (int)floor((double)OutputTimeInterval / dt);	//Number of timesteps between recalculating dispersion coefficient or propagator graph

	cout << "Simulation timestep dt = " << dt << endl;

	Random.TimeSeed();			//Seed random number generator based on local time

	InitialiseParticles();		//Initialise particles positions

	WriteOutParameters();		//Write out diffusion coefficient and average velocities to console

	int ThreadNum = WinMultiThreading_TNum;

	WinThread* Threads = new WinThread[ThreadNum];

	int* ThreadParticleCount = new int[ThreadNum];

	int ModParticles = NParticles%ThreadNum;						//Remainder of particles after dividing up
	int NParticlesThread = (NParticles - ModParticles)/ThreadNum;	//Divide up particles among processors

	int i=0;
	while(i!=ThreadNum){
		ThreadParticleCount[i] = NParticlesThread;
		i++;
	}

	while(ModParticles>0){
		i--;

		ThreadParticleCount[i]++;	//Divide up remainder among threads
		ModParticles--;

		if(i==0){
			i = ThreadNum;
		}
	}

	Threads[0].ID = 0;
	Threads[0].ParticlesIndex0 = 0;
	Threads[0].ParticlesIndex1 = ThreadParticleCount[0] - 1;

	Threads[0].t = t;
	Threads[0].dt = dt;
	Threads[0].tmax = tmax;
	Threads[0].LogInterval = LogInterval;
	Threads[0].RelativeBinWidth = RelativePropagatorBinWidth;
	Threads[0].TimeSteps = TimeSteps;

	i=1;
	while(i!=ThreadNum){
		
		Threads[i].ID = i;
		Threads[i].ParticlesIndex0 = Threads[i-1].ParticlesIndex1 + 1;
		Threads[i].ParticlesIndex1 = Threads[i].ParticlesIndex0 + ThreadParticleCount[i] - 1;

		Threads[i].t = t;
		Threads[i].dt = dt;
		Threads[i].tmax = tmax;
		Threads[i].LogInterval = LogInterval;
		Threads[i].RelativeBinWidth = RelativePropagatorBinWidth;
		Threads[i].TimeSteps = TimeSteps;

		Threads[i].Status = 1;	//Paused state
		
		_beginthreadex(0, 0, ThreadMain, &Threads[i], 0, 0);	//Threads

		i++;
	}

	ThreadMain((void*)(&Threads[0]));	//Thread 0

	delete[] ThreadParticleCount;
	delete[] Threads;


#elif (MPI_Mode!=1)	//////////////////////////////////////// MPI Full Lattice /////////////////////////////////////////////////

	CreateLattice();			//Allocate lattice and read in all data

	SetLatticeInfo();			//Fills in global LatticeInfo struct

	if(SimulationTimestep == 0){
		dt = GetTimestep();				//Optimum timestep
	}

	int LogInterval = (int)floor((double)OutputTimeInterval / dt);	//Number of timesteps between recalculating dispersion coefficient or propagator graph

	cout << "Simulation timestep dt = " << dt << endl;

	int ThreadID, ThreadNum;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ThreadNum);	//Get number of threads running
	MPI_Comm_rank(MPI_COMM_WORLD, &ThreadID);	//Get thread ID (index from 0 -> ThreadNum-1)

	//Multithread zone

	int* ThreadParticleCount = new int[ThreadNum];					//Number of particles for each thread

	int ModParticles = ParticleNum%ThreadNum;						//Remainder of particles after dividing up
	int NParticlesThread = (ParticleNum - ModParticles)/ThreadNum;	//Divide up particles among processors

	int i=0;
	while(i!=ThreadNum){
		ThreadParticleCount[i] = NParticlesThread;
		i++;
	}

	while(ModParticles>0){
		i--;

		ThreadParticleCount[i]++;	//Divide up remainder among threads
		ModParticles--;

		if(i==0){
			i = ThreadNum;
		}
	}

	NParticlesThread = ThreadParticleCount[ThreadID];

	int ThreadParticleIndex = 0;	//Gives the index of the first particle this thread deals with

	i=0;
	while(i!=ThreadID){
		ThreadParticleIndex += ThreadParticleCount[i];
		i++;
	}

	Random.ThreadSeed(ThreadID,ThreadNum);		//Random number seed depending on thread and local time

	InitialiseParticlesMPI(NParticlesThread, ThreadID, ThreadNum, ThreadParticleCount, ThreadParticleIndex);

	ofstream outfile;	//FileIO dealt with by Thread 0

	if(ThreadID==0){
		WriteOutParameters();	//Write out diffusion coefficient and average velocities to console
		
		outfile.open(OutputFilePath(DispersionOutputFile), ios::trunc);		//Write output file data header
		outfile << "Time\tAv X\tSigma^2 X\tAv Y\tSigma^2 Y\tAv Z\tSigma^2 Z\tSigma^2\tDispersion Coefficient" << endl;
	}

	printf("Greetings from thread %i dealing with %i particles starting at index %i\n", ThreadID, NParticlesThread, ThreadParticleIndex);

	//Main thread loop

	while(t<=tmax){

		RunTimeStep(dt);

		TimeSteps++;
		t+=dt;

		if(TimeSteps == LogInterval){

			double dt_remainder = (double)OutputTimeInterval - (double)TimeSteps*dt;		//Carry out any remaining time
			if(dt_remainder > (10E-16)){
				RunTimeStep(dt_remainder);
				t += dt_remainder;
			}

			if(ThreadID==0){
				Thread0CalculateDispersionCoefficient(ThreadNum,&Variance,t,dt,LogInterval,outfile);
				Thread0OutputPropagatorGraph(RelativePropagatorBinWidth,t,ThreadNum);
			}else{
				SendDispersionStatsToThread0();
				SendPropagatorStatsToThread0(RelativePropagatorBinWidth,t);
			}

			TimeSteps = 0;
		}

	}

MPIFinalisation:			//Jump to here from inside loop to finalise and finish

	if(ThreadID==0){
		outfile.close();
	}
	delete[] Particles;
	MPI_Finalize();		//End of multithreading

#else	//////////////////////////////////////// MPI Lattice Division /////////////////////////////////////////////////

	int ThreadID, ThreadNum;
	MPI_Status Stat;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ThreadNum);	//Get number of threads running
	MPI_Comm_rank(MPI_COMM_WORLD, &ThreadID);	//Get thread ID (index from 0 -> ThreadNum-1)

	//Multithread zone

	Random.ThreadSeed(ThreadID,ThreadNum);		//Random number seed depending on thread and local time

	ofstream outfile;	//FileIO dealt with by Thread 0

	int NParticlesThread;	//Number of particles dealt with by thread;

	if(ThreadID==0){

		if(LatticeDecomposition(ThreadNum, ThreadID)!=0){		//Decomposes geometry and sends all initialisation data to threads
			goto MPIFinalisation;	//Decomposition failed
		}	

		outfile.open(OutputFilePath(DispersionOutputFile), ios::trunc);		//Write output file data header
		outfile << "Time\tAv X\tSigma^2 X\tAv Y\tSigma^2 Y\tAv Z\tSigma^2 Z\tSigma^2\tDispersion Coefficient" << endl;
		outfile.close();

		//Initialise particles
		DomainInitialiseParticlesUniform(ThreadNum, ThreadID, &NParticlesThread);

		//Lattice information
		WriteOutParameters();

		//Wait for threads to signal ready
		int Flag = 0;
		int c = 1;
		while(c!=ThreadNum){
			MPI_Recv(&Flag, 1, MPI_INT, c, 0, MPI_COMM_WORLD, &Stat);
			c++;
		}
		//Inform threads to begin simulation
		c = 1;
		while(c!=ThreadNum){
			MPI_Send(&Flag, 1, MPI_INT, c, 0, MPI_COMM_WORLD);
			c++;
		}

	}else{

		if(ReceiveLatticeDecompositionData(ThreadNum, ThreadID)!=0){
			goto MPIFinalisation;	//Decomposition failed
		}

		//Initialise particles
		DomainInitialiseParticlesUniform(ThreadNum, ThreadID, &NParticlesThread);

		//Wait to receive begin flag from thread 0
		int Flag = 0;
		MPI_Send(&Flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);				//Send ready flag
		MPI_Recv(&Flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);		//Recieve start simulation flag

	}

	//

	if(SimulationTimestep == 0){
		dt = GetTimestep();				//Optimum timestep
	}

	int LogInterval = (int)floor((double)OutputTimeInterval / dt);	//Number of timesteps between recalculating dispersion coefficient or propagator graph

	//Simulation

	while(t<=tmax){

		DomainRunTimeStep(dt);

		TimeSteps++;
		t+=dt;

		if(TimeSteps == LogInterval){

			double dt_remainder = (double)OutputTimeInterval - (double)TimeSteps*dt;		//Carry out any remaining time
			if(dt_remainder > (10E-16)){
				DomainRunTimeStep(dt_remainder);
				t += dt_remainder;
			}

			if(ThreadID==0){
				
			}else{

			}

			TimeSteps = 0;
		}

	}


MPIFinalisation:			//Jump to here from inside loop to finalise and finish

	delete[] LatticeBase;
	delete[] Lattice;
	delete[] Particles;
	MPI_Finalize();		//End of multithreading

#endif

	return 0;
}

void UpdateParticleHistory(int i){

	Particle* P = &(Particles[i]);
	
	const short LatticeQX = NLattice_x / 4;
	const short LatticeQY = NLattice_y / 4;
	const short LatticeQZ = NLattice_z / 4;
	
	short dVX = P->VoxelX - P->LastVoxelX;
	if(dVX > LatticeQX){
		P->XCycles--;
	}else if(dVX < -LatticeQX){
		P->XCycles++;
	}

	short dVY = P->VoxelY - P->LastVoxelY;
	if(dVY > LatticeQY){
		P->YCycles--;
	}else if(dVY < -LatticeQY){
		P->YCycles++;
	}

	short dVZ = P->VoxelZ - P->LastVoxelZ;
	if(dVZ > LatticeQZ){
		P->ZCycles--;
	}else if(dVZ < -LatticeQZ){
		P->ZCycles++;
	}

	P->LastVoxelX = P->VoxelX;
	P->LastVoxelY = P->VoxelY;
	P->LastVoxelZ = P->VoxelZ;

}

void GetCoordinates(short x, short y, short z, Coords* C){

	C->x=x;
	C->y=y;
	C->z=z;
	
	C->xp=C->x+1;				//Voxel to the right
	C->xn=C->x-1;				//Voxel to the left
	C->yp=C->y+1;				
	C->yn=C->y-1;					
	C->zp=C->z+1;				
	C->zn=C->z-1;				

	if(C->xp==NLattice_x){C->xp=0;}	//If voxel out of bounds, apply periodic boundary conditions
	if(C->xn<0){C->xn=NLattice_x-1;}
	if(C->yp==NLattice_y){C->yp=0;}	
	if(C->yn<0){C->yn=NLattice_y-1;}
	if(C->zp==NLattice_z){C->zp=0;}
	if(C->zn<0){C->zn=NLattice_z-1;}

}

void GetCoordinates(int i, Coords* C){

	C->x = Particles[i].VoxelX;	//voxel that particle is in
	C->y = Particles[i].VoxelY;
	C->z = Particles[i].VoxelZ;
	
	C->xp = C->x+1;				//Voxel to the right
	C->xn = C->x-1;				//Voxel to the left
	C->yp = C->y+1;				
	C->yn = C->y-1;					
	C->zp = C->z+1;				
	C->zn = C->z-1;				

	if(C->xp==NLattice_x){C->xp=0;}	//If voxel out of bounds, apply periodic boundary conditions
	if(C->xn<0){C->xn=NLattice_x-1;}
	if(C->yp==NLattice_y){C->yp=0;}	
	if(C->yn<0){C->yn=NLattice_y-1;}
	if(C->zp==NLattice_z){C->zp=0;}
	if(C->zn<0){C->zn=NLattice_z-1;}

}

void DomainGetCoordinates(short x, short y, short z, Coords* C){

	C->x = x;
	C->y = y;
	C->z = z;
	
	C->xp = C->x + 1;				//Voxel to the right
	C->xn = C->x - 1;				//Voxel to the left
	C->yp = C->y + 1;				
	C->yn = C->y - 1;					
	C->zp = C->z + 1;				
	C->zn = C->z - 1;				
	
	//Allow C->x - 1 to be negative
}

#if Using_MPI	//MPI Functions

void DomainAdvectParticle(int i, double dt, int* ParticleStatus){

	Particle* P = &(Particles[i]);

	int cyclecount = 0;

	while(dt>0.0){

	cyclecount++;

	short x = P->VoxelX;		//Domain X, Y, Z
	short y = P->VoxelY;
	short z = P->VoxelZ;

	Coords C;
	DomainGetCoordinates(x,y,z,&C);		//Gets coordinates of voxel and its neighbours

	Voxel* LatticeElement = GetLattice(x,y,z);

	short NSolids = 0;

	int c=0;
	unsigned char Solids = (unsigned char)(LatticeElement->Solid);
	while(c!=6){
		if((Solids>>=1)&0x01){			//Count number of neighbour solid flag bits set
			NSolids++;
		}
		c++;
	}

	switch(NSolids){

		case 0:
			CaseNoSolids(i,&C,LatticeElement, &dt);
			break;

		case 1:
			CaseOneSolid(i,&C,LatticeElement, &dt);
			break;

		case 2:
			CaseTwoSolids(i,&C,LatticeElement, &dt);
			break;
			
		case 3:
			CaseThreeSolids(i,&C,LatticeElement, &dt);
			break;

		case 4:
			CaseFourSolids(i,&C,LatticeElement, &dt);
			break;

		case 5:
		case 6:
			dt = 0;
			break;

		}

	//Check if particle moved cross-domains
	x = P->VoxelX;
	y = P->VoxelY;
	z = P->VoxelZ;

	if(x==-1){									//Transfer into -ve x domain
		*ParticleStatus = 1;
		return;
	}else if(x==ProcessorDomain.DomainWidthX){	//Transfer into +ve x domain
		*ParticleStatus = 2;
		return;
	}

	if(y==-1){									//Transfer into -ve y domain
		*ParticleStatus = 3;
		return;
	}else if(y==ProcessorDomain.DomainWidthY){	//Transfer into +ve y domain
		*ParticleStatus = 4;
		return;
	}

	if(z==-1){									//Transfer into -ve z domain
		*ParticleStatus = 5;
		return;
	}else if(z==ProcessorDomain.DomainWidthZ){	//Transfer into +ve z domain
		*ParticleStatus = 6;
		return;
	}

	if(cyclecount==8){
		cout << "Particle " << i << " with " << NSolids << " solids exceeded 8 cycles, terminating." << endl;
		cout << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;
		dt = 0;
		break;
	}
	
	}

	*ParticleStatus = 0;
}

#endif

void AdvectParticle(int i,double dt){

	Particle* P = &(Particles[i]);

#if DebugAdvectionMode>=2
	double Restoredt;
	Particle RestoreParticle;
	
	Restoredt = dt;
	memcpy(&RestoreParticle,P,sizeof(Particle));

AdvectDebugRestart:			//Jump here to restart advection with newly restored particle in debug mode
#endif

	int cyclecount = 0;

	while(dt>0.0){

	cyclecount++;

	short x = P->VoxelX;
	short y = P->VoxelY;
	short z = P->VoxelZ;

	Coords C;
	GetCoordinates(x,y,z,&C);		//Gets coordinates of voxel and its neighbours

	Voxel* LatticeElement = GetLattice(x,y,z);

	short NSolids = 0;

	int c=0;
	unsigned char Solids = (unsigned char)(LatticeElement->Solid);
	while(c!=6){
		if((Solids>>=1)&0x01){			//Count number of neighbour solid flag bits set
			NSolids++;
		}
		c++;
	}

#if DebugAdvectionMode>=2
	if(DebugParticle){
		AdvectDebug << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ") with " << NSolids << " solids" << endl;
		AdvectDebug << "Voxel Information:" << endl << "(Vx, Vy, Vz) = " << "(" << LatticeElement->Vx << ", " << LatticeElement->Vy << ", " << LatticeElement->Vz << ")" << endl;
		AdvectDebug << "u1 = " << LatticeElement->u1 << endl << "u2 = " << LatticeElement->u2 << endl << "v1 = " << LatticeElement->v1 << endl << "v2 = " << LatticeElement->v2 << endl << "w1 = " << LatticeElement->w1 << endl << "w2 = " << LatticeElement->w2 << endl;
	}
#endif

	switch(NSolids){

		case 0:
			CaseNoSolids(i,&C,LatticeElement, &dt);
			break;

		case 1:
			CaseOneSolid(i,&C,LatticeElement, &dt);
			break;

		case 2:
			CaseTwoSolids(i,&C,LatticeElement, &dt);
			break;
			
		case 3:
			CaseThreeSolids(i,&C,LatticeElement, &dt);
			break;

		case 4:
			CaseFourSolids(i,&C,LatticeElement, &dt);
			break;

		case 5:
		case 6:
			dt = 0;
			break;

		}

#if DebugAdvectionMode==1

	if(cyclecount==8){
		cout << "Particle " << i << " with " << NSolids << " solids exceeded 8 cycles, terminating." << endl;
		cout << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;
		AdvectDebug << "Particle " << i << " with " << NSolids << " solids exceeded 8 cycles, terminating." << endl;
		AdvectDebug << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;

		dt = 0;
		break;
	}

#endif

#if DebugAdvectionMode>=2
	///////////////////////////// Advection Debugging /////////////////////////////////////////////////////////////////////////////////////////////////

	if(cyclecount == 8){
		if(!DebugParticle){
			cout << "Trouble with particle " << i << " with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Trouble with particle " << i << " with " << NSolids << " solids. Debugging:" << endl;

			dt = Restoredt;
			memcpy(P,&RestoreParticle,sizeof(Particle));
			DebugParticle = true;
		
			goto AdvectDebugRestart;
		}else{
			AdvectDebug << "Exceeded 8 cycles. Terminating." << endl << endl;
			DebugParticle = false;
			memcpy(P,&RestoreParticle,sizeof(Particle));	//Restore for randomwalk
			dt = 0;
			break;
		}
	}

	if(DebugAdvectionMode==3&&!DebugParticle){	//Very slow, catches #inf errors

	char DebugChar[128];
	int _n = _snprintf(DebugChar,128,"%f%f%f",P->X,P->Y,P->Z);
	int _i=0;
	while(_i!=_n){
		if(DebugChar[_i]=='#'){
			cout  << "Indef. fault after advection: Particle " << i << " with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Indef. fault after advection: Particle " << i << " with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Error in position: (" << P->X << ", " << P->Y << ", " << P->Z << ")" << endl;

			dt = Restoredt;
			memcpy(P,&RestoreParticle,sizeof(Particle));
			DebugParticle = true;
		
			goto AdvectDebugRestart;
		}
		_i++;
	}

	}

	if((P->VoxelX != (int)floor(P->X) && P->VoxelX+1 != (int)floor(P->X)) || (P->VoxelY != (int)floor(P->Y) && P->VoxelY+1 != (int)floor(P->Y)) || (P->VoxelZ != (int)floor(P->Z) && P->VoxelZ+1 != (int)floor(P->Z))){
		if(!DebugParticle){
			cout << "Numerical fault after advection: Particle " << i << " with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Numerical fault after advection: Particle " << i << " with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Error in position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;

			dt = Restoredt;
			memcpy(P,&RestoreParticle,sizeof(Particle));
			DebugParticle = true;
		
			goto AdvectDebugRestart;
		}else{
			AdvectDebug << endl << endl;
			DebugParticle = false;
			memcpy(P,&RestoreParticle,sizeof(Particle));	//Restore for randomwalk
			dt = 0;
			break;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

	}

#if DebugAdvection
	DebugParticle = false;
#endif
	
}

void CycleCoordinates(Vector* V){	//Applies boundary conditions to coordinates outside of lattice range

	const double LSizeX = (double)NLattice_x;
	const double LSizeY = (double)NLattice_y;
	const double LSizeZ = (double)NLattice_z;

	if(V->X<0){
		V->X=LSizeX + V->X;	//If below 0, apply periodic boundary condition
	}else if(V->X>=LSizeX){
		V->X-=LSizeX;
	}

	if(V->Y<0){
		V->Y=LSizeY + V->Y;	//If below 0, apply periodic boundary condition
	}else if(V->Y>=LSizeY){
		V->Y-=LSizeY;
	}

	if(V->Z<0){
		V->Z=LSizeZ + V->Z;	//If below 0, apply periodic boundary condition
	}else if(V->Z>=LSizeZ){
		V->Z-=LSizeZ;
	}

}

void CycleCoordinates(double* Px, double* Py, double* Pz){	//Applies boundary conditions to coordinates outside of lattice range

	const double LSizeX = (double)NLattice_x;
	const double LSizeY = (double)NLattice_y;
	const double LSizeZ = (double)NLattice_z;

	if((*Px)<0){
		(*Px)=LSizeX + (*Px);	//If below 0, apply periodic boundary condition
	}else if((*Px)>=LSizeX){
		(*Px)-=LSizeX;
	}

	if((*Py)<0){
		(*Py)=LSizeY + (*Py);	//If below 0, apply periodic boundary condition
	}else if((*Py)>=LSizeY){
		(*Py)-=LSizeY;
	}

	if((*Pz)<0){
		(*Pz)=LSizeZ + (*Pz);	//If below 0, apply periodic boundary condition
	}else if((*Pz)>=LSizeZ){
		(*Pz)-=LSizeZ;
	}

}

void CycleVoxelCoordinates(int* Vx, int* Vy, int* Vz){	//Applies boundary conditions to voxel indices out of lattice range

	if((*Vx)<0){
		(*Vx)+=NLattice_x;
	}else if((*Vx)>=NLattice_x){
		(*Vx)-=NLattice_x;
	}

	if((*Vy)<0){
		(*Vy)+=NLattice_y;
	}else if((*Vy)>=NLattice_y){
		(*Vy)-=NLattice_y;
	}

	if((*Vz)<0){
		(*Vz)+=NLattice_z;
	}else if((*Vz)>=NLattice_z){
		(*Vz)-=NLattice_z;
	}

}

void VoxelAfterDisplacement(Particle* P, Vector* V,Coords* C){
	
	double Px = P->X + V->X;
	double Py = P->Y + V->Y;
	double Pz = P->Z + V->Z;

	double dcx = V->X + (P->X - (double)P->VoxelX);
	int cx = (int)floor(dcx);

	double dcy = V->Y + (P->Y - (double)P->VoxelY);
	int cy = (int)floor(dcy);

	double dcz = V->Z + (P->Z - (double)P->VoxelZ);
	int cz = (int)floor(dcz);

	short x = P->VoxelX;
	short y = P->VoxelY;
	short z = P->VoxelZ;

	int Vx=0;
	int Vy=0;
	int Vz=0;

	GetCoordinates(x,y,z,C);

	while(Vx!=cx){

		if(cx<0){
			GetCoordinates(C->xn,C->y,C->z,C);
		}else{
			GetCoordinates(C->xp,C->y,C->z,C);
		}

		if(cx<0){
			Vx--;
		}else{
			Vx++;
		}

	}

	while(Vy!=cy){

		if(cy<0){
			GetCoordinates(C->x,C->yn,C->z,C);
		}else{
			GetCoordinates(C->x,C->yp,C->z,C);
		}


		if(cy<0){
			Vy--;
		}else{
			Vy++;
		}

	}

	while(Vz!=cz){

		if(cz<0){
			GetCoordinates(C->x,C->y,C->zn,C);
		}else{
			GetCoordinates(C->x,C->y,C->zp,C);
		}


		if(cz<0){
			Vz--;
		}else{
			Vz++;
		}

	}

}

void DomainVoxelAfterDisplacement(Particle* P, Vector* V,Coords* C){
	
	double Px = P->X + V->X;
	double Py = P->Y + V->Y;
	double Pz = P->Z + V->Z;

	double dcx = V->X + (P->X - (double)P->VoxelX);
	int cx = (int)floor(dcx);

	double dcy = V->Y + (P->Y - (double)P->VoxelY);
	int cy = (int)floor(dcy);

	double dcz = V->Z + (P->Z - (double)P->VoxelZ);
	int cz = (int)floor(dcz);

	short x = P->VoxelX;
	short y = P->VoxelY;
	short z = P->VoxelZ;
	
	DomainGetCoordinates(x+cx, y+cy, z+cz, C);

}

void RandomWalkVector(double L, Vector* V){		//Returns a random walk vector of length L

	double zp = 2.0*Random.GetUniform() - 1.0;

	double phi = acos(zp);
	double theta = 2*PI*Random.GetUniform();

	double LSinPhi = L*sin(phi);

	V->X = LSinPhi*cos(theta);
	V->Y = LSinPhi*sin(theta);
	V->Z = L*zp;
}

void RandomWalkParticle(int i, double dt, double RandomWalkLength){

	Vector V;

	RandomWalkVector(RandomWalkLength,&V); //Get spherically isotropic random vector

	Particle* P = &(Particles[i]);

	while(true){	//Loops until the walk vector has been carried out entirely

	int CurrentVoxelX = P->VoxelX;
	int CurrentVoxelY = P->VoxelY;
	int CurrentVoxelZ = P->VoxelZ;

	Voxel* LatticeElement = GetLattice(CurrentVoxelX, CurrentVoxelY, CurrentVoxelZ);	//Voxel particle is currently in

	Coords C;					
	GetCoordinates(CurrentVoxelX,CurrentVoxelY,CurrentVoxelZ,&C);	//Gets coordinates of voxel and its neighbours
	
	double Px = P->X + V.X;
	double Py = P->Y + V.Y;
	double Pz = P->Z + V.Z;

	Coords Cafter;
	VoxelAfterDisplacement(P,&V,&Cafter);

	int NewX = Cafter.x;	//Voxel that particle is in after random walk
	int NewY = Cafter.y;
	int NewZ = Cafter.z;

	if(NewX == CurrentVoxelX && NewY == CurrentVoxelY && NewZ == CurrentVoxelZ){	//Still in same voxel, no solid boundaries to worry about

		CycleCoordinates(&Px,&Py,&Pz);	//Apply boundary conditions to coordinates

		P->X = Px;		//Update particle position
		P->Y = Py;
		P->Z = Pz;

		return;		//Exits loop here
	}

	//Particle has walked into a neighbouring voxel

	double Xdist,Ydist,Zdist;	//distance to nearest side of voxel in each direction

	if(V.X<0){
		Xdist = P->X - CurrentVoxelX;		//distance to -x side from original particle position
	}else{
		Xdist = (CurrentVoxelX + 1) - P->X;	//distance to +x side
	}

	if(V.Y<0){
		Ydist = P->Y - CurrentVoxelY;		//distance to -y side from original particle position
	}else{
		Ydist = (CurrentVoxelY + 1) - P->Y;	//distance to +y side
	}

	if(V.Z<0){
		Zdist = P->Z - CurrentVoxelZ;		//distance to -z side from original particle position
	}else{
		Zdist = (CurrentVoxelZ + 1) - P->Z;	//distance to +z side
	}

	double Xtime = Xdist/fabs(V.X);	//Relative time for component to reach nearest voxel boundary
	double Ytime = Ydist/fabs(V.Y);	//	and ratios of (distance to travel) / (component length),
	double Ztime = Zdist/fabs(V.Z);	//	ie proportion of walk which would be carried out by moving to border

	//Which component x,y or z reaches the boundary first

	if(Ztime<=Ytime && Ztime<=Xtime){		//Reaches z boundary first

		bool VSolid;					//Is bounding voxel solid

		if((V.Z<0 && (LatticeElement->Solid&NegativeZSolid)) || (V.Z>=0 && (LatticeElement->Solid&PositiveZSolid))){
			VSolid=true;
		}else{
			VSolid=false;
		}

		if(V.Z<0 && !VSolid){					//if moving in -z and neighbour is not solid

			P->Z = (double)(C.zn+1);				//Move particle to voxel boundary and inside lower voxel
			P->VoxelZ = C.zn;

		}else if(V.Z<0 && VSolid){				//if moving in -z and neighbour is solid

			P->Z = (double)CurrentVoxelZ;				//Move particle to border but inside current voxel

		}else if(V.Z>0 && !VSolid){				//if moving in +z and neighbour is not solid

			P->Z = (double)C.zp;			//Move particle to voxel boundary and into upper voxel
			P->VoxelZ = C.zp;

		}else{									//if moving in +z and neighbour is solid

			P->Z = (double)(CurrentVoxelZ+1);			//Move particle to voxel boundary but inside current voxel

		}

			P->X += V.X*Ztime;	//Move in proportion to Vz moved
			P->Y += V.Y*Ztime;

			V.X *= (1-Ztime);				//Reduce walk vector by amount moved
			V.Y *= (1-Ztime);
			V.Z *= (1-Ztime);
		
		if(VSolid){
			V.Z = -V.Z;		//If boundary is solid, reverse component (reflection)
		}
					
	}else if(Ytime<=Xtime && Ytime<=Ztime){		//Reaches y boundary first

		bool VSolid;					//Is bounding voxel solid

		if((V.Y<0 && (LatticeElement->Solid&NegativeYSolid)) || (V.Y>=0 && (LatticeElement->Solid&PositiveYSolid))){
			VSolid=true;
		}else{
			VSolid=false;
		}

		if(V.Y<0 && !VSolid){					//if moving in -y and neighbour is not solid

			P->Y = (double)(C.yn+1);				//Move particle to voxel boundary and inside left voxel
			P->VoxelY = C.yn;

		}else if(V.Y<0 && VSolid){				//if moving in -y and neighbour is solid

			P->Y = (double)CurrentVoxelY;				//Move particle to border but inside current voxel

		}else if(V.Y>0 && !VSolid){				//if moving in +y and neighbour is not solid

			P->Y = (double)C.yp;			//Move particle to voxel boundary and into right voxel
			P->VoxelY = C.yp;

		}else{									//if moving in +y and neighbour is solid

			P->Y = (double)(CurrentVoxelY+1);					//Move particle to voxel boundary but inside current voxel

		}

			P->X += V.X*Ytime;				//Move in proportion to Vy moved
			P->Z += V.Z*Ytime;

			V.X *= (1-Ytime);				//Reduce walk vector by amount moved
			V.Y *= (1-Ytime);
			V.Z *= (1-Ytime);
		
		if(VSolid){
			V.Y = -V.Y;					//If boundary is solid, reverse component (reflection)
		}

	}else{								//Reaches x boundary first

		bool VSolid;					//Is bounding voxel solid

		if((V.X<0 && (LatticeElement->Solid&NegativeXSolid)) || (V.X>=0 && (LatticeElement->Solid&PositiveXSolid))){
			VSolid=true;
		}else{
			VSolid=false;
		}

		if(V.X<0 && !VSolid){					//if moving in -x and neighbour is not solid

			P->X = (double)(C.xn+1);				//Move particle to voxel boundary and inside left voxel
			P->VoxelX = C.xn;

		}else if(V.X<0 && VSolid){				//if moving in -x and neighbour is solid

			P->X = (double)CurrentVoxelX;			//Move particle to border but inside current voxel

		}else if(V.X>0 && !VSolid){				//if moving in +x and neighbour is not solid

			P->X = (double)C.xp;					//Move particle to voxel boundary and into right voxel
			P->VoxelX = C.xp;

		}else{									//if moving in +x and neighbour is solid

			P->X = (double)(CurrentVoxelX+1);		//Move particle to voxel boundary but inside current voxel

		}

			P->Y += V.Y*Xtime;	//Move in proportion to Vx moved
			P->Z += V.Z*Xtime;

			V.X *= (1-Xtime);				//Reduce walk vector by amount moved
			V.Y *= (1-Xtime);
			V.Z *= (1-Xtime);
		
		if(VSolid){
			V.X = -V.X;	//If boundary is solid, reverse component (reflection)
		}

	}

	}	//end of while loop

}

#if Using_MPI	//MPI Functions

void DomainRandomWalkParticle(int i, double dt, double RandomWalkLength, int* ParticleStatus){

	Vector V;

	RandomWalkVector(RandomWalkLength,&V); //Get spherically isotropic random vector

	Particle* P = &(Particles[i]);

	while(true){	//Loops until the walk vector has been carried out entirely

	short CurrentVoxelX = P->VoxelX;
	short CurrentVoxelY = P->VoxelY;
	short CurrentVoxelZ = P->VoxelZ;

	Voxel* LatticeElement = GetLattice(CurrentVoxelX, CurrentVoxelY, CurrentVoxelZ);	//Voxel particle is currently in

	Coords C;					
	DomainGetCoordinates(CurrentVoxelX,CurrentVoxelY,CurrentVoxelZ,&C);	//Gets coordinates of voxel and its neighbours
	
	double Px = P->X + V.X;
	double Py = P->Y + V.Y;
	double Pz = P->Z + V.Z;
	
	Coords Cafter;
	DomainVoxelAfterDisplacement(P,&V,&Cafter);

	int NewX = Cafter.x;	//Voxel that particle is in after random walk
	int NewY = Cafter.y;
	int NewZ = Cafter.z;

	if(NewX == CurrentVoxelX && NewY == CurrentVoxelY && NewZ == CurrentVoxelZ){	//Still in same voxel, no solid boundaries to worry about
		
		P->X = Px;		//Update particle position
		P->Y = Py;
		P->Z = Pz;

		return;		//Exits loop here
	}

	//Particle has walked into a neighbouring voxel

	double Xdist,Ydist,Zdist;	//distance to nearest side of voxel in each direction

	if(V.X<0){
		Xdist = P->X - (double)CurrentVoxelX;		//distance to -x side from original particle position
	}else{
		Xdist = (double)(CurrentVoxelX + 1) - P->X;	//distance to +x side
	}

	if(V.Y<0){
		Ydist = P->Y - (double)CurrentVoxelY;		//distance to -y side from original particle position
	}else{
		Ydist = (double)(CurrentVoxelY + 1) - P->Y;	//distance to +y side
	}

	if(V.Z<0){
		Zdist = P->Z - (double)CurrentVoxelZ;		//distance to -z side from original particle position
	}else{
		Zdist = (double)(CurrentVoxelZ + 1) - P->Z;	//distance to +z side
	}

	double Xtime = Xdist/fabs(V.X);	//Relative time for component to reach nearest voxel boundary
	double Ytime = Ydist/fabs(V.Y);	//	and ratios of (distance to travel) / (component length),
	double Ztime = Zdist/fabs(V.Z);	//	ie proportion of walk which would be carried out by moving to border

	//Which component x,y or z reaches the boundary first

	if(Ztime<=Ytime && Ztime<=Xtime){		//Reaches z boundary first

		bool VSolid;					//Is bounding voxel solid

		if((V.Z<0 && (LatticeElement->Solid&NegativeZSolid)) || (V.Z>=0 && (LatticeElement->Solid&PositiveZSolid))){
			VSolid=true;
		}else{
			VSolid=false;
		}

		if(V.Z<0 && !VSolid){					//if moving in -z and neighbour is not solid

			P->Z = (double)(C.zn+1);			//Move particle to voxel boundary and inside lower voxel
			P->VoxelZ = C.zn;

		}else if(V.Z<0 && VSolid){				//if moving in -z and neighbour is solid

			P->Z = (double)CurrentVoxelZ;		//Move particle to border but inside current voxel

		}else if(V.Z>0 && !VSolid){				//if moving in +z and neighbour is not solid

			P->Z = (double)C.zp;				//Move particle to voxel boundary and into upper voxel
			P->VoxelZ = C.zp;

		}else{									//if moving in +z and neighbour is solid

			P->Z = (double)(CurrentVoxelZ+1);	//Move particle to voxel boundary but inside current voxel

		}

			P->X += V.X*Ztime;				//Move in proportion to Vz moved
			P->Y += V.Y*Ztime;

			V.X *= (1-Ztime);				//Reduce walk vector by amount moved
			V.Y *= (1-Ztime);
			V.Z *= (1-Ztime);
		
		if(VSolid){
			V.Z = -V.Z;		//If boundary is solid, reverse component (reflection)
		}
					
	}else if(Ytime<=Xtime && Ytime<=Ztime){		//Reaches y boundary first

		bool VSolid;					//Is bounding voxel solid

		if((V.Y<0 && (LatticeElement->Solid&NegativeYSolid)) || (V.Y>=0 && (LatticeElement->Solid&PositiveYSolid))){
			VSolid=true;
		}else{
			VSolid=false;
		}

		if(V.Y<0 && !VSolid){					//if moving in -y and neighbour is not solid

			P->Y = (double)(C.yn+1);			//Move particle to voxel boundary and inside left voxel
			P->VoxelY = C.yn;

		}else if(V.Y<0 && VSolid){				//if moving in -y and neighbour is solid

			P->Y = (double)CurrentVoxelY;		//Move particle to border but inside current voxel

		}else if(V.Y>0 && !VSolid){				//if moving in +y and neighbour is not solid

			P->Y = (double)C.yp;				//Move particle to voxel boundary and into right voxel
			P->VoxelY = C.yp;

		}else{									//if moving in +y and neighbour is solid

			P->Y = (double)(CurrentVoxelY+1);	//Move particle to voxel boundary but inside current voxel

		}

			P->X += V.X*Ytime;				//Move in proportion to Vy moved
			P->Z += V.Z*Ytime;

			V.X *= (1-Ytime);				//Reduce walk vector by amount moved
			V.Y *= (1-Ytime);
			V.Z *= (1-Ytime);
		
		if(VSolid){
			V.Y = -V.Y;					//If boundary is solid, reverse component (reflection)
		}

	}else{								//Reaches x boundary first

		bool VSolid;					//Is bounding voxel solid

		if((V.X<0 && (LatticeElement->Solid&NegativeXSolid)) || (V.X>=0 && (LatticeElement->Solid&PositiveXSolid))){
			VSolid=true;
		}else{
			VSolid=false;
		}

		if(V.X<0 && !VSolid){					//if moving in -x and neighbour is not solid

			P->X = (double)(C.xn+1);			//Move particle to voxel boundary and inside left voxel
			P->VoxelX = C.xn;

		}else if(V.X<0 && VSolid){				//if moving in -x and neighbour is solid

			P->X = (double)CurrentVoxelX;		//Move particle to border but inside current voxel

		}else if(V.X>0 && !VSolid){				//if moving in +x and neighbour is not solid

			P->X = (double)C.xp;				//Move particle to voxel boundary and into right voxel
			P->VoxelX = C.xp;

		}else{									//if moving in +x and neighbour is solid

			P->X = (double)(CurrentVoxelX+1);	//Move particle to voxel boundary but inside current voxel

		}

			P->Y += V.Y*Xtime;				//Move in proportion to Vx moved
			P->Z += V.Z*Xtime;

			V.X *= (1-Xtime);				//Reduce walk vector by amount moved
			V.Y *= (1-Xtime);
			V.Z *= (1-Xtime);
		
		if(VSolid){
			V.X = -V.X;		//If boundary is solid, reverse component (reflection)
		}

	}

	int x = P->VoxelX;
	int y = P->VoxelY;
	int z = P->VoxelZ;

	if(x==-1){									//Transfer into -ve x domain
		*ParticleStatus = 1;
		return;
	}else if(x==ProcessorDomain.DomainWidthX){	//Transfer into +ve x domain
		*ParticleStatus = 2;
		return;
	}

	if(y==-1){									//Transfer into -ve y domain
		*ParticleStatus = 3;
		return;
	}else if(y==ProcessorDomain.DomainWidthY){	//Transfer into +ve y domain
		*ParticleStatus = 4;
		return;
	}

	if(z==-1){									//Transfer into -ve z domain
		*ParticleStatus = 5;
		return;
	}else if(z==ProcessorDomain.DomainWidthZ){	//Transfer into +ve z domain
		*ParticleStatus = 6;
		return;
	}

	}	//end of while loop


}

#endif

void CaseNoSolids(int i, Coords* C, Voxel* LatticeElement, double* dt){

	int x = C->x;
	int y = C->y;
	int z = C->z;

	double u1 = LatticeElement->u1;
	double u2 = LatticeElement->u2;
	double v1 = LatticeElement->v1;
	double v2 = LatticeElement->v2;
	double w1 = LatticeElement->w1;
	double w2 = LatticeElement->w2;

	double x1 = (double)x;
	double x2 = (double)(x1+1);
	double y1 = (double)y;
	double y2 = (double)(y1+1);
	double z1 = (double)z;
	double z2 = (double)(z1+1);

	Particle* P = &(Particles[i]);

	double Px = P->X;					//X coordinate of particle
	double Py = P->Y;					//Y coordinate of particle
	double Pz = P->Z;					//Z coordinate of particle

	//cout << "Position: (" << Px << ", " << Py << ", " << Pz << "), Voxel: (" << x << ", " << y << ", " << z << ")" << endl;

	double Vx = (Px-x1)*(u2-u1) + u1;		//Velocity in x = ((Px-x1)/dx)*(u2-u1) + u1;
	double Vy = (Py-y1)*(v2-v1) + v1;		//Velocity in y = ((Py-y1)/dy)*(v2-v1) + v1;
	double Vz = (Pz-z1)*(w2-w1) + w1;		//Velocity in z = ((Pz-z1)/dz)*(w2-w1) + w1;

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 0 Solids:" << endl;
		AdvectDebug << "u1 = " << u1 << endl << "u2 = " << u2 << endl << "v1 = " << v1 << endl << "v2 = " << v2 << endl << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(Vx,Vy,Vz) = (" << Vx << ", " << Vy << ", " << Vz << ")" << endl;
	}

	//Calculate tau, time to exit voxel in each direction

	//cout << "Px: " << Px << ", x1: " << x1 << endl;
	//cout << "u1 = " << u1 << "; u2 = " << u2 << endl;
	//cout << "(fabs(u2)-fabs(u1)) = " << (fabs(u2)-fabs(u1)) << endl;

	double tau_x;
	if(Vx==0 || (u2<0 && u1>0)){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 1" << endl; }
			tau_x = DBL_MAX;
	}else if(fabs(u2-u1)>Min_dV){		
		if(Vx>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 2" << endl; }
			tau_x = log(u2/(u1 + (u2-u1)*(Px-x1)))/(u2-u1);		// = (dx/(u2-u1))*log((u2*dx)/(u1*dx+(u2-u1)*(Px-x1)));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 3" << endl; }
			tau_x = log(u1/(u2 + (u1-u2)*(x2-Px)))/(u2-u1);		// = (dx/(u2-u1))*log((u1*dx)/(u2*dx+(u1-u2)*(x2-Px)));
		}
	}else{
		if(Vx>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 4" << endl; }
			tau_x =	(x2-Px)/Vx;									// = 2*(x2-Px)/(u1+u2);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 5" << endl; }
			tau_x = (x1-Px)/Vx;									// = 2*(x1-Px)/(u1+u2);
		}
	}

	double tau_y;

	if(Vy==0 || (v2<0 && v1>0)){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 1" << endl; }
			tau_y = DBL_MAX;
	}else if(fabs(v2-v1)>Min_dV){
		if(Vy>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 2" << endl; }
			tau_y = log(v2/(v1 + (v2-v1)*(Py-y1)))/(v2-v1);	// = (dy/(v2-v1))*log((v2*dy)/(v1*dy+(v2-v1)*(Py-y1)));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 3" << endl; }
			tau_y = log(v1/(v2 + (v1-v2)*(y2-Py)))/(v2-v1);	// = (dy/(v2-v1))*log((v1*dy)/(v2*dy+(v1-v2)*(y2-Py)));
		}
	}else{
		if(Vy>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 4" << endl; }
			tau_y = (y2-Py)/Vy;								// = 2*(y2-Py)/(v1+v2);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 5" << endl; }
			tau_y = (y1-Py)/Vy;								// = 2*(y1-Py)/(v1+v2);
		}
	}
		
	double tau_z;
	if(Vz==0 || (w2<0 && w1>0)){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
			tau_z = DBL_MAX;
	}else if(fabs(w2-w1)>Min_dV){
		if(Vz>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
			tau_z = log(w2/(w1 + (w2-w1)*(Pz-z1)))/(w2-w1);	// = (dz/(w2-w1))*log((w2*dz)/(w1*dz+(w2-w1)*(Pz-z1)));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 3" << endl; }
			tau_z = log(w1/(w2 + (w1-w2)*(z2-Pz)))/(w2-w1);	// = dz/(w2-w1))*log((w1*dz)/(w2*dz+(w1-w2)*(z2-Pz)));
		}
	}else{
		if(Vz>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 4" << endl; }
			tau_z = (z2-Pz)/Vz;								// = 2*(z2-Pz)/(w1+w2);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 5" << endl; }
			tau_z = (z1-Pz)/Vz;								// = 2*(z1-Pz)/(w1+w2);
		}
	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: (" << tau_x << ", " << tau_y << ", " << tau_z << ")" << endl;
	}

	//cout << "Taus: (" << tau_x << ", " << tau_y << ", " << tau_z << ")" << endl;

	double time;								//Time to advect
	bool ReachesBorderX=false;
	bool ReachesBorderY=false;
	bool ReachesBorderZ=false;

	if(tau_x<=tau_y && tau_x<=tau_z){			//Reaches x boundary first

		if(tau_x>(*dt)){						//Remains in same voxel

			time = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			time = tau_x;						//Advect for tau_x
			(*dt) -= tau_x;						//Reduce dt by tau_x time used

			if(Vx<0){							//Reaches border with -ve neighbour
				P->VoxelX = C->xn;				//Set particle at border and in -ve voxel
				P->X = (double)(C->xn+1);
			}else{
				P->VoxelX = C->xp;				//Set particle at border and in +ve voxel
				P->X = (double)(C->xp);
			}

			ReachesBorderX=true;

		}

	}else if(tau_y<=tau_x && tau_y<=tau_z){		//Reaches y boundary first

		if(tau_y>(*dt)){						//Remains in same voxel

			time = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			time = tau_y;						//Advect for tau_y
			(*dt) -= tau_y;						//Reduce dt by tau_y time used

			if(Vy<0){							//Reaches border with -ve neighbour
				P->VoxelY = C->yn;				//Set particle at border and in -ve voxel
				P->Y = (double)(C->yn+1);
			}else{
				P->VoxelY = C->yp;				//Set particle at border and in +ve voxel
				P->Y = (double)(C->yp);
			}

			ReachesBorderY=true;

		}

	}else if(tau_z<=tau_x && tau_z<=tau_y){		//Reaches z boundary first

		if(tau_z>(*dt)){						//Remains in same voxel

			time = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			time = tau_z;						//Advect for tau_z
			(*dt) -= tau_z;						//Reduce dt by tau_z time used

			if(Vz<0){							//Reaches border with -ve neighbour
				P->VoxelZ = C->zn;	//Set particle at border and in -ve voxel
				P->Z = (double)(C->zn+1);
			}else{
				P->VoxelZ = C->zp;	//Set particle at border and in +ve voxel
				P->Z = (double)(C->zp);
			}

			ReachesBorderZ=true;

		}

	}

	double xe, ye, ze;	//New positions

	if(!ReachesBorderX){		//If reached boundary, new position has already been set

		if(Vx==0){
			xe = Px;
		}else if(fabs(u2-u1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 1" << endl; }
			double u1_u2u1 = u1/(u2-u1);
			xe = x1 - u1_u2u1 + (u1_u2u1 + (Px-x1)) * exp((u2-u1)*time);	// = x1 - (u1*dx/(u2-u1)) + ((u1*dx/(u2-u1))+(Px-x1)) * exp((u2-u1)*time/dx);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 2" << endl; }
			xe = Px + time*Vx;
		}

		P->X = xe;
	}

	if(!ReachesBorderY){		//If reached boundary, new position has already been set

		if(Vy==0){
			ye = Py;
		}else if(fabs(v2-v1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 1" << endl; }
			double v1_v2v1 = v1/(v2-v1);
			ye = y1 - v1_v2v1 + (v1_v2v1 + (Py-y1)) * exp((v2-v1)*time);	// = y1 - (v1*dy/(v2-v1)) + ((v1*dy/(v2-v1))+(Py-y1)) * exp((v2-v1)*time/dy);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 2" << endl; }
			ye = Py + time*Vy;
		}

		P->Y = ye;
	}

	if(!ReachesBorderZ){		//If reached boundary, new position has already been set

		if(Vz==0){
			ze = Pz;
		}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = z1 - w1_w2w1 + (w1_w2w1 + (Pz-z1)) * exp((w2-w1)*time);	// = z1 - (w1*dz/(w2-w1)) + ((w1*dz/(w2-w1))+(Pz-z1)) * exp((w2-w1)*time/dz);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
			ze = Pz + time*Vz;
		}

		P->Z = ze;
	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "New Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl << endl;
	}

}

struct Case1SolidCoords{

double Velocity_a1;				//Inward velocity opposite solid
double Velocity_b1;				//Inward velocity to voxel in b direction
double Velocity_b2;				//Outward velocity to voxel in b direction
double Velocity_c1;				//Inward velocity to voxel in c direction
double Velocity_c2;				//Outward velocity to voxel in c direction

double Pa_a1;					//Position in a direction within voxel (Pa-a1)
double Pb_b1;					//Position in b direction within voxel (Pb-b1)
double Pc_c1;					//Position in c direction within voxel (Pc-c1)

bool ReachedNegativeABorder;
bool ReachedPositiveBBorder;	//Function returns whether particle reached border
bool ReachedNegativeBBorder;
bool ReachedPositiveCBorder;
bool ReachedNegativeCBorder;

double* dt;						//Advection time pointer
};


//Advects particles in a voxel with two opposite solids on the a1 and a2 faces, used as a general case
//Dimensions referred to as x,y,z in the function, but correspond to a,b,c in the struct to reduce
//  confusion when transforming into the general case in Case2Solids()
void Case1Solid(Case1SolidCoords* Case){

	double u1 = Case->Velocity_a1;

	double v1 = Case->Velocity_b1;
	double v2 = Case->Velocity_b2;

	double w1 = Case->Velocity_c1;
	double w2 = Case->Velocity_c2;

	double Px = Case->Pa_a1;
	double Py = Case->Pb_b1;
	double Pz = Case->Pc_c1;

	double x1 = 0;	//Imaginary voxel at 0,0,0
	double x2 = 1;
	double y1 = 0;
	double y2 = 1;
	double z1 = 0;
	double z2 = 1;

	double Tx;
	double Ty;
	double Tz;

	double* dt = Case->dt;

	double v_face = (v1 + (v2-v1)*Py);
	double w_face = (w1 + (w2-w1)*Pz);

	double u = u1*(1-Px)*(1-Px);			// = u1*((x2-Px)*(x2-Px)/(dx*dx));
	double v = 2*(1-Px)*v_face;				// = 2*v1*(x2-Px)/dx + 2*(v2-v1)*(x2-Px)*(Py-y1)/(dx*dy);
	double w = 2*(1-Px)*w_face;				// = 2*w1*(x2-Px)/dx + 2*(w2-w1)*(x2-Px)*(Pz-z1)/(dx*dz);

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 1 Solid:" << endl;
		AdvectDebug << "u1 = " << u1 << endl << "v1 = " << v1 << endl << "v2 = " << v2 << endl << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	if(u<0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 1" << endl; }
			Tx = (1 - (1/(1-Px)))/u1;					//Calculate tau_x = (dx*dx/u1)*((1/dx) - (1/(x2-Px)));
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 2" << endl; }
			Tx = DBL_MAX;
	}

	if(v==0 || (v2<0 && v1>0)){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 1" << endl; }
				Ty = DBL_MAX;
	}else if(fabs(v2-v1)>Min_dV){						//Calculate tau_y
		if(v>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 2" << endl; }
				Ty = (1-Py)/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 3" << endl; }
				Ty = pow( v2/(v1 + (v2-v1)*Py) , 0.5*u1/(v2-v1) )/(u1*(1-Px)) - 1/(u1*(1-Px));	// = fabs((dx*dx/(u1*(x2-Px))) * pow( v2*dy/(v1*dy+(v2-v1)*(Py-y1)), u1*dy/(2*dx*(v2-v1)) ) - dx*dx/(u1*(x2-Px)));
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 4" << endl; }
				Ty = -Py/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 5" << endl; }
				Ty = pow( v1/(v1 + (v2-v1)*Py) , 0.5*u1/(v2-v1) )/(u1*(1-Px)) - 1/(u1*(1-Px));	// = fabs((dx*dx/(u1*(x2-Px))) * pow( v2*dy/(v1*dy+(v2-v1)*(y2-Py)), u1*dy/(2*dx*(v2-v1)) ) - dx*dx/(u1*(x2-Px)));
				//cout << "Pow( " << v1/(v1 + (v2-v1)*Py) << " , " << 0.5*u1/(v2-v1) << ") / " << (u1*(1-Px)) - 1/(u1*(1-Px)) << endl;
			}
		}
	}else{
		if(v>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 6" << endl; }
				Ty = (1-Py)/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 7" << endl; }
				Ty = (exp( (u1*(1-Py))/(2*v_face) ) - 1)/(u1*(1-Px));	//Modification: (v1+v2) -> v_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(y2-Py))/(dx*(v1+v2))) - 1);
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 8" << endl; }
				Ty = -Py/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 9" << endl; }
				Ty = (exp( -(u1*Py)/(2*v_face) ) - 1)/(u1*(1-Px));	//Modification: (v1+v2) -> v_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(y1-Py))/(dx*(v1+v2))) - 1);
			}
		}
	}

	if(w==0 || (w2<0 && w1>0)){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
				Tz = DBL_MAX;
	}else if(fabs(w2-w1)>Min_dV){								//Calculate tau_z
		if(w>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
				Tz = (1-Pz)/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 3" << endl; }
				Tz = fabs(pow( w2/(w1 + (w2-w1)*Pz) , 0.5*u1/(w2-w1) )/(u1*(1-Px)) - 1/(u1*(1-Px)));		// = fabs((dx*dx/(u1*(x2-Px))) * pow( w2*dz/(w1*dz+(w2-w1)*(Pz-z1)), u1*dz/(2*dx*(w2-w1)) ) - dx*dx/(u1*(x2-Px)));
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 4" << endl; }
				Tz = -Pz/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 5" << endl; }
				Tz = fabs(pow( w1/(w1 + (w2-w1)*Pz) , 0.5*u1/(w2-w1) )/(u1*(1-Px)) - 1/(u1*(1-Px)));	// = fabs((dx*dx/(u1*(x2-Px))) * pow( w2*dz/(w1*dz+(w2-w1)*(z2-Pz)), u1*dz/(2*dx*(w2-w1)) ) - dx*dx/(u1*(x2-Px)));
			}
		}
	}else{
		if(w>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 6" << endl; }
				Tz = (1-Pz)/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 7" << endl; }
				Tz = (exp( (u1*(1-Pz))/(2*w_face) ) - 1)/(u1*(1-Px));	//Modification: (w1+w2) -> w_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(z2-Pz))/(dx*(w1+w2))) - 1);
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 8" << endl; }
				Tz = -Pz/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 9" << endl; }
				Tz = (exp( -(u1*Pz)/(2*w_face) ) - 1)/(u1*(1-Px));	//Modification: (w1+w2) -> w_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(z1-Pz))/(dx*(w1+w2))) - 1);
			}
		}
	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

	// cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	double t;

	if(Tx<=Ty && Tx<=Tz){			//Reaches x boundary first

		//cout<<"Min tau=x"<<endl;

		if(Tx>(*dt)){							//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
	
		}else{									//Advect to the border
	
			t = Tx;							//Advect for tau_x
			(*dt) -= Tx;						//Reduce dt by tau_x time used

			Case->ReachedNegativeABorder=true;

		}

	}else if(Ty<=Tx && Ty<=Tz){		//Reaches y boundary first

	//	cout<<"Min tau=y"<<endl;

		if(Ty>(*dt)){							//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
	
		}else{									//Advect to the border
	
			t = Ty;							//Advect for tau_y
			(*dt) -= Ty;						//Reduce dt by tau_y time used

			if(v<0){							//Reaches border with -ve neighbour
				Case->ReachedNegativeBBorder=true;
			}else{
				Case->ReachedPositiveBBorder=true;
			}

		}

	}else if(Tz<=Tx && Tz<=Ty){		//Reaches z boundary first
	
		//cout<<"Min tau=z"<<endl;

		if(Tz>(*dt)){						//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
	
		}else{									//Advect to the border
	
			t = Tz;							//Advect for tau_z
			(*dt) -= Tz;						//Reduce dt by tau_z time used

			if(w<0){							//Reaches border with -ve neighbour
				Case->ReachedNegativeCBorder=true;
			}else{
				Case->ReachedPositiveCBorder=true;
			}

		}

	}
		

	//cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" <<endl;


	double xe;		//Positions after advection
	double ye;
	double ze;	

	if(Case->ReachedNegativeABorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 1" << endl; }
			xe = 0;
	}else if(u1==0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 2" << endl; }
			xe = Px;
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 3" << endl; }
			xe = 1 - 1/(t*u1 + 1/(1-Px));		// = x2 - 1/(t*u1/(dx*dx)+1/(x2-Px));
	}

	if(Case->ReachedNegativeBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 1" << endl; }
			ye = 0;
	}else if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 2" << endl; }
			ye = 1;
	}else if(fabs(v2-v1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 3" << endl; }
			double v1_v2v1 = v1/(v2-v1);
			ye = -v1_v2v1 + (v1_v2v1 + Py) * pow( 1+(u1*(1-Px)*t) , 2*(v2-v1)/u1 );		// = y1 - (v1*dy/(v2-v1)) + ((v1*dy+(v2-v1)*(Py-y1))/(v2-v1)) * pow( 1+(u1*t*(x2-Px)/(dx*dx)) , 2*dx*(v2-v1)/(u1*dy) );
	}else{
		if(u1==0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 4" << endl; }
			ye = Py + v*t;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 5" << endl; }
			ye = Py + log( 1 + (u1*(1-Px)*t) )*(v1+v2)/u1;		// = Py + (dx*(v1+v2)/u1) * log((u1*(x2-Px)*t)/(dx*dx) + 1);
		}
	}

	if(Case->ReachedNegativeCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
			ze = 0;
	}else if(Case->ReachedPositiveCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
			ze = 1;
	}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 3" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = -w1_w2w1 + (w1_w2w1 + Pz) * pow( 1+(u1*(1-Px)*t) , 2*(w2-w1)/u1 );	// = z1 - (w1*dz/(w2-w1)) + ((w1*dz+(w2-w1)*(Pz-z1))/(w2-w1)) * pow((1+(u1*t*(x2-Px)/(dx*dx))) , (2*dx*(w2-w1)/(u1*dz)));
	}else{
		if(u1==0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 4" << endl; }
			ze = Pz + w*t;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 5" << endl; }
			ze = Pz + log( 1 + (u1*(1-Px)*t) )*(w1+w2)/u1;		// = Pz + (dx*(w1+w2)/u1) * log((u1*(x2-Px)*t)/(dx*dx) + 1);
		}
	}

	/*

	if(xe<0||xe>1.00000001||ye<0||ye>1.00000001||ze<0||ze>1.00000001){
		cout << "1 Solid Fault - Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}
	*/

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

	Case->Pa_a1 = xe;
	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

}

//Advection function used when one neighbouring voxel is solid (6 cases)
void CaseOneSolid(int i,Coords* C, Voxel* LatticeElement, double* dt){

	int x = C->x;
	int y = C->y;
	int z = C->z;

	double u1 = LatticeElement->u1;
	double u2 = LatticeElement->u2;
	double v1 = LatticeElement->v1;
	double v2 = LatticeElement->v2;
	double w1 = LatticeElement->w1;
	double w2 = LatticeElement->w2;

	double x1 = (double)x;
	double x2 = (double)(x1+1);
	double y1 = (double)y;
	double y2 = (double)(y1+1);
	double z1 = (double)z;
	double z2 = (double)(z1+1);

	Particle* P = &(Particles[i]);

	double Px = P->X;					//X coordinate of particle
	double Py = P->Y;					//Y coordinate of particle
	double Pz = P->Z;					//Z coordinate of particle
	
	bool ReachesBorderX=false;
	bool ReachesBorderY=false;
	bool ReachesBorderZ=false;

	double xe, ye, ze;			//New positions after advection


	Case1SolidCoords coords;	//Struct for solids function

	coords.dt = dt;
	coords.ReachedNegativeABorder=false;
	coords.ReachedPositiveBBorder=false;
	coords.ReachedNegativeBBorder=false;
	coords.ReachedPositiveCBorder=false;
	coords.ReachedNegativeCBorder=false;

	//cout << "Position: (" << Px << ", " << Py << ", " << Pz << "), Voxel: (" << x << ", " << y << ", " << z << ")" << endl;

	//calculate tau, time to exit voxel in each direction

	unsigned char Solids = (unsigned char)(LatticeElement->Solid);

	if(Solids&PositiveXSolid){		//Solid on the positive x side
		//a=x, b=y, c=z

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid" << endl;
		}

		coords.Velocity_a1 = u1;				//Inward velocity opposite solid
		coords.Velocity_b1 = v1;				//Inward velocity to voxel in b direction
		coords.Velocity_b2 = v2;				//Outward velocity to voxel in b direction
		coords.Velocity_c1 = w1;				//Inward velocity to voxel in c direction
		coords.Velocity_c2 = w2;				//Outward velocity to voxel in c direction
		
		coords.Pa_a1 = Px - x1;
		coords.Pb_b1 = Py - y1;
		coords.Pc_c1 = Pz - z1;

		Case1Solid(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}

		if(coords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		if(coords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = coords.Pa_a1 + x1;
		ye = coords.Pb_b1 + y1;
		ze = coords.Pc_c1 + z1;
		
		
	}else if(Solids&PositiveYSolid){//Solid on the positive y side
			//a=y, b=-x, c=z

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Y Solid" << endl;
		}

		coords.Velocity_a1 = v1;				//Inward velocity opposite solid
		coords.Velocity_b1 = -u2;				//Inward velocity to voxel in b direction
		coords.Velocity_b2 = -u1;				//Outward velocity to voxel in b direction
		coords.Velocity_c1 = w1;				//Inward velocity to voxel in c direction
		coords.Velocity_c2 = w2;				//Outward velocity to voxel in c direction
		
		coords.Pa_a1 = Py - y1;
		coords.Pb_b1 = x2 - Px;
		coords.Pc_c1 = Pz - z1;

		Case1Solid(&coords);

		if(coords.ReachedNegativeABorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}

		if(coords.ReachedPositiveBBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		if(coords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = x2 - coords.Pb_b1;
		ye = coords.Pa_a1 + y1;
		ze = coords.Pc_c1 + z1;

	}else if(Solids&PositiveZSolid){//Solid on the positive z side
		//a=z, b=y, c=-x

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid" << endl;
		}

		coords.Velocity_a1 = w1;				//Inward velocity opposite solid
		coords.Velocity_b1 = v1;				//Inward velocity to voxel in b direction
		coords.Velocity_b2 = v2;				//Outward velocity to voxel in b direction
		coords.Velocity_c1 = -u2;				//Inward velocity to voxel in c direction
		coords.Velocity_c2 = -u1;				//Outward velocity to voxel in c direction
		
		coords.Pa_a1 = Pz - z1;
		coords.Pb_b1 = Py - y1;
		coords.Pc_c1 = x2 - Px;

		Case1Solid(&coords);

		if(coords.ReachedNegativeABorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		if(coords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}

		if(coords.ReachedNegativeCBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
		}


		xe = x2 - coords.Pc_c1;
		ye = coords.Pb_b1 + y1;
		ze = coords.Pa_a1 + z1;

	}else if(Solids&NegativeXSolid){//Solid on the negative x side
		//a=-x, b=-y, c=z

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid" << endl;
		}

		coords.Velocity_a1 = -u2;				//Inward velocity opposite solid
		coords.Velocity_b1 = -v2;				//Inward velocity to voxel in b direction
		coords.Velocity_b2 = -v1;				//Outward velocity to voxel in b direction
		coords.Velocity_c1 = w1;				//Inward velocity to voxel in c direction
		coords.Velocity_c2 = w2;				//Outward velocity to voxel in c direction
		
		coords.Pa_a1 = x2 - Px;
		coords.Pb_b1 = y2 - Py;
		coords.Pc_c1 = Pz - z1;

		Case1Solid(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
		}

		if(coords.ReachedPositiveBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		if(coords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = x2 - coords.Pa_a1;
		ye = y2 - coords.Pb_b1;
		ze = coords.Pc_c1 + z1;
		


	}else if(Solids&NegativeYSolid){//Solid on the negative y side

		//a=-y, b=x, c=z

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid" << endl;
		}

		coords.Velocity_a1 = -v2;				//Inward velocity opposite solid
		coords.Velocity_b1 = u1;				//Inward velocity to voxel in b direction
		coords.Velocity_b2 = u2;				//Outward velocity to voxel in b direction
		coords.Velocity_c1 = w1;				//Inward velocity to voxel in c direction
		coords.Velocity_c2 = w2;				//Outward velocity to voxel in c direction
		
		coords.Pa_a1 = y2 - Py;
		coords.Pb_b1 = Px - x1;
		coords.Pc_c1 = Pz - z1;

		Case1Solid(&coords);

		if(coords.ReachedNegativeABorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}

		if(coords.ReachedPositiveBBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
						
		}

		if(coords.ReachedNegativeBBorder){
			
			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		if(coords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = coords.Pb_b1 + x1;
		ye = y2 - coords.Pa_a1;
		ze = coords.Pc_c1 + z1;
	

	}else if(Solids&NegativeZSolid){//Solid on the negative z side
		//a=-z, b=y, c=x

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid" << endl;
		}

		coords.Velocity_a1 = -w2;				//Inward velocity opposite solid
		coords.Velocity_b1 = v1;				//Inward velocity to voxel in b direction
		coords.Velocity_b2 = v2;				//Outward velocity to voxel in b direction
		coords.Velocity_c1 = u1;				//Inward velocity to voxel in c direction
		coords.Velocity_c2 = u2;				//Outward velocity to voxel in c direction
		
		coords.Pa_a1 = z2 - Pz;
		coords.Pb_b1 = Py - y1;
		coords.Pc_c1 = Px - x1;

		Case1Solid(&coords);

		if(coords.ReachedNegativeABorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		if(coords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
			
		}

		if(coords.ReachedNegativeCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}


		xe = coords.Pc_c1 + x1;
		ye = coords.Pb_b1 + y1;
		ze = z2 - coords.Pa_a1;

	}		
		

	if(!ReachesBorderX){		//If reached boundary, new position has already been set
		P->X = xe;
	}
	if(!ReachesBorderY){
		P->Y = ye;
	}
	if(!ReachesBorderZ){
		P->Z = ze;
	}

	//cout << "New position: (" << Particles[i].X << ", " << Particles[i].Y << ", " << Particles[i].Z << ")" << endl << endl;

}

struct Case2OppositeSolidsCoords{
double Velocity_b1;				//Inward velocity to voxel in b direction
double Velocity_b2;				//Outward velocity to voxel in b direction
double Velocity_c1;				//Inward velocity to voxel in c direction
double Velocity_c2;				//Outward velocity to voxel in c direction

double Pa_a1;					//Position in a direction within voxel (Pa-a1)
double Pb_b1;					//Position in b direction within voxel (Pb-b1)
double Pc_c1;					//Position in c direction within voxel (Pc-c1)

bool ReachedPositiveBBorder;	//Function returns whether particle reached border
bool ReachedNegativeBBorder;
bool ReachedPositiveCBorder;
bool ReachedNegativeCBorder;

double* dt;						//Advection time pointer
};


//Advects particles in a voxel with two opposite solids on the a1 and a2 faces, used as a general case
//Dimensions referred to as x,y,z in the function, but correspond to a,b,c in the struct to reduce
//  confusion when transforming into the general case in Case2Solids()
void Case2OppositeSolids(Case2OppositeSolidsCoords* Case){

	double v1 = Case->Velocity_b1;
	double v2 = Case->Velocity_b2;

	double w1 = Case->Velocity_c1;
	double w2 = Case->Velocity_c2;

	double Px = Case->Pa_a1;
	double Py = Case->Pb_b1;
	double Pz = Case->Pc_c1;

	double x1 = 0;	//Imaginary voxel at 0,0,0
	double x2 = 1;
	double y1 = 0;
	double y2 = 1;
	double z1 = 0;
	double z2 = 1;

	double Tx;
	double Ty;
	double Tz;

	double* dt = Case->dt;

	double _1_PxPx = (1-Px)*Px;					//Ubiquitous quantity, precalculate

	double u = 0;
	double v = 6*_1_PxPx*(v1 + (v2-v1)*Py);	// = 6*v1*(x2-Px)*(Px-x1)/(dx*dx) + 6*(v2-v1)*(x2-Px)*(Px-x1)*(Py-y1)/(dx*dx*dy);
	double w = 6*_1_PxPx*(w1 + (w2-w1)*Pz);	// = 6*w1*(x2-Px)*(Px-x1)/(dx*dx) + 6*(w2-w1)*(x2-Px)*(Px-x1)*(Pz-z1)/(dx*dx*dz);

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 2 Opposite Solids:" << endl;
		AdvectDebug << "v1 = " << v1 << endl << "v2 = " << v2 << endl << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 1" << endl; }
			Tx = DBL_MAX;			//X will not reach border
	
	if(v==0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 1" << endl; }
			Ty = DBL_MAX;
	}else if(fabs(v2-v1)>Min_dV){	//Expression fails if denominator too small, treat as a special case if velocities close
		if(v>0 && v2<=0){			//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 2" << endl; }
			Ty = DBL_MAX;
		}else if(v>0){				//Exit time to v2 face
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 3" << endl; }
			Ty = log( v2/(v1 + (v2-v1)*Py) )/(6*(v2-v1)*_1_PxPx);	// = (dx*dx*dy/(6*(v2-v1)*(x2-Px)*(Px-x1))) * log(v2*dy/(v1*dy+(v2-v1)*(Py-y1)));
		}else if(v<0 && v1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 4" << endl; }
			Ty = DBL_MAX;
		}else{						//Exit time to v1 face if v<0
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 5" << endl; }
			Ty = log( v1/(v1 + (v2-v1)*Py) )/(6*(v2-v1)*_1_PxPx);	// = (dx*dx*dy/(6*(v2-v1)*(x2-Px)*(Px-x1))) * -log(v2*dy/(v1*dy+(v2-v1)*(y2-Py)));
		}
	}else{
		if(v>0 && v2<=0){			//Take care of bizarre cases for rigour's sake - outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 6" << endl; }
			Ty = DBL_MAX;
		}else if(v>0){				//Exit time to v2 face
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 7" << endl; }
			Ty = (1-Py)/v;
		}else if(v<0 && v1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 8" << endl; }
			Ty = DBL_MAX;
		}else{						//Exit time to v1 face if v<0
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 9" << endl; }
			Ty = -Py/v;
		}
	}


	if(w==0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
			Tz = DBL_MAX;
	}else if(fabs(w2-w1)>Min_dV){
		if(w>0 && w2<=0){			//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
			Tz = DBL_MAX;
		}else if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 3" << endl; }
			Tz = log( w2/(w1 + (w2-w1)*Pz) )/(6*(w2-w1)*_1_PxPx);	// = (dx*dx*dz/(6*(w2-w1)*(x2-Px)*(Px-x1))) * log(w2*dz/(w1*dz+(w2-w1)*(Pz-z1)));
		}else if(w<0 && w1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 4" << endl; }
			Tz = DBL_MAX;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 5" << endl; }
			Tz = log( w1/(w1 + (w2-w1)*Pz) )/(6*(w2-w1)*_1_PxPx);	// = (dx*dx*dz/(6*(w2-w1)*(x2-Px)*(Px-x1))) * -log(w2*dz/(w1*dz+(w2-w1)*(z2-Pz)));
		}
	}else{
		if(w>0 && w2<=0){			//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 6" << endl; }
			Tz = DBL_MAX;
		}else if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 7" << endl; }
			Tz = (1-Pz)/w;
		}else if(w<0 && w1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 8" << endl; }
			Tz = DBL_MAX;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 9" << endl; }
			Tz = -Pz/w;
		}
	}


	double t;

	if(Ty<=Tx && Ty<=Tz){			//Reaches y boundary first

		if(Ty>(*dt)){						//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Ty;							//Advect for tau_y
			(*dt) -= Ty;						//Reduce dt by tau_y time used

			if(v<0){							//Reaches border with -ve neighbour
				Case->ReachedNegativeBBorder=true;
			}else{
				Case->ReachedPositiveBBorder=true;
			}

		}

	}else if(Tz<=Tx && Tz<=Ty){		//Reaches z boundary first

		if(Tz>(*dt)){						//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Tz;							//Advect for tau_z
			(*dt) -= Tz;						//Reduce dt by tau_z time used

			if(w<0){							//Reaches border with -ve neighbour
				Case->ReachedNegativeCBorder=true;
			}else{
				Case->ReachedPositiveCBorder=true;
			}

		}

	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

	//cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 1" << endl; }

	double xe = Px;		//Positions after advection
	double ye;
	double ze;	

	if(Case->ReachedNegativeBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 1" << endl; }
		ye = 0;

	}else if(Case->ReachedPositiveBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 2" << endl; }
		ye = 1;

	}else if(fabs(v2-v1)>Min_dV){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 3" << endl; }
		double v1_v2v1 = v1/(v2-v1);
		ye = y1 - v1_v2v1 + (v1_v2v1 + Py) * exp( 6*(v2-v1)*_1_PxPx*t );	// = y1 - (v1*dy)/(v2-v1) + ( (v1*dy + (v2-v1)*(Py-y1))/((v2-v1)*dy) ) * exp( (6*(v2-v1)*(Px-x1)*(x2-Px)/(dx*dx*dy))*t );

	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 4" << endl; }	
		ye = Py + 3*(v1+v2)*_1_PxPx*t;		// = Py + ( 3*(v1+v2)*(x2-Px)*(Px-x1)/(dx*dx) )*t;

	}

	if(Case->ReachedNegativeCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
		ze = 0;

	}else if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
		ze = 1;

	}else if(fabs(w2-w1)>Min_dV){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 3" << endl; }
		double w1_w2w1 = w1/(w2-w1);
		ze = z1 - w1_w2w1 + (w1_w2w1 + Pz) * exp( 6*(w2-w1)*_1_PxPx*t );	// = z1 - (w1*dz)/(w2-w1) + ( (w1*dz + (w2-w1)*(Pz-z1))/((w2-w1)*dz) ) * exp( (6*(w2-w1)*(Px-x1)*(x2-Px)/(dx*dx*dz))*t );

	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 4" << endl; }
		ze = Pz + 3*(w1+w2)*_1_PxPx*t;	// = Pz + ( 3*(w1+w2)*(x2-Px)*(Px-x1)/(dx*dx) )*t;

	}

	/*

	if(ye<0||ye>1.00000001||ze<0||ze>1.00000001){
		cout << "2 Opposite Solids Fault - Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
		cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" <<endl;
		//Sleep(5000);
	}
	*/

	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

}

struct Case2AdjacentSolidsCoords{
double Velocity_a2;				//Outward velocity opposite solid a direction
double Velocity_b2;				//Outward velocity opposite solid b direction
double Velocity_c1;				//Inward velocity to c1 face
double Velocity_c2;				//Outward velocity to c2 face

double Pa_a1;					//Position in a direction within voxel (Pa-a1)
double Pb_b1;					//Position in b direction within voxel (Pb-b1)
double Pc_c1;					//Position in c direction within voxel (Pc-c1)

bool ReachedPositiveABorder;	//Function returns whether particle reached border
bool ReachedPositiveBBorder;
bool ReachedPositiveCBorder;
bool ReachedNegativeCBorder;

double* dt;						//Advection time pointer
};

//Advects particles in a voxel with two adjacent solids on the a1 and b1 faces, used as a general case
//Dimensions referred to as x,y,z in the function, but correspond to a,b,c in the struct to reduce
//  confusion when transforming into the general case in Case2Solids()
void Case2AdjacentSolids(Case2AdjacentSolidsCoords* Case){

	//cout << "Case 2 Adjacent Solids" << endl;

	double u2 = Case->Velocity_a2;

	double v2 = Case->Velocity_b2;

	double w1 = Case->Velocity_c1;
	double w2 = Case->Velocity_c2;

	double Px = Case->Pa_a1;
	double Py = Case->Pb_b1;
	double Pz = Case->Pc_c1;

	double x1 = 0;	//Imaginary voxel at 0,0,0
	double x2 = 1;
	double y1 = 0;
	double y2 = 1;
	double z1 = 0;
	double z2 = 1;

	double* dt = Case->dt;

	double _2PxPy = 2*Px*Py;

	double u = _2PxPy*u2*Px;						// = (2*u2)*(Px-x1)*(Px-x1)*(Py-y1)/(dx*dx*dy);
	double v = _2PxPy*v2*Py;						// = (2*v2)*(Px-x1)*(Py-y1)*(Py-y1)/(dx*dy*dy);
	double w = 2*_2PxPy*(w1 + (w2-w1)*Pz);			// = (4*w1)*(Px-x1)*(Py-y1)/(dx*dy) + 4*(w2-w1)*(Px-x1)*(Py-y1)*(Pz-z1)/(dx*dy*dz);

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 2 Adjacent Solids:" << endl;
		AdvectDebug << "u2 = " << u2 << endl << "v2 = " << v2 << endl << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	//cout << "Velocity: (" << u << ", " << v << ", " << w << ")" <<endl;

	double D = (u2 + v2)*_2PxPy;					// = 2*(u2*dy + v2*dx)*(Px-x1)*(Py-y1)/(dx*dx*dy*dy);
	double invD = 1/D;

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "D = " << D << endl << "invD = " << invD << endl;
	}

	double Tx;
	double Ty;
	double Tz;

	//Calculate time to reach the border, Tau
	
	//Tau x
	
	if(u > 0){
		if(u2 < Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 1" << endl; }
			Tx = (1-Px)/u;
		}else if(fabs(u2 + v2)>Min_dV){									// fabs(u2*dy + v2*dx)>Min_dV
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 2" << endl; }
			Tx = invD * (1 - pow( Px , (u2 + v2)/u2 ));				// = (1/D) * (1 - pow( dx/(Px-x1) , -(u2*dy + v2*dx)/(u2*dy) ));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 3" << endl; }
			Tx = -log( Px )/(u2*_2PxPy);							// = ((dx*dx*dy*dy)/(2*u2*(Px-x1)*(Py-y1))) * log(dx/(Px-x1));
		}
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 4" << endl; }
			Tx = DBL_MAX;
	}

	//Tau y

	if(v > 0){
		if(v2 < Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 1" << endl; }
			Ty = (1-Py)/v;
		}else if(fabs(u2 + v2)>Min_dV){									// fabs(u2*dy + v2*dx)>Min_dV
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 2" << endl; }
			Ty = invD * (1 - pow( Py , (u2 + v2)/v2 ));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 3" << endl; }
			Ty = -log( Py )/(v2*_2PxPy);							// = ((dx*dx*dy*dy)/(2*v2*(Px-x1)*(Py-y1))) * log(dy/(Py-y1));
		}
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 4" << endl; }
			Ty = DBL_MAX;
	}

	//Tau z
	if(w==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
				Tz = DBL_MAX;
	}else if(w>0){
		if(w2<=0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
				Tz = DBL_MAX;
		}else if(fabs(w2-w1)>Min_dV){
			if(fabs(u2 + v2)>Min_dV){														// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 3" << endl; }
				Tz = invD * (1 - pow( w2/(Pz*(w2-w1) + w1) , -(u2+v2)/(2*(w2-w1)) ));		// = (1/D)*(1 - pow( (dx*(w2-w1) + w1*dz)/((Pz-z1)*(w2-w1) + w1*dz) , -(dz*(u2*dy+v2*dx)/(2*(w2-w1)*dx*dy)) ));
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 4" << endl; }
				Tz = log( w2/(Pz*(w2-w1) + w1) )/(2*(w2-w1)*_2PxPy);						// = ((dx*dy*dz)/(4*(w2-w1)*(Px-x1)*(Py-y1))) * log((dx*(w2-w1)+(w1*dz))/((Pz-z1)*(w2-w1)+w1*dz));
			}
		}else{
			if(fabs(u2 + v2)>Min_dV){													// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 5" << endl; }
				Tz = invD * (1 - exp( -(1-Pz)*(u2+v2)/(2*w2) ));							// = (1/D) * (1 - exp( -(z2-Pz)*(u2*dy+v2*dx)/(2*w1*dx*dy) ));
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 6" << endl; }
				Tz = (1-Pz)/(2*w2*_2PxPy);													// = (x2-Pz)*dx*dy/(4*w1*(Px-x1)*(Py-y1));
			}
		}
	}else{
		if(w1>=0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 7" << endl; }
				Tz = DBL_MAX;
		}else if(fabs(w2-w1)>Min_dV){
			if(fabs(u2 + v2)>Min_dV){														// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 8" << endl; }
				Tz = invD * (1 - pow( w1/(Pz*(w2-w1) + w1) , -(u2+v2)/(2*(w2-w1)) ));		// = (1/D)*(1 - pow( (w1*dz)/((Pz-z1)*(w2-w1) + w1*dz) , -(dz*(u2*dy+v2*dx)/(2*(w2-w1)*dx*dy)) ));
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 9" << endl; }
				Tz = log( w1/(Pz*(w2-w1) + w1) )/(2*(w2-w1)*_2PxPy);						// = ((dx*dy*dz)/(4*(w2-w1)*(Px-x1)*(Py-y1))) * log((w1*dz)/((Pz-z1)*(w2-w1)+w1*dz));
			}
		}else{
			if(fabs(u2 + v2)>Min_dV){														// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 10" << endl; }
				Tz = invD * (1 - exp( Pz*(u2+v2)/(2*w1) ));									// = (1 - exp( -(z1-Pz)*(u2*dy+v2*dx)/(2*w1*dx*dy) ))/D;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 11" << endl; }
				Tz = -Pz/(2*w1*_2PxPy);														// = (z1-Pz)*dx*dy/(4*w1*(Px-x1)*(Py-y1));
			}
		}
	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}
	//cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	
	double t;

	if(Tx<=Ty && Tx<=Tz){			//Reaches x boundary first

		if(Tx>(*dt)){							//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Tx;							//Advect for tau_x
			(*dt) -= Tx;						//Reduce dt by tau_x time used

			Case->ReachedPositiveABorder=true;

		}

	}else if(Ty<=Tx && Ty<=Tz){		//Reaches y boundary first

		if(Ty>(*dt)){							//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Ty;							//Advect for tau_y
			(*dt) -= Ty;						//Reduce dt by tau_y time used

			Case->ReachedPositiveBBorder=true;

		}

	}else if(Tz<=Tx && Tz<=Ty){		//Reaches z boundary first

		if(Tz>(*dt)){						//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Tz;							//Advect for tau_z
			(*dt) -= Tz;						//Reduce dt by tau_z time used

			if(w<0){							//Reaches border with -ve neighbour
				Case->ReachedNegativeCBorder=true;
			}else{
				Case->ReachedPositiveCBorder=true;
			}

		}

	}

	//cout << "Analytical Tau:      (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	//Positions

	double xe, ye, ze;

	if(fabs(u2 + v2)>Min_dV){															// fabs(u2*dy + v2*dx)>Min_dV

		if(Case->ReachedPositiveABorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 1" << endl; }
			xe = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 2" << endl; }
			xe = Px * pow( 1-D*t , -u2/(u2 + v2) );								// = x1 + (Px-x1) * pow( 1-D*t , -(u2*dy)/(u2*dy + v2*dx) );
		}

		if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 1" << endl; }
			ye = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 2" << endl; }
			ye = Py * pow( 1-D*t , -v2/(u2 + v2) );								// = y1 + (Py-y1) * pow( 1-D*t , -(v2*dy)/(u2*dy + v2*dx) );
		}

		if(Case->ReachedNegativeCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
			ze = 0;
		}else if(Case->ReachedPositiveCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
			ze = 1;
		}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 3" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = -w1_w2w1 + (w1_w2w1 + Pz) * pow( 1-D*t , -2*(w2-w1)/(u2 + v2) );	// = z1 - ((w1*dz)/(w2-w1)) + ((w1*dz)/(w2-w1) + (Pz-z1)) * pow( 1-D*t , -(2*(w2-w1)*dx*dy)/(dz*(u2*dy + v2*dx)) );
			
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 4" << endl; }
			ze = Pz - ((2*w1)/(u2+v2)) * log( 1 - D*t );								// = Pz - ((2*w1*dx*dy)/(u2*dy+v2*dx)) * log( 1 - D*t );

		}

	}else{

		if(Case->ReachedPositiveABorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 3" << endl; }
			xe = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 4" << endl; }
			xe = Px * exp( u2*_2PxPy*t );											// = x1 + (Px-x1) * exp( (2*u2*(Px-x1)*(Py-y1)/(dx*dx*dy))*t );
		}

		if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 3" << endl; }
			ye = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 4" << endl; }
			ye = Py * exp( v2*_2PxPy*t );											// = y1 + (Py-y1) * exp( (2*v2*(Px-x1)*(Py-y1)/(dx*dy*dy))*t );
		}

		if(Case->ReachedNegativeCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 5" << endl; }
			ze = 0;
		}else if(Case->ReachedPositiveCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 6" << endl; }
			ze = 1;
		}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 7" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = -w1_w2w1 + (w1_w2w1 + Pz) * exp( 2*(w2-w1)*_2PxPy*t );				// = z1 - (w1*dz/(w2-w1)) + ( w1*dz/(w2-w1) + (Pz-z1) ) * exp( (4*(w2-w1)*(Px-x1)*(Py-y1)/(dx*dy*dz))*t );

		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 8" << endl; }
			ze = Pz + 2*w1*_2PxPy*t;												// = Pz + (4*w1*(Px-x1)*(Py-y1)/(dx*dy))*t;

		}

	}

	/*
	double numt=0;
	double numdt=0.00001;
	double numtmax=t;
	
	double numPx=Px;
	double numPy=Py;
	double numPz=Pz;

	while(numt<numtmax){

		double numu = (2*u2)*(numPx-x1)*(numPx-x1)*(numPy-y1);
		double numv = (2*v2)*(numPx-x1)*(numPy-y1)*(numPy-y1);
		double numw = (4*w1)*(numPx-x1)*(numPy-y1) + 4*(w2-w1)*(numPx-x1)*(numPy-y1)*(numPz-z1);

		numPx+=numu*numdt;
		numPy+=numv*numdt;
		numPz+=numw*numdt;

		numt+=numdt;
	}
	

	cout << "New Transformed Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	cout << "New  Numerical  Position: (" << numPx << ", " << numPy << ", " << numPz << ")" << endl;

	Sleep(20000);
	*/

	/*

	if(xe<0||xe>1.00000001||ye<0||ye>1.00000001||ze<0||ze>1.00000001){
		cout << "2 Adjacent Solids Fault - Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
		//Sleep(3000);
	}

	*/

	Case->Pa_a1 = xe;
	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

}

//Advection function used when two neighbouring voxels are solid (15 cases)
void CaseTwoSolids(int i,Coords* C, Voxel* LatticeElement, double* dt){

	int x = C->x;
	int y = C->y;
	int z = C->z;

	double u1 = LatticeElement->u1;
	double u2 = LatticeElement->u2;
	double v1 = LatticeElement->v1;
	double v2 = LatticeElement->v2;
	double w1 = LatticeElement->w1;
	double w2 = LatticeElement->w2;

	double x1 = (double)x;
	double x2 = (double)(x1+1);
	double y1 = (double)y;
	double y2 = (double)(y1+1);
	double z1 = (double)z;
	double z2 = (double)(z1+1);

	Particle* P = &(Particles[i]);

	double Px = P->X;					//X coordinate of particle
	double Py = P->Y;					//Y coordinate of particle
	double Pz = P->Z;					//Z coordinate of particle

	bool ReachesBorderX=false;
	bool ReachesBorderY=false;
	bool ReachesBorderZ=false;

	double xe, ye, ze;

	Case2OppositeSolidsCoords OppCoords;	//Struct for opposite solids function

	OppCoords.dt = dt;
	OppCoords.ReachedPositiveBBorder=false;
	OppCoords.ReachedNegativeBBorder=false;
	OppCoords.ReachedPositiveCBorder=false;
	OppCoords.ReachedNegativeCBorder=false;

	Case2AdjacentSolidsCoords AdjCoords;	//Struct for adjacent solids function

	AdjCoords.dt = dt;
	AdjCoords.ReachedPositiveABorder=false;
	AdjCoords.ReachedPositiveBBorder=false;
	AdjCoords.ReachedPositiveCBorder=false;
	AdjCoords.ReachedNegativeCBorder=false;

	unsigned char Solids = (unsigned char)(LatticeElement->Solid);

	//Opposing solids

	if(Solids&PositiveXSolid && Solids&NegativeXSolid){			//Solid on both x sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Negative X Solid" << endl;
		}

		OppCoords.Velocity_b1 = v1;
		OppCoords.Velocity_b2 = v2;
		OppCoords.Velocity_c1 = w1;
		OppCoords.Velocity_c2 = w2;
		
		OppCoords.Pa_a1 = Px - x1;
		OppCoords.Pb_b1 = Py - y1;
		OppCoords.Pc_c1 = Pz - z1;

		Case2OppositeSolids(&OppCoords);

		if(OppCoords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;			

		}

		if(OppCoords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(OppCoords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		if(OppCoords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}

		xe = OppCoords.Pa_a1 + x1;
		ye = OppCoords.Pb_b1 + y1;
		ze = OppCoords.Pc_c1 + z1;


	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid){		//Solid on both y sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Y Solid, Negative Y Solid" << endl;
		}

		OppCoords.Velocity_b1 = u1;
		OppCoords.Velocity_b2 = u2;
		OppCoords.Velocity_c1 = w1;
		OppCoords.Velocity_c2 = w2;
		
		OppCoords.Pa_a1 = y2 - Py;
		OppCoords.Pb_b1 = Px - x1;
		OppCoords.Pc_c1 = Pz - z1;

		Case2OppositeSolids(&OppCoords);

		if(OppCoords.ReachedPositiveBBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;			

		}

		if(OppCoords.ReachedNegativeBBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}
		
		if(OppCoords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		if(OppCoords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}

		xe = OppCoords.Pb_b1 + x1;
		ye = y2 - OppCoords.Pa_a1;
		ze = OppCoords.Pc_c1 + z1;		

	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid){		//Solid on both z sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid, Negative Z Solid" << endl;
		}

		OppCoords.Velocity_b1 = v1;
		OppCoords.Velocity_b2 = v2;
		OppCoords.Velocity_c1 = u1;
		OppCoords.Velocity_c2 = u2;
		
		OppCoords.Pa_a1 = z2 - Pz;
		OppCoords.Pb_b1 = Py - y1;
		OppCoords.Pc_c1 = Px - x1;

		Case2OppositeSolids(&OppCoords);

		if(OppCoords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;			

		}

		if(OppCoords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(OppCoords.ReachedPositiveCBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		if(OppCoords.ReachedNegativeCBorder){
			
			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
			
		}

		xe = OppCoords.Pc_c1 + x1;
		ye = OppCoords.Pb_b1 + y1;
		ze = z2 - OppCoords.Pa_a1;


	//Adjacent solids

	}else if(Solids&NegativeXSolid && Solids&NegativeYSolid){		//(1) Solid on x1 and y1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Negative Y Solid" << endl;
		}

		AdjCoords.Velocity_a2 = u2;
		AdjCoords.Velocity_b2 = v2;
		AdjCoords.Velocity_c1 = w1;
		AdjCoords.Velocity_c2 = w2;
		
		AdjCoords.Pa_a1 = Px - x1;
		AdjCoords.Pb_b1 = Py - y1;
		AdjCoords.Pc_c1 = Pz - z1;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}

		xe = AdjCoords.Pa_a1 + x1;
		ye = AdjCoords.Pb_b1 + y1;
		ze = AdjCoords.Pc_c1 + z1;
		

	}else if(Solids&NegativeXSolid && Solids&NegativeZSolid){		//(2) Solid on x1 and z1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Negative Z Solid" << endl;
		}
		
		AdjCoords.Velocity_a2 = u2;
		AdjCoords.Velocity_b2 = w2;
		AdjCoords.Velocity_c1 = -v2;
		AdjCoords.Velocity_c2 = -v1;
		
		AdjCoords.Pa_a1 = Px - x1;
		AdjCoords.Pb_b1 = Pz - z1;
		AdjCoords.Pc_c1 = y2 - Py;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}

		xe = AdjCoords.Pa_a1 + x1;
		ye = y2 - AdjCoords.Pc_c1;
		ze = AdjCoords.Pb_b1 + z1;


	}else if(Solids&NegativeXSolid && Solids&PositiveYSolid){		//(3) Solid on x1 and y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive Y Solid" << endl;
		}

		AdjCoords.Velocity_a2 = u2;
		AdjCoords.Velocity_b2 = -v1;
		AdjCoords.Velocity_c1 = -w2;
		AdjCoords.Velocity_c2 = -w1;
		
		AdjCoords.Pa_a1 = Px - x1;
		AdjCoords.Pb_b1 = y2 - Py;
		AdjCoords.Pc_c1 = z2 - Pz;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}

		xe = AdjCoords.Pa_a1 + x1;
		ye = y2 - AdjCoords.Pb_b1;
		ze = z2 - AdjCoords.Pc_c1;

	}else if(Solids&NegativeXSolid && Solids&PositiveZSolid){		//(4) Solid on x1 and z2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive Z Solid" << endl;
		}

		AdjCoords.Velocity_a2 = u2;
		AdjCoords.Velocity_b2 = -w1;
		AdjCoords.Velocity_c1 = v1;
		AdjCoords.Velocity_c2 = v2;
		
		AdjCoords.Pa_a1 = Px - x1;
		AdjCoords.Pb_b1 = z2 - Pz;
		AdjCoords.Pc_c1 = Py - y1;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
			
		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		xe = AdjCoords.Pa_a1 + x1;
		ye = AdjCoords.Pc_c1 + y1;
		ze = z2 - AdjCoords.Pb_b1;

	}else if(Solids&NegativeYSolid && Solids&NegativeZSolid){		//(5) Solid on y1 and z1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid, Negative Z Solid" << endl;
		}

		AdjCoords.Velocity_a2 = w2;
		AdjCoords.Velocity_b2 = v2;
		AdjCoords.Velocity_c1 = -u2;
		AdjCoords.Velocity_c2 = -u1;
		
		AdjCoords.Pa_a1 = Pz - z1;
		AdjCoords.Pb_b1 = Py - y1;
		AdjCoords.Pc_c1 = x2 - Px;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;

			ReachesBorderZ=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
			
		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		xe = x2 - AdjCoords.Pc_c1;
		ye = AdjCoords.Pb_b1 + y1;
		ze = AdjCoords.Pa_a1 + z1;

	}else if(Solids&NegativeYSolid && Solids&PositiveXSolid){		//(6) Solid on y1 and x2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid, Positive X Solid" << endl;
		}

		AdjCoords.Velocity_a2 = -u1;
		AdjCoords.Velocity_b2 = v2;
		AdjCoords.Velocity_c1 = -w2;
		AdjCoords.Velocity_c2 = -w1;
		
		AdjCoords.Pa_a1 = x2 - Px;
		AdjCoords.Pb_b1 = Py - y1;
		AdjCoords.Pc_c1 = z2 - Pz;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
			
		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;

			ReachesBorderZ=true;

		}

		xe = x2 - AdjCoords.Pa_a1;
		ye = AdjCoords.Pb_b1 + y1;
		ze = z2 - AdjCoords.Pc_c1;

	}else if(Solids&NegativeYSolid && Solids&PositiveZSolid){		//(7) Solid on y1 and z2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid, Positive Z Solid" << endl;
		}

		AdjCoords.Velocity_a2 = -w1;
		AdjCoords.Velocity_b2 = v2;
		AdjCoords.Velocity_c1 = u1;
		AdjCoords.Velocity_c2 = u2;
		
		AdjCoords.Pa_a1 = z2 - Pz;
		AdjCoords.Pb_b1 = Py - y1;
		AdjCoords.Pc_c1 = Px - x1;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
			
		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		xe = AdjCoords.Pc_c1 + x1;
		ye = AdjCoords.Pb_b1 + y1;
		ze = z2 - AdjCoords.Pa_a1;

	}else if(Solids&NegativeZSolid && Solids&PositiveXSolid){		//(8) Solid on z1 and x2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid, Positive X Solid" << endl;
		}

		AdjCoords.Velocity_a2 = -u1;
		AdjCoords.Velocity_b2 = w2;
		AdjCoords.Velocity_c1 = v1;
		AdjCoords.Velocity_c2 = v2;
		
		AdjCoords.Pa_a1 = x2 - Px;
		AdjCoords.Pb_b1 = Pz - z1;
		AdjCoords.Pc_c1 = Py - y1;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;			

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		xe = x2 - AdjCoords.Pa_a1;
		ye = AdjCoords.Pc_c1 + y1;
		ze = AdjCoords.Pb_b1 + z1;

	}else if(Solids&NegativeZSolid && Solids&PositiveYSolid){		//(9) Solid on z1 and y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid, Positive Y Solid" << endl;
		}

		AdjCoords.Velocity_a2 = w2;
		AdjCoords.Velocity_b2 = -v1;
		AdjCoords.Velocity_c1 = u1;
		AdjCoords.Velocity_c2 = u2;
		
		AdjCoords.Pa_a1 = Pz - z1;
		AdjCoords.Pb_b1 = y2 - Py;
		AdjCoords.Pc_c1 = Px - x1;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		xe = AdjCoords.Pc_c1 + x1;
		ye = y2 - AdjCoords.Pb_b1;
		ze = AdjCoords.Pa_a1 + z1;

	}else if(Solids&PositiveXSolid && Solids&PositiveYSolid){		//(10) Solid on x2 and y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Positive Y Solid" << endl;
		}

	//	cout << "Positive x2, y2" << endl;

		AdjCoords.Velocity_a2 = -u1;
		AdjCoords.Velocity_b2 = -v1;
		AdjCoords.Velocity_c1 = w1;
		AdjCoords.Velocity_c2 = w2;

	//	cout << "Velocities: a2 = " << -u1 << "; b2 = " << -v1 << "; c1 = " << w1 << "; c2 = " << w2 << ";" << endl;
		
		AdjCoords.Pa_a1 = x2 - Px;
		AdjCoords.Pb_b1 = y2 - Py;
		AdjCoords.Pc_c1 = Pz - z1;

		//cout << "Transformed Position: (" << AdjCoords.Pa_a1 << ", " << AdjCoords.Pb_b1 << ", " << AdjCoords.Pc_c1 << ")" << endl;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		xe = x2 - AdjCoords.Pa_a1;
		ye = y2 - AdjCoords.Pb_b1;
		ze = AdjCoords.Pc_c1 + z1;

	}else if(Solids&PositiveZSolid && Solids&PositiveYSolid){		//(11) Solid on z2 and y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid, Positive Y Solid" << endl;
		}

		AdjCoords.Velocity_a2 = -v1;
		AdjCoords.Velocity_b2 = -w1;
		AdjCoords.Velocity_c1 = u1;
		AdjCoords.Velocity_c2 = u2;
		
		AdjCoords.Pa_a1 = y2 - Py;
		AdjCoords.Pb_b1 = z2 - Pz;
		AdjCoords.Pc_c1 = Px - x1;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		if(AdjCoords.ReachedNegativeCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		xe = AdjCoords.Pc_c1 + x1;
		ye = y2 - AdjCoords.Pa_a1;
		ze = z2 - AdjCoords.Pb_b1;

	}else if(Solids&PositiveZSolid && Solids&PositiveXSolid){		//(12) Solid on z2 and x2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid, Positive X Solid" << endl;
		}

		AdjCoords.Velocity_a2 = -w1;
		AdjCoords.Velocity_b2 = -u1;
		AdjCoords.Velocity_c1 = v1;
		AdjCoords.Velocity_c2 = v2;
		
		AdjCoords.Pa_a1 = z2 - Pz;
		AdjCoords.Pb_b1 = x2 - Px;
		AdjCoords.Pc_c1 = Py - y1;

		Case2AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedPositiveABorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		if(AdjCoords.ReachedPositiveBBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}
		
		if(AdjCoords.ReachedPositiveCBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}

		if(AdjCoords.ReachedNegativeCBorder){
			
			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		xe = x2 - AdjCoords.Pb_b1;
		ye = AdjCoords.Pc_c1 + y1;
		ze = z2 - AdjCoords.Pa_a1;

	}

	
	if(!ReachesBorderX){		//If reached boundary, new position has already been set
		P->X = xe;
	}
	if(!ReachesBorderY){
		P->Y = ye;
	}
	if(!ReachesBorderZ){
		P->Z = ze;
	}

}

struct Case3AdjacentSolidsCoords{
double Velocity_a1;				//Inward velocity opposite solid a direction
double Velocity_b1;				//Inward velocity opposite solid b direction
double Velocity_c2;				//Outward velocity opposite solid c direction

double Pa_a1;					//Position in a direction within voxel (Pa-a1)
double Pb_b1;					//Position in b direction within voxel (Pb-b1)
double Pc_c1;					//Position in c direction within voxel (Pc-c1)

bool ReachedNegativeABorder;	//Function returns whether particle reached border
bool ReachedNegativeBBorder;
bool ReachedPositiveCBorder;

double* dt;						//Advection time pointer
};


//Advects particles in a voxel with three adjacent solids on the a2, b2 and c1 faces, used as a general case
//Dimensions referred to as x,y,z in the function, but correspond to a,b,c in the struct to reduce
//  confusion when transforming into the general case in Case3Solids()
void Case3AdjacentSolids(Case3AdjacentSolidsCoords* Case){

	//cout<<"3 Adjacent Solids"<<endl;

	double u1 = Case->Velocity_a1;
	double v1 = Case->Velocity_b1;
	double w2 = Case->Velocity_c2;

	double Px = Case->Pa_a1;
	double Py = Case->Pb_b1;
	double Pz = Case->Pc_c1;

	double x1 = 0;	//Imaginary voxel at 0,0,0
	double x2 = 1;
	double y1 = 0;
	double y2 = 1;
	double z1 = 0;
	double z2 = 1;

	double* dt = Case->dt;

	double _4PxPyPz = 4*(x2-Px)*(y2-Py)*(Pz-z1);

	double u = u1*(1-Px)*_4PxPyPz;			// = (4*u1/(dx*dx*dy*dz))*(x2-Px)*(x2-Px)*(y2-Py)*(Pz-z1);
	double v = v1*(1-Py)*_4PxPyPz;			// = (4*v1/(dx*dy*dy*dz))*(x2-Px)*(y2-Py)*(y2-Py)*(Pz-z1);
	double w = w2*Pz*_4PxPyPz;				// = (4*w2/(dx*dy*dz*dz))*(x2-Px)*(y2-Py)*(Pz-z1)*(Pz-z1);

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 3 Adjacent Solids:" << endl;
		AdvectDebug << "u1 = " << u1 << endl << "v1 = " << v1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	//cout << "Velocity: (" << u << ", " << v << ", " << w << ")" <<endl;
	//cout<<"Face Velocities: (" << u1 << ", " << v1 << ", " << w2 << ")" <<endl;

	double Tx, Ty, Tz;

	//Calculate time to reach the border, Tau
	
	//Tau x

	if(u<0){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 1" << endl; }
		Tx = log( 1-Px )/(_4PxPyPz*u1);		// = (dx*dx*dy*dz*log((x2-Px)/dx)/(4*(x2-Px)*(y2-Py)*(Pz-z1)*u1));
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 2" << endl; }
		Tx = DBL_MAX;
	}

	//Tau y

	if(v<0){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 1" << endl; }
		Ty = log( 1-Py )/(_4PxPyPz*v1);		// = (dy*dy*dx*dz*log((y2-Py)/dy)/(4*(x2-Px)*(y2-Py)*(Pz-z1)*v1));
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 2" << endl; }
		Ty = DBL_MAX;
	}

	//Tau z

	if(w>0){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
		Tz = -log( Pz )/(_4PxPyPz*w2);		// = (dz*dz*dx*dy*log(dz/(Pz-z1))/(4*(x2-Px)*(y2-Py)*(Pz-z1)*w2));
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
		Tz = DBL_MAX;
	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}
	//cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	double t;

	if(Tx<=Ty && Tx<=Tz){			//Reaches x boundary first
		//cout<<"Min tau x"<<endl;

		if(Tx>(*dt)){							//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Tx;							//Advect for tau_x
			(*dt) -= Tx;						//Reduce dt by tau_x time used

			Case->ReachedNegativeABorder=true;
		//	cout<<"Reaches negative A border"<<endl;

		}

	}else if(Ty<=Tx && Ty<=Tz){		//Reaches y boundary first
		//cout<<"Min tau y"<<endl;

		if(Ty>(*dt)){							//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Ty;							//Advect for tau_y
			(*dt) -= Ty;						//Reduce dt by tau_y time used

			Case->ReachedNegativeBBorder=true;
			//cout<<"Reaches negative B border"<<endl;

		}

	}else if(Tz<=Tx && Tz<=Ty){		//Reaches z boundary first
		//cout<<"Min tau z"<<endl;

		if(Tz>(*dt)){						//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			//cout<<"Tz > dt"<<endl;
			//cout<<(*dt)<<endl;
			
		}else{									//Advect to the border
			
			t = Tz;							//Advect for tau_z
			(*dt) -= Tz;						//Reduce dt by tau_z time used

			Case->ReachedPositiveCBorder=true;
		//	cout<<"Reaches positive C border"<<endl;
			
		}

	}

	//cout << "Analytical Tau:      (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	//Positions

	double xe, ye, ze;
	if(Case->ReachedNegativeABorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 1" << endl; }
		xe = 0;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 2" << endl; }
		xe = 1 - (1-Px) * exp(-u1*_4PxPyPz*t);		// = x2-(x2-Px)*exp(-4*u1*(x2-Px)*(y2-Py)*(Pz-z1)*t/(dx*dx*dy*dz));
	}

	if(Case->ReachedNegativeBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 1" << endl; }
		ye = 0;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 2" << endl; }
		ye = 1 - (1-Py) * exp(-v1*_4PxPyPz*t);		// = y2-(y2-Py)*exp(-4*v1*(x2-Px)*(y2-Py)*(Pz-z1)*t/(dx*dy*dy*dz));
	}

	if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
		ze = 1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
		ze = Pz * exp(w2*_4PxPyPz*t);			// = z1+(Pz-z1)*exp(4*w2*(x2-Px)*(y2-Py)*(Pz-z1)*t/(dx*dy*dz*dz));
	}

	/*
	double numt=0;
	double numdt=0.00001;
	double numtmax=t;
	
	double numPx=Px;
	double numPy=Py;
	double numPz=Pz;

	while(numt<numtmax){

		double numu = (2*u2)*(numPx-x1)*(numPx-x1)*(numPy-y1)/(dx*dx*dy);
		double numv = (2*v2)*(numPx-x1)*(numPy-y1)*(numPy-y1)/(dx*dy*dy);
		double numw = (4*w1)*(numPx-x1)*(numPy-y1)/(dx*dy) + 4*(w2-w1)*(numPx-x1)*(numPy-y1)*(numPz-z1)/(dx*dy*dz);

		numPx+=numu*numdt;
		numPy+=numv*numdt;
		numPz+=numw*numdt;

		numt+=numdt;
	}
	*/

	//cout << "New Transformed Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	//cout << "New  Numerical  Position: (" << numPx << ", " << numPy << ", " << numPz << ")" << endl;

	/*
	if(xe<0||xe>1.00000001||ye<0||ye>1.00000001||ze<0||ze>1.00000001){
		cout << "3 Adjacent Solids Fault - Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
		//Sleep(3000);
	}
	*/

	Case->Pa_a1 = xe;
	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

	
	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

}

struct Case3OppositeSolidsCoords{
double Velocity_b1;				//Inward velocity in b direction, free of solids
double Velocity_b2;				//Outward velocity in b direction, free of solids
double Velocity_c2;				//Outward velocity opposite solid c direction

double Pa_a1;					//Position in a direction within voxel (Pa-a1)
double Pb_b1;					//Position in b direction within voxel (Pb-b1)
double Pc_c1;					//Position in c direction within voxel (Pc-c1)

bool ReachedNegativeBBorder;	//Function returns whether particle reached border
bool ReachedPositiveBBorder;
bool ReachedPositiveCBorder;

double* dt;						//Advection time pointer
};


//Advects particles in a voxel with three opposite solids on the a1, a2 and c1 faces, used as a general case
//Dimensions referred to as x,y,z in the function, but correspond to a,b,c in the struct to reduce
//  confusion when transforming into the general case in Case3Solids()
void Case3OppositeSolids(Case3OppositeSolidsCoords* Case){

	//cout<<"3 Opposite Solids"<<endl;

	double v1 = Case->Velocity_b1;
	double v2 = Case->Velocity_b2;
	double w2 = Case->Velocity_c2;

	//cout<<"Velocity: "<<Case->Velocity_b1<<'\t'<<Case->Velocity_b2<<endl;

	double Px = Case->Pa_a1;
	double Py = Case->Pb_b1;
	double Pz = Case->Pc_c1;

	double x1 = 0;	//Imaginary voxel at 0,0,0
	double x2 = 1;
	double y1 = 0;
	double y2 = 1;
	double z1 = 0;
	double z2 = 1;

	double* dt = Case->dt;

	double _6PxPxPz = 6*(1-Px)*Px*Pz;

	double v_face = (v1 + (v2-v1)*Py);

	double u = 0;
	double v = 2*_6PxPxPz*v_face;		// = (12*v1/(dx*dx*dz))*(x2-Px)*(Px-x1)*(Pz-z1)+(12*(v2-v1)/(dx*dx*dy*dz))*(x2-Px)*(Px-x1)*(Pz-z1)*(Py-y1);
	double w = w2*Pz*_6PxPxPz;						// = (6*w2/(dx*dx*dz*dz))*(x2-Px)*(Px-x1)*(Pz-z1)*(Pz-z1);

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 3 Opposite Solids:" << endl;
		AdvectDebug << "v1 = " << v1 << endl << "v2 = " << v2 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	//cout << "Velocity: (" << u << ", " << v << ", " << w << ")" <<endl;

	double Tx;
	double Ty;
	double Tz;

	//Calculate time to reach the border, Tau
	
	//Tau x

		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tx 1" << endl; }
		Tx = DBL_MAX;

	//Tau y

	if(v==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 0" << endl; }
				Ty = DBL_MAX;
	}else if(v>0){
		if(v2<=0 || (v2<0 && v1>0)){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 1" << endl; }
				Ty = DBL_MAX;
		}else if(fabs(v2-v1)>Min_dV){
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 2" << endl; }
				Ty = (1-Py)/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 3" << endl; }
				Ty = (1 - pow( v2/v_face , -w2/(2*(v2-v1)) ))/(w2*_6PxPxPz);	// = (dx*dx*dz*dz/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1))) * (1 - pow( (v1*dy+(v2-v1)*dy)/(v1*dy+(v2-v1)*(Py-y1)) , -w2*dy/(2*(v2-v1)*dz) ));
			}
		}else{
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 4" << endl; }
				Ty = (1-Py)/v ;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 5" << endl; }
				Ty = (1 - exp( -w2*(1-Py)/(2*v_face) ))/(w2*_6PxPxPz);						// = ((dx*dx*dz*dz)/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1))) * (1 - exp( -w2*(y2-Py)/(2*v1*dz) ));
			}
		}
	}else{
		if(v1>=0 || (v2<0 && v1>0)){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 6" << endl; }
				Ty = DBL_MAX;
		}else if(fabs(v2-v1)>Min_dV){
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 7" << endl; }
				Ty = -Py/v ;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 8" << endl; }
				Ty = (1 - pow( v1/v_face , -w2/(2*(v2-v1)) ))/(w2*_6PxPxPz);	// = (dx*dx*dz*dz/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1)))*(1-pow((v1*dy)/(v1*dy+(v2-v1)*(Py-y1)), -w2*dy/(2*(v2-v1)*dz)));
			}
		}else{
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 9" << endl; }
				Ty = -Py/v ;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 10" << endl; }
				Ty = (1 - exp( w2*Py/(2*v_face) ))/(w2*_6PxPxPz);							// = ((dx*dx*dz*dz)/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1)))*(1-exp(-w2*(y1-Py)/(2*v1*dz)));
			}
		}
	}

	//Tau z

	if(w>0){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
		Tz = ( 1/(Pz-z1) - 1 )/( 6*(1-Px)*Px*w2 );										// = (dx*dx*dz*dz/(6*(x2-Px)*(Px-x1)*w2))*(1/(Pz-z1)-1/dz);
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
		Tz = DBL_MAX;
	}


	//cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" <<endl;

	double t;

	if(Ty<=Tx && Ty<=Tz){		//Reaches y boundary first
	//	cout<<"Min tau = y"<<endl;

		if(Ty>(*dt)){							//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Ty;							//Advect for tau_y
			(*dt) -= Ty;						//Reduce dt by tau_y time used

			if(v<0){							//Reaches border with -ve neighbour
				Case->ReachedNegativeBBorder=true;
			}else{
				Case->ReachedPositiveBBorder=true;
			}

		}

	}else if(Tz<=Tx && Tz<=Ty){		//Reaches z boundary first
		//cout<<"Min tau = z"<<endl;

		if(Tz>(*dt)){						//Remains in same voxel

			t = (*dt);						//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Tz;							//Advect for tau_z
			(*dt) -= Tz;						//Reduce dt by tau_z time used

			Case->ReachedPositiveCBorder=true;
			
		}

	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

	//cout << "Analytical Tau:      (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	//Positions

	double xe, ye, ze;

	if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "xe 1" << endl; }
	xe = Px;

	if(Case->ReachedNegativeBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 1" << endl; }
			ye = 0;
	}else if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 2" << endl; }
			ye = 1;
	}else if(fabs(v2-v1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 3" << endl; }
			double v1_v2v1 = v1/(v2-v1);
			ye = -v1_v2v1 + (v1_v2v1 + Py) * pow( 1-w2*_6PxPxPz*t , -2*(v2-v1)/w2 );		// = y1 - (v1*dy/(v2-v1)) + pow( 1-6*w2*(x2-Px)*(Px-x1)*(Pz-z1)*t/(dx*dx*dz*dz) , -2*(v2-v1)*dz/(w2*dy))*(v1*dy/(v2-v1)+(Py-y1));
	}else{
		if(w2==0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 4" << endl; }
			ye = Py + v*t;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 5" << endl; }
			ye = Py - (2*v1/w2) * log( 1 - w2*_6PxPxPz*t );											// = Py - (2*v1*dz/w2) * log( 1 - (6*w2*(x2-Px)*(Px-x1)*(Pz-z1)*t/(dx*dx*dz*dz)) );
		}
	}

	if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
		ze = 1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
		ze = Pz/(1 - w2*_6PxPxPz*t);																// = z1 + dx*dx*dz*dz*(Pz-z1)/(dx*dx*dz*dz-6*w2*(x2-Px)*(Px-x1)*(Pz-z1)*t);
	}

	//cout << "New Transformed Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;

	/*
	if(xe<0||xe>1.00000001||ye<0||ye>1.00000001||ze<0||ze>1.00000001){
		cout << "3 Opposite Solids Fault - Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
		//Sleep(3000);
	}
	*/

	Case->Pa_a1 = xe;
	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

}

//Advection function used when three neighbouring voxels are solid (20 cases)
void CaseThreeSolids(int i,Coords* C, Voxel* LatticeElement, double* dt){

	int x = C->x;
	int y = C->y;
	int z = C->z;

	double u1 = LatticeElement->u1;
	double u2 = LatticeElement->u2;
	double v1 = LatticeElement->v1;
	double v2 = LatticeElement->v2;
	double w1 = LatticeElement->w1;
	double w2 = LatticeElement->w2;

	double x1 = (double)x;
	double x2 = (double)(x1+1);
	double y1 = (double)y;
	double y2 = (double)(y1+1);
	double z1 = (double)z;
	double z2 = (double)(z1+1);

	Particle* P = &(Particles[i]);

	double Px = P->X;					//X coordinate of particle
	double Py = P->Y;					//Y coordinate of particle
	double Pz = P->Z;					//Z coordinate of particle

	bool ReachesBorderX=false;
	bool ReachesBorderY=false;
	bool ReachesBorderZ=false;

	double xe, ye, ze;			//New positions after advection

	Case3AdjacentSolidsCoords coords;	//Struct for adjacent solids function
	Case3OppositeSolidsCoords oppcoords;	//Struct for opposite solids function

	coords.dt = dt;
	coords.ReachedNegativeABorder=false;
	coords.ReachedNegativeBBorder=false;
	coords.ReachedPositiveCBorder=false;

	oppcoords.dt = dt;
	oppcoords.ReachedNegativeBBorder=false;
	oppcoords.ReachedPositiveBBorder=false;
	oppcoords.ReachedPositiveCBorder=false;

	unsigned char Solids = (unsigned char)(LatticeElement->Solid);
	
	if(Solids&PositiveXSolid && Solids&PositiveYSolid && Solids&NegativeZSolid){		//(1) Solid on x2, y2 and z1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Positive Y Solid, Negative Z Solid" << endl; 
		}

		//a=x, b=y, c=z
		coords.Velocity_a1 = u1;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = v1;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = w2;//Outward velocity opposite solid c direction
		
		coords.Pa_a1 = Px - x1;
		coords.Pb_b1 = Py - y1;
		coords.Pc_c1 = Pz - z1;


		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;


			ReachesBorderX=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}


		xe = coords.Pa_a1 + x1;
		ye = coords.Pb_b1 + y1;
		ze = coords.Pc_c1 + z1;
		
		

	}else if(Solids&PositiveXSolid && Solids&PositiveYSolid && Solids&PositiveZSolid){		//Solid on x2,y2 &z2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Positive Y Solid, Positive Z Solid" << endl; 
		}

		//a=z, b=y, c=-x
		coords.Velocity_a1 = w1;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = v1;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = -u1;//Outward velocity opposite solid c direction

	//	cout<<"Velocities before transform: u1 = "<<u1<<" , v1 = "<<v1<<" , w1 = "<<w1<<endl;
		
		coords.Pa_a1 = Pz - z1;
		coords.Pb_b1 = Py - y1;
		coords.Pc_c1 = x2 - Px;

		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;

			ReachesBorderZ=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}


		xe = x2 - coords.Pc_c1;
		ye = coords.Pb_b1 + y1;
		ze = coords.Pa_a1 + z1;

	}else if(Solids&NegativeXSolid && Solids&PositiveYSolid && Solids&PositiveZSolid){//Solid on x1,y2 &z2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive Y Solid, Positive Z Solid" << endl; 
		}

		//a=-x, b=y, c=-z

		coords.Velocity_a1 = -u2;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = v1;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = -w1;//Outward velocity opposite solid c direction
		
		coords.Pa_a1 = x2 - Px;
		coords.Pb_b1 = Py - y1;
		coords.Pc_c1 = z2 - Pz;


		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = x2 - coords.Pa_a1;
		ye = coords.Pb_b1 + y1;
		ze = z2 - coords.Pc_c1;
		
	}else if(Solids&PositiveXSolid && Solids&NegativeYSolid && Solids&PositiveZSolid){//Solid on x2,y1 &z2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Negative Y Solid, Positive Z Solid" << endl; 
		}

		//a=x, b=-y, c=-z

		coords.Velocity_a1 = u1;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = -v2;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = -w1;//Outward velocity opposite solid c direction
		
		coords.Pa_a1 = Px - x1;
		coords.Pb_b1 = y2 - Py;
		coords.Pc_c1 = z2 - Pz;


		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;

			ReachesBorderX=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = coords.Pa_a1 + x1;
		ye = y2 - coords.Pb_b1;
		ze = z2 - coords.Pc_c1;
		

	}else if(Solids&NegativeXSolid &&Solids&NegativeYSolid && Solids&NegativeZSolid){//Solid on x1,y1 &z1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Negative Y Solid, Negative Z Solid" << endl; 
		}

		//a=-x, b=-y, c=z

		coords.Velocity_a1 = -u2;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = -v2;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = w2;//Outward velocity opposite solid c direction
		
		coords.Pa_a1 = x2 - Px;
		coords.Pb_b1 = y2 - Py;
		coords.Pc_c1 = Pz - z1;

		//cout<<"Coords: "<<coords.Pa_a1<<'\t'<<coords.Pb_b1<<'\t'<<coords.Pc_c1<<endl;



		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}


		xe = x2 - coords.Pa_a1;
		ye = y2 - coords.Pb_b1;
		ze = coords.Pc_c1 + z1;

	}else if(Solids&NegativeXSolid && Solids&NegativeYSolid && Solids&PositiveZSolid){//Solid on x1,y1 &z2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Negative Y Solid, Positive Z Solid" << endl; 
		}

		//a=-x, b=z, c=y

		coords.Velocity_a1 = -u2;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = w1;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = v2;//Outward velocity opposite solid c direction
		
		coords.Pa_a1 = x2 - Px;
		coords.Pb_b1 = Pz - z1;
		coords.Pc_c1 = Py - y1;


		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}


		xe = x2 - coords.Pa_a1;
		ye = coords.Pc_c1 + y1;
		ze = coords.Pb_b1 + z1;

	}else if(Solids&NegativeXSolid &&Solids&PositiveYSolid && Solids&NegativeZSolid){//Solid on x1,y2 &z1 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive Y Solid, Negative Z Solid" << endl; 
		}

		//a=-z, b=y, c=x
		coords.Velocity_a1 = -w2;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = v1;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = u2;//Outward velocity opposite solid c direction
		
		coords.Pa_a1 = z2 - Pz;
		coords.Pb_b1 = Py - y1;
		coords.Pc_c1 = Px - x1;

		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;

			ReachesBorderZ=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
		}


		xe = coords.Pc_c1 + x1;
		ye = coords.Pb_b1 + y1;
		ze = z2 - coords.Pa_a1;

		

	}else if(Solids&PositiveXSolid && Solids&NegativeYSolid && Solids&NegativeZSolid){//Solid on x2,y1 &z1 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Negative Y Solid, Negative Z Solid" << endl; 
		}
		
		//a=x, b=-z, c=y

		coords.Velocity_a1 = u1;//Inward velocity opposite solid a direction
		coords.Velocity_b1 = -w2;//Inward velocity opposite solid b direction
		coords.Velocity_c2 = v2;//Outward velocity opposite solid c direction
		
		coords.Pa_a1 = Px - x1;
		coords.Pb_b1 = z2 - Pz;
		coords.Pc_c1 = Py - y1;


		Case3AdjacentSolids(&coords);

		if(coords.ReachedNegativeABorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;

			ReachesBorderX=true;			

		}

		if(coords.ReachedNegativeBBorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}
		
		if(coords.ReachedPositiveCBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}


		xe = coords.Pa_a1 + x1;
		ye = coords.Pc_c1 + y1;
		ze = z2 - coords.Pb_b1;



		//Opposite Solids
	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&NegativeZSolid){//Solid on x2,x1 &z1 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Negative X Solid, Negative Z Solid" << endl; 
		}

		//a=x, b=y, c=z
		oppcoords.Velocity_b1 = v1;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = v2;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = w2;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Px - x1;
		oppcoords.Pb_b1 = Py - y1;
		oppcoords.Pc_c1 = Pz - z1;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;

			ReachesBorderY=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}


		xe = oppcoords.Pa_a1 + x1;
		ye = oppcoords.Pb_b1 + y1;
		ze = oppcoords.Pc_c1 + z1;
		

		

	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&NegativeZSolid){//Solid on y2,y1 &z1 sides
	
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Y Solid, Negative Y Solid, Negative Z Solid" << endl; 
		}
		
		//a=y, b=-x, c=z
		oppcoords.Velocity_b1 = -u2;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = -u1;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = w2;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Py - y1;
		oppcoords.Pb_b1 = x2 - Px;
		oppcoords.Pc_c1 = Pz - z1;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}


		xe = x2 - oppcoords.Pb_b1;
		ye = oppcoords.Pa_a1 + y1;
		ze = oppcoords.Pc_c1 + z1;
		

	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&NegativeYSolid){//Solid on x2,x1 &y1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Negative X Solid, Negative Y Solid" << endl; 
		}

		//a=x, b=-z, c=y
		oppcoords.Velocity_b1 = -w2;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = -w1;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = v2;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Px - x1;
		oppcoords.Pb_b1 = z2 - Pz;
		oppcoords.Pc_c1 = Py - y1;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;

			ReachesBorderZ=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}


		xe = oppcoords.Pa_a1 + x1;
		ye = oppcoords.Pc_c1 + y1;
		ze = z2 - oppcoords.Pb_b1;

		

	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&NegativeXSolid){//Solid on y2,y1 &x1 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Y Solid, Negative Y Solid, Negative X Solid" << endl; 
		}

		//a=y, b=z, c=x
		oppcoords.Velocity_b1 = w1;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = w2;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = u2;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Py - y1;
		oppcoords.Pb_b1 = Pz - z1;
		oppcoords.Pc_c1 = Px - x1;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;

			ReachesBorderZ=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
		}


		xe = oppcoords.Pc_c1 + x1;
		ye = oppcoords.Pa_a1 + y1;
		ze = oppcoords.Pb_b1 + z1;
		
		

	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&NegativeXSolid){//Solid on z2,z1 &x1 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid, Negative Z Solid, Negative X Solid" << endl; 
		}

		//a=z, b=-y, c=x
		oppcoords.Velocity_b1 = -v2;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = -v1;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = u2;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Pz - z1;
		oppcoords.Pb_b1 = y2 - Py;
		oppcoords.Pc_c1 = Px - x1;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;

			ReachesBorderY=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}


		xe = oppcoords.Pc_c1 + x1;
		ye = y2 - oppcoords.Pb_b1;
		ze = oppcoords.Pa_a1 + z1;
		

	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&NegativeYSolid){//Solid on z2,z1 &y1 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid, Negative Z Solid, Negative Y Solid" << endl; 
		}

		//a=z, b=x, c=y
		oppcoords.Velocity_b1 = u1;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = u2;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = v2;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Pz - z1;
		oppcoords.Pb_b1 = Px - x1;
		oppcoords.Pc_c1 = Py - y1;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;

			ReachesBorderX=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}


		xe = oppcoords.Pb_b1 + x1;
		ye = oppcoords.Pc_c1 + y1;
		ze = oppcoords.Pa_a1 + z1;
		
	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveZSolid){//Solid on x2,x1 &z2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Negative X Solid, Positive Z Solid" << endl; 
		}

		//a=x, b=-y, c=-z
		oppcoords.Velocity_b1 = -v2;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = -v1;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = -w1;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Px - x1;
		oppcoords.Pb_b1 = y2 - Py;
		oppcoords.Pc_c1 = z2 - Pz;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;

			ReachesBorderY=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = oppcoords.Pa_a1 + x1;
		ye = y2 - oppcoords.Pb_b1;
		ze = z2 - oppcoords.Pc_c1;

	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&PositiveZSolid){//Solid on y2,y1 &z2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Y Solid, Negative Y Solid, Positive Z Solid" << endl; 
		}

		//a=y, b=x, c=-z
		oppcoords.Velocity_b1 = u1;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = u2;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = -w1;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Py - y1;
		oppcoords.Pb_b1 = Px - x1;
		oppcoords.Pc_c1 = z2 - Pz;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;

			ReachesBorderX=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}


		xe = oppcoords.Pb_b1 + x1;
		ye = oppcoords.Pa_a1 + y1;
		ze = z2 - oppcoords.Pc_c1;


	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveYSolid){//Solid on x2,x1 &y2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive X Solid, Negative X Solid, Positive Z Solid" << endl; 
		}

		//a=x, b=z, c=-y
		oppcoords.Velocity_b1 = w1;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = w2;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = -v1;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Px - x1;
		oppcoords.Pb_b1 = Pz - z1;
		oppcoords.Pc_c1 = y2 - Py;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;

			ReachesBorderZ=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}


		xe = oppcoords.Pa_a1 + x1;
		ye = y2 - oppcoords.Pc_c1;
		ze = oppcoords.Pb_b1 + z1;


	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&PositiveXSolid){//Solid on y2,y1 &x2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Y Solid, Negative Y Solid, Positive X Solid" << endl; 
		}

		//a=y, b=-z, c=-x
		oppcoords.Velocity_b1 = -w2;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = -w1;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = -u1;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Py - y1;
		oppcoords.Pb_b1 = z2 - Pz;
		oppcoords.Pc_c1 = x2 - Px;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Z = (double)(C->zp);
			P->VoxelZ = C->zp;

			ReachesBorderZ=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}


		xe = x2 - oppcoords.Pc_c1;
		ye = oppcoords.Pa_a1 + y1;
		ze = z2 - oppcoords.Pb_b1;


	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&PositiveXSolid){	//Solid on z2,z1 &x2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid, Negative Z Solid, Positive X Solid" << endl; 
		}

		//a=z, b=y, c=-x

		oppcoords.Velocity_b1 = v1;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = v2;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = -u1;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Pz - z1;
		oppcoords.Pb_b1 = Py - y1;
		oppcoords.Pc_c1 = x2 - Px;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;

			ReachesBorderY=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->Y = (double)(C->yp);
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}


		xe = x2 - oppcoords.Pc_c1;
		ye = oppcoords.Pb_b1 + y1;
		ze = oppcoords.Pa_a1 + z1;

	}else if(Solids&PositiveZSolid &&Solids&NegativeZSolid && Solids&PositiveYSolid){//Solid on z2,z1 &y2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Positive Z Solid, Negative Z Solid, Positive Y Solid" << endl; 
		}

		//a=z, b=-x, c=-y
		oppcoords.Velocity_b1 = -u2;//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = -u1;//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = -v1;//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = Pz - z1;
		oppcoords.Pb_b1 = x2 - Px;
		oppcoords.Pc_c1 = y2 - Py;


		Case3OppositeSolids(&oppcoords);

		if(oppcoords.ReachedNegativeBBorder){

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;

			ReachesBorderX=true;			

		}

		if(oppcoords.ReachedPositiveBBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;
		}
		
		if(oppcoords.ReachedPositiveCBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;
		}


		xe = x2 - oppcoords.Pb_b1;
		ye = y2 - oppcoords.Pc_c1;
		ze = oppcoords.Pa_a1 + z1;

	}

	if(!ReachesBorderX){		//If reached boundary, new position has already been set
		P->X = xe;
	}
	if(!ReachesBorderY){
		P->Y = ye;
	}
	if(!ReachesBorderZ){
		P->Z = ze;
	}

	//cout<<"Voxel: "<<P->VoxelX<<'\t'<<P->VoxelY<<'\t'<<P->VoxelZ<<endl;

	//cout<<"Final Position: "<<P->X <<'\t'<<P->Y<<'\t'<<P->Z<<endl<<endl;

}

struct Case4OppositeSolidsCoords{
double Velocity_c1;				//Inward velocity to voxel in c direction
double Velocity_c2;				//Outward velocity to voxel in c direction

double Pa_a1;					//Position in a direction within voxel (Pa-a1)
double Pb_b1;					//Position in b direction within voxel (Pb-b1)
double Pc_c1;					//Position in c direction within voxel (Pc-c1)

bool ReachedPositiveCBorder;	//Function returns whether particle reached border
bool ReachedNegativeCBorder;

double* dt;						//Advection time pointer
};


//Advects particles in a voxel with two opposite solids on the a1 and a2 faces, used as a general case
//Dimensions referred to as x,y,z in the function, but correspond to a,b,c in the struct to reduce
//  confusion when transforming into the general case in Case2Solids()
void Case4OppositeSolids(Case4OppositeSolidsCoords* Case){

	double w1 = Case->Velocity_c1;
	double w2 = Case->Velocity_c2;

	double Px = Case->Pa_a1;
	double Py = Case->Pb_b1;
	double Pz = Case->Pc_c1;

	double x1 = 0;	//Imaginary voxel at 0,0,0
	double x2 = 1;
	double y1 = 0;
	double y2 = 1;
	double z1 = 0;
	double z2 = 1;

	double Tz;

	double* dt = Case->dt;

	double u = 0;
	double v = 0;
	double w = (w1 + (w2-w1)*Pz)*36*(1-Px)*Px*(1-Py)*Py;	//double w = 36*w1*(1-Px)*Px*(1-Py)*Py;			// = 36*w1*(x2-Px)*(Px-x1)*(y2-Py)*(Py-y1)/(dx*dx*dy*dy);

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 4 Opposite Solids:" << endl;
		AdvectDebug << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	
	if(w==0){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
		Tz = DBL_MAX;
	}else if(fabs(w2-w1) > Min_dV){	//Extension to the model on account of LB flux non-conservation
		if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
			Tz = log( w2/(w1 + (w2-w1)*Pz) ) / (36*(w2-w1)*(1-Px)*Px*(1-Py)*Py);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 3" << endl; }
			Tz = log( w1/(w1 + (w2-w1)*Pz) ) / (36*(w2-w1)*(1-Px)*Px*(1-Py)*Py);
		}
	}else{
		if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 4" << endl; }
			Tz = (1-Pz)/w;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 5" << endl; }
			Tz = -Pz/w;
		}
	}

	double t;

	if(Tz>(*dt)){									//Remains in same voxel

		t = (*dt);									//Advect for dt
		(*dt) = 0;									//Set dt to 0, entire dt walk carried out
			
	}else{											//Advect to the border
		
		t = Tz;										//Advect for tau_z
		(*dt) -= Tz;								//Reduce dt by tau_z time used

		if(w<0){									//Reaches border with -ve neighbour
			Case->ReachedNegativeCBorder=true;
		}else{
			Case->ReachedPositiveCBorder=true;
		}

	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: ( [DBL_MAX], [DBL_MAX], " << Tz << ")" << endl;
	}

	double ze;

	if(Case->ReachedNegativeCBorder){
		ze = 0;
	}else if(Case->ReachedPositiveCBorder){
		ze = 1;
	}else if(fabs(w2-w1) > Min_dV){	//Extension to the model on account of LB flux non-conservation
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
		double _w1_w2w1 = w1/(w2-w1);
		ze = (Pz + _w1_w2w1)*exp( 36*(w2-w1)*Px*(1-Px)*Py*(1-Py)*t ) - _w1_w2w1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
		ze = Pz + w*t;
	}

	Case->Pc_c1 = ze;

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Raw Position: (" << Px << ", " << Py << ", " << ze << ")" << endl;
	}

}

struct Case4AdjacentSolidsCoords{
double Velocity_b1;				//Inward velocity into b1 face
double Velocity_c2;				//Outward velocity from c2 face

double Pa_a1;					//Position in a direction within voxel (Pa-a1)
double Pb_b1;					//Position in b direction within voxel (Pb-b1)
double Pc_c1;					//Position in c direction within voxel (Pc-c1)
	
bool ReachedPositiveCBorder;	//Function returns whether particle reached border
bool ReachedNegativeBBorder;

double* dt;						//Advection time pointer
};

//Advects particles in a voxel with two adjacent solids on the a1 and b1 faces, used as a general case
//Dimensions referred to as x,y,z in the function, but correspond to a,b,c in the struct to reduce
//  confusion when transforming into the general case in Case2Solids()
void Case4AdjacentSolids(Case4AdjacentSolidsCoords* Case){

	double v1 = Case->Velocity_b1;

	double w2 = Case->Velocity_c2;

	double Px = Case->Pa_a1;
	double Py = Case->Pb_b1;
	double Pz = Case->Pc_c1;

	double x1 = 0;	//Imaginary voxel at 0,0,0
	double x2 = 1;
	double y1 = 0;
	double y2 = 1;
	double z1 = 0;
	double z2 = 1;

	double* dt = Case->dt;

	double A = 12*(x2-Px)*(Px-x1)*(y2-Py)*(Pz-z1);		// = 12*(y2-Py)*(Pz-z1)*(x2-Px)*(Px-x1) / (dx*dx*dy*dz);

	double u = 0;
	double v = A*v1*(y2-Py);							// = 12*v1*(x2-Px)*(Px-x1)*(y2-Py)*(y2-Py)*(Pz-z1) / (dx*dx*dy*dy*dz);
	double w = A*w2*(Pz-z1);							// = 12*w2*(x2-Px)*(Px-x1)*(y2-Py)*(Pz-z1)*(Pz-z1) / (dx*dx*dy*dz*dz);

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Case 4 Adjacent Solids:" << endl;
		AdvectDebug << "v1 = " << v1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
		AdvectDebug << "A = " << A << endl;
	}

	double Ty;
	double Tz;

	//Calculate time to reach the border, Tau

	//Tau y
	
	if(v1<0){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 1" << endl; }
		Ty = log( 1-Py )/(A*v1);						// = ( (dx*dx*dy*dy*dz)/(12*v1*(y2-Py)*(Pz-z1)*(x2-Px)*(Px-x1)) ) * log( (y2-Py)/dy );
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Ty 2" << endl; }
		Ty = DBL_MAX;
	}

	//Tau z
	if(w2>0){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 1" << endl; }
		Tz = -log( Pz )/(A*w2);							// = ( (dx*dx*dy*dz*dz)/(12*w2*(y2-Py)*(Pz-z1)*(x2-Px)*(Px-x1)) ) * log( dz/(Pz-z1) );
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "Tz 2" << endl; }
		Tz = DBL_MAX;
	}

	double t;

	if(Ty<=Tz){		//Reaches y boundary first

		if(Ty>(*dt)){						//Remains in same voxel

			t = (*dt);							//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Ty;								//Advect for tau_y
			(*dt) -= Ty;						//Reduce dt by tau_y time used

			Case->ReachedNegativeBBorder=true;

		}

	}else if(Tz<=Ty){		//Reaches z boundary first

		if(Tz>(*dt)){						//Remains in same voxel

			t = (*dt);							//Advect for dt
			(*dt) = 0;							//Set dt to 0, entire dt walk carried out
			
		}else{									//Advect to the border
			
			t = Tz;								//Advect for tau_z
			(*dt) -= Tz;						//Reduce dt by tau_z time used

			Case->ReachedPositiveCBorder=true;	//Reaches w2 border

		}

	}

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Taus: ( [DBL_MAX], " << Ty << ", " << Tz << ")" << endl;
	}

	//Positions
	double ye, ze;

	if(Case->ReachedNegativeBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 1" << endl; }
		ye = 0;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ye 2" << endl; }
		ye = 1 - (1 - Py) * exp( -A*v1*t );	// = y2 - (y2 - Py) * exp( -((12*v1*(Px-x1)*(x2-Px)*(y2-Py)*(Pz-z1))/(dx*dx*dy*dy*dz))*t );
	}

	if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 1" << endl; }
		ze = 1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle){ AdvectDebug << "ze 2" << endl; }
		ze = Pz * exp( A*w2*t );				// = z1 + (Pz - z1) * exp( ((12*w2*(Px-x1)*(x2-Px)*(y2-Py)*(Pz-z1))/(dx*dx*dy*dz*dz))*t );
	}

	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

	if(DebugAdvectionMode>=2 && DebugParticle){
		AdvectDebug << "Raw Position: (" << Px << ", " << ye << ", " << ze << ")" << endl;
	}

}

//Advection function used when four neighbouring voxels are solid (15 cases)
void CaseFourSolids(int i,Coords* C, Voxel* LatticeElement, double* dt){

	int x = C->x;
	int y = C->y;
	int z = C->z;

	double u1 = LatticeElement->u1;
	double u2 = LatticeElement->u2;
	double v1 = LatticeElement->v1;
	double v2 = LatticeElement->v2;
	double w1 = LatticeElement->w1;
	double w2 = LatticeElement->w2;

	double x1 = (double)x;
	double x2 = (double)(x1+1);
	double y1 = (double)y;
	double y2 = (double)(y1+1);
	double z1 = (double)z;
	double z2 = (double)(z1+1);

	Particle* P = &(Particles[i]);

	double Px = P->X;					//X coordinate of particle
	double Py = P->Y;					//Y coordinate of particle
	double Pz = P->Z;					//Z coordinate of particle

	bool ReachesBorderX=false;
	bool ReachesBorderY=false;
	bool ReachesBorderZ=false;

	double xe, ye, ze;			//New positions after advection

	//cout << "Position: (" << Px << ", " << Py << ", " << Pz << "), Voxel: (" << x << ", " << y << ", " << z << ")" << endl;

	//calculate tau, time to exit voxel in each direction

	Case4OppositeSolidsCoords OppCoords;

	OppCoords.dt=dt;
	OppCoords.ReachedNegativeCBorder=false;
	OppCoords.ReachedPositiveCBorder=false;

	Case4AdjacentSolidsCoords AdjCoords;

	AdjCoords.dt=dt;
	AdjCoords.ReachedNegativeBBorder=false;
	AdjCoords.ReachedPositiveCBorder=false;

	unsigned char Solids = (unsigned char)(LatticeElement->Solid);

	if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveYSolid && Solids&NegativeZSolid ){	//(1) Solid on x1, x2, z1 & y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive X Solid, Negative Z Solid, Positive Y Solid" << endl;
		}

		AdjCoords.Velocity_b1 = v1;
		AdjCoords.Velocity_c2 = w2;

		AdjCoords.Pa_a1 = Px - x1;
		AdjCoords.Pb_b1 = Py - y1;
		AdjCoords.Pc_c1 = Pz - z1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		xe = AdjCoords.Pa_a1 + x1;
		ye = AdjCoords.Pb_b1 + y1;
		ze = AdjCoords.Pc_c1 + z1;

	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveZSolid && Solids&NegativeYSolid ){	//(2) Solid on x1, x2, z2 & y1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive X Solid, Positive Z Solid, Negative Y Solid" << endl;
		}

		AdjCoords.Velocity_b1 = w1;
		AdjCoords.Velocity_c2 = v2;

		AdjCoords.Pa_a1 = x2 - Px;
		AdjCoords.Pb_b1 = Pz - z1;
		AdjCoords.Pc_c1 = Py - y1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}

		xe = x2 - AdjCoords.Pa_a1;
		ye = AdjCoords.Pc_c1 + y1;
		ze = AdjCoords.Pb_b1 + z1;

	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveZSolid && Solids&PositiveYSolid ){	//(3) Solid on x1, x2, z2 & y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive X Solid, Positive Z Solid, Positive Y Solid" << endl;
		}
	
		AdjCoords.Velocity_b1 = w1;
		AdjCoords.Velocity_c2 = -v1;

		AdjCoords.Pa_a1 = Px - x1;
		AdjCoords.Pb_b1 = Pz - z1;
		AdjCoords.Pc_c1 = y2 - Py;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		xe = AdjCoords.Pa_a1 + x1;
		ye = y2 - AdjCoords.Pc_c1;
		ze = AdjCoords.Pb_b1 + z1;

	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&NegativeZSolid && Solids&NegativeYSolid ){	//(4) Solid on x1, x2, z1 & y1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive X Solid, Negative Z Solid, Negative Y Solid" << endl;
		}
	
		AdjCoords.Velocity_b1 = -w2;
		AdjCoords.Velocity_c2 = v2;

		AdjCoords.Pa_a1 = Px - x1;
		AdjCoords.Pb_b1 = z2 - Pz;
		AdjCoords.Pc_c1 = Py - y1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){

			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}

		xe = AdjCoords.Pa_a1 + x1;
		ye = AdjCoords.Pc_c1 + y1;
		ze = z2 - AdjCoords.Pb_b1;

	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&PositiveXSolid && Solids&NegativeZSolid ){	//(5) Solid on y1, y2, z1 & x2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid, Positive Y Solid, Negative Z Solid, Positive X Solid" << endl;
		}

		AdjCoords.Velocity_b1 = u1;
		AdjCoords.Velocity_c2 = w2;

		AdjCoords.Pa_a1 = y2 - Py;
		AdjCoords.Pb_b1 = Px - x1;
		AdjCoords.Pc_c1 = Pz - z1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		xe = AdjCoords.Pb_b1 + x1;
		ye = y2 - AdjCoords.Pa_a1;
		ze = AdjCoords.Pc_c1 + z1;

	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&PositiveZSolid && Solids&NegativeXSolid ){	//(6) Solid on y1, y2, z2 & x1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid, Positive Y Solid, Positive Z Solid, Negative X Solid" << endl;
		}
		
		AdjCoords.Velocity_b1 = w1;
		AdjCoords.Velocity_c2 = u2;

		AdjCoords.Pa_a1 = Py - y1;
		AdjCoords.Pb_b1 = Pz - z1;
		AdjCoords.Pc_c1 = Px - x1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		xe = AdjCoords.Pc_c1 + x1;
		ye = AdjCoords.Pa_a1 + y1;
		ze = AdjCoords.Pb_b1 + z1;

	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&PositiveZSolid && Solids&PositiveXSolid ){	//(7) Solid on y1, y2, z2 & x2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid, Positive Y Solid, Positive Z Solid, Positive X Solid" << endl;
		}
	
		AdjCoords.Velocity_b1 = u1;
		AdjCoords.Velocity_c2 = -w1;

		AdjCoords.Pa_a1 = Py - y1;
		AdjCoords.Pb_b1 = Px - x1;
		AdjCoords.Pc_c1 = z2 - Pz;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		xe = AdjCoords.Pb_b1 + x1;
		ye = AdjCoords.Pa_a1 + y1;
		ze = z2 - AdjCoords.Pc_c1;

	}else if(Solids&PositiveYSolid && Solids&NegativeYSolid && Solids&NegativeZSolid && Solids&NegativeXSolid ){	//(8) Solid on y1, y2, z1 & x1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Y Solid, Positive Y Solid, Negative Z Solid, Negative X Solid" << endl;
		}
	
		AdjCoords.Velocity_b1 = -u2;
		AdjCoords.Velocity_c2 = w2;

		AdjCoords.Pa_a1 = Py - y1;
		AdjCoords.Pb_b1 = x2 - Px;
		AdjCoords.Pc_c1 = Pz - z1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		xe = x2 - AdjCoords.Pb_b1;
		ye = AdjCoords.Pa_a1 + y1;
		ze = AdjCoords.Pc_c1 + z1;

	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&PositiveYSolid && Solids&NegativeXSolid ){	//(9) Solid on z1, z2, x1 & y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid, Positive Z Solid, Negative X Solid, Positive Y Solid" << endl;
		}

		AdjCoords.Velocity_b1 = v1;
		AdjCoords.Velocity_c2 = u2;

		AdjCoords.Pa_a1 = z2 - Pz;
		AdjCoords.Pb_b1 = Py - y1;
		AdjCoords.Pc_c1 = Px - x1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		xe = AdjCoords.Pc_c1 + x1;
		ye = AdjCoords.Pb_b1 + y1;
		ze = z2 - AdjCoords.Pa_a1;
		
	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&PositiveXSolid && Solids&NegativeYSolid ){	//(10) Solid on z1, z2, x2 & y1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid, Positive Z Solid, Positive X Solid, Negative Y Solid" << endl;
		}

		AdjCoords.Velocity_b1 = u1;
		AdjCoords.Velocity_c2 = v2;

		AdjCoords.Pa_a1 = Pz - z1;
		AdjCoords.Pb_b1 = Px - x1;
		AdjCoords.Pc_c1 = Py - y1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}

		xe = AdjCoords.Pb_b1 + x1;
		ye = AdjCoords.Pc_c1 + y1;
		ze = AdjCoords.Pa_a1 + z1;
		
	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&PositiveXSolid && Solids&PositiveYSolid ){	//(11) Solid on z1, z2, x2 & y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid, Positive Z Solid, Positive X Solid, Positive Y Solid" << endl;
		}
	
		AdjCoords.Velocity_b1 = v1;
		AdjCoords.Velocity_c2 = -u1;

		AdjCoords.Pa_a1 = Pz - z1;
		AdjCoords.Pb_b1 = Py - y1;
		AdjCoords.Pc_c1 = x2 - Px;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){
		
			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		xe = x2 - AdjCoords.Pc_c1;
		ye = AdjCoords.Pb_b1 + y1;
		ze = AdjCoords.Pa_a1 + z1;

	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&NegativeXSolid && Solids&NegativeYSolid ){	//(12) Solid on z1, z2, x1 & y1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid, Positive Z Solid, Negative X Solid, Negative Y Solid" << endl;
		}
	
		AdjCoords.Velocity_b1 = -v2;
		AdjCoords.Velocity_c2 = u2;

		AdjCoords.Pa_a1 = Pz - z1;
		AdjCoords.Pb_b1 = y2 - Py;
		AdjCoords.Pc_c1 = Px - x1;

		Case4AdjacentSolids(&AdjCoords);

		if(AdjCoords.ReachedNegativeBBorder){

			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}

		if(AdjCoords.ReachedPositiveCBorder){
			
			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		xe = AdjCoords.Pc_c1 + x1;
		ye = y2 - AdjCoords.Pb_b1;
		ze = AdjCoords.Pa_a1 + z1;

	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveYSolid && Solids&NegativeYSolid ){	//Solid on x1, x2, y2 & y1 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive X Solid, Positive Y Solid, Negative Y Solid" << endl;
		}

		OppCoords.Velocity_c1 = w1;
		OppCoords.Velocity_c2 = w2;

		OppCoords.Pa_a1 = Px - x1;
		OppCoords.Pb_b1 = Py - y1;
		OppCoords.Pc_c1 = Pz - z1;

		Case4OppositeSolids(&OppCoords);

		if(OppCoords.ReachedNegativeCBorder){

			P->Z = (double)(C->zn+1);
			P->VoxelZ = C->zn;
			
			ReachesBorderZ=true;

		}

		if(OppCoords.ReachedPositiveCBorder){
			
			P->Z = (double)C->zp;
			P->VoxelZ = C->zp;
			
			ReachesBorderZ=true;

		}

		xe = OppCoords.Pa_a1 + x1;
		ye = OppCoords.Pb_b1 + y1;
		ze = OppCoords.Pc_c1 + z1;

	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveZSolid && Solids&NegativeZSolid ){	//Solid on x1, x2, z1 & z2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative X Solid, Positive X Solid, Negative Z Solid, Positive Z Solid" << endl;
		}

		OppCoords.Velocity_c1 = v1;
		OppCoords.Velocity_c2 = v2;

		OppCoords.Pa_a1 = Px - x1;
		OppCoords.Pb_b1 = z2 - Pz;
		OppCoords.Pc_c1 = Py - y1;

		Case4OppositeSolids(&OppCoords);

		if(OppCoords.ReachedNegativeCBorder){

			P->Y = (double)(C->yn+1);
			P->VoxelY = C->yn;
			
			ReachesBorderY=true;

		}

		if(OppCoords.ReachedPositiveCBorder){
			
			P->Y = (double)C->yp;
			P->VoxelY = C->yp;
			
			ReachesBorderY=true;

		}

		xe = OppCoords.Pa_a1 + x1;
		ye = OppCoords.Pc_c1 + y1;
		ze = z2 - OppCoords.Pb_b1;

	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&PositiveYSolid && Solids&NegativeYSolid ){	//Solid on z1, z2, y1 & y2 sides

		if(DebugAdvectionMode>=2 && DebugParticle){
			AdvectDebug << "Negative Z Solid, Positive Z Solid, Negative Y Solid, Positive Y Solid" << endl;
		}
		
		OppCoords.Velocity_c1 = u1;
		OppCoords.Velocity_c2 = u2;

		OppCoords.Pa_a1 = z2 - Pz;
		OppCoords.Pb_b1 = Py - y1;
		OppCoords.Pc_c1 = Px - x1;

		Case4OppositeSolids(&OppCoords);

		if(OppCoords.ReachedNegativeCBorder){

			P->X = (double)(C->xn+1);
			P->VoxelX = C->xn;
			
			ReachesBorderX=true;

		}

		if(OppCoords.ReachedPositiveCBorder){
			
			P->X = (double)C->xp;
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;

		}

		xe = OppCoords.Pc_c1 + x1;
		ye = OppCoords.Pb_b1 + y1;
		ze = z2 - OppCoords.Pa_a1;

	}	

	
	if(!ReachesBorderX){		//If reached boundary, new position has already been set
		P->X = (double)xe;
	}
	if(!ReachesBorderY){
		P->Y = (double)ye;
	}
	if(!ReachesBorderZ){
		P->Z = (double)ze;
	}

	//cout << "New position: (" << Particles[i].X << ", " << Particles[i].Y << ", " << Particles[i].Z << ")" << endl << endl;

}

struct DistributionDatum{
double XPos;
double XCount;
double YPos;
double YCount;
double ZPos;
double ZCount;
};

//Non-MPI outputs particle distribution on the lattice (not propagator)
void OutputParticleDistribution(int XBinsPerVoxel, int YBinsPerVoxel, int ZBinsPerVoxel){

	ofstream outfile(OutputFilePath(DistributionOutputFile));

	if(outfile==0){
		return;
	}
	
	outfile << "X\tCount\tY\tCount\tZ\tCount" << endl;

	int XCount = XBinsPerVoxel*NLattice_x;
	int YCount = YBinsPerVoxel*NLattice_y;
	int ZCount = ZBinsPerVoxel*NLattice_z;

	int ArrCount = max(XCount,max(YCount,XCount));
	DistributionDatum* Data = new DistributionDatum[ArrCount];

	int i=0;
	while(i!=ArrCount){
		
		Data[i].XPos=-1;
		Data[i].XCount=0;
		Data[i].YPos=-1;
		Data[i].YCount=0;
		Data[i].ZPos=-1;
		Data[i].ZCount=0;

		i++;
	}

	double XPos=0;
	double dX = 1.0/((double)XBinsPerVoxel);
	int x=0;

	while(x!=XCount){

		Data[x].XPos=XPos;

		i=0;
		while(i!=NParticles){

			if(Particles[i].X>=XPos && Particles[i].X<(XPos+dX)){
				Data[x].XCount++;
			}

			i++;
		}

		x++;
		XPos+=dX;
	}

	double YPos=0;
	double dY = 1.0/((double)YBinsPerVoxel);
	int y=0;

	while(y!=YCount){

		Data[y].YPos=YPos;

		i=0;
		while(i!=NParticles){

			if(Particles[i].Y>=YPos && Particles[i].Y<(YPos+dY)){
				Data[y].YCount++;
			}

			i++;
		}

		y++;
		YPos+=dY;
	}

	double ZPos=0;
	double dZ = 1.0/((double)ZBinsPerVoxel);
	int z=0;

	while(z!=ZCount){

		Data[z].ZPos=ZPos;

		i=0;
		while(i!=NParticles){

			if(Particles[i].Z>=ZPos && Particles[i].Z<(ZPos+dZ)){
				Data[z].ZCount++;
			}

			i++;
		}

		z++;
		ZPos+=dZ;
	}

	i=0;
	while(i!=ArrCount){

		if(Data[i].XPos>=0){
			outfile << Data[i].XPos << '\t' << Data[i].XCount << '\t';
		}else{
			outfile << "\t\t";
		}

		if(Data[i].YPos>=0){
			outfile << Data[i].YPos << '\t' << Data[i].YCount << '\t';
		}else{
			outfile << "\t\t";
		}

		if(Data[i].ZPos>=0){
			outfile << Data[i].ZPos << '\t' << Data[i].ZCount;
		}else{
			outfile << "\t";
		}

		outfile << endl;

		i++;
	}

	delete[] Data;

}

void _DecompositionReadSolids(){	//Reads solid nodes from file

	ifstream InFile;
	InFile.open(InputFilePath(SolidsFileName));
	
	if(InFile!=0){
		cout << "Reading in solids from '" << SolidsFileName << "'" << endl;
	}else{
		cout << "Could not open solids file" << endl;
		return;
	}

	InFile.unsetf(ios_base::skipws);

	InFile.seekg(0,ios::end);
	__int64 flength = InFile.tellg();	//File length
	InFile.seekg(0,ios::beg);

	const __int64 buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	
	__int64 nRead = 0;				//Number of chars read from file in total
	__int64 readl = 0;				//Number of chars read from file into buffer
	__int64 readpos = 0;			//Read position in buffer

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				cout << "Error: reached end of geometry file before lattice entirely read" << endl;
				delete[] buf;
				InFile.close();
				return;
			}
			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]=='1'||buf[readpos]=='0'){

			_DecompositionSetLattice(x,y,z,(buf[readpos]=='1'));

			x++;
			if(x==NLattice_x){
				x=0;
				y++;
				if(y==NLattice_y){
					y=0;
					z++;
					if(z==NLattice_z){
						delete[] buf;
						InFile.close();
						return;
					}
				}
			}

		}

		readpos++;
	}

}

#if Using_MPI && MPI_Mode==1

//Decomposes lattice and sends data to all threads
int LatticeDecomposition(int NProcessors, int ThreadID){

	ofstream DecompositionLog;
	DecompositionLog.open(OutputFilePath(DecompositionLogFile));

	int NArray = ((NLattice_x*NLattice_y*NLattice_z)>>3) + 1;
	_DecompositionLattice = new char[NArray];

	memset(_DecompositionLattice,0,NArray);
	_DecompositionReadSolids();

	double ycoefficient = (double)NLattice_y / (double)NLattice_x;
	double zcoefficient = (double)NLattice_z / (double)NLattice_x;

	double N = pow((double)NProcessors / (ycoefficient*zcoefficient) , 1.0/3.0);

	int nZ = (int)ceil( zcoefficient*N );				//Number of CPUs to divide up Z direction

	if(!MPIMode1_DecomposeZ){	//No Z decomposition
		nZ = 1;
	}

	int* NonSolidsZSlices = new int[NLattice_z];		//Count number of non solids in each slice in z
	int* NonSolidsZCumulative = new int[NLattice_z];	//Cumulative number of non solids in z
	int TotalNonSolids = 0;

	int z=0;
	while(z!=NLattice_z){
		NonSolidsZSlices[z]=0;

		int x=0;
		while(x!=NLattice_x){
			int y=0;
			while(y!=NLattice_y){
				if(!_DecompositionGetLattice(x,y,z)){
					NonSolidsZSlices[z]++;
					TotalNonSolids++;
				}
				y++;
			}
			x++;
		}

		if(z!=0){
			NonSolidsZCumulative[z] = NonSolidsZCumulative[z-1] + NonSolidsZSlices[z];
		}else{
			NonSolidsZCumulative[0] = NonSolidsZSlices[0];
		}

		z++;
	}

	double Porosity = ((double)TotalNonSolids)/((double)(NLattice_x*NLattice_y*NLattice_z));
	double NonSolidsPerCPU = ((double)TotalNonSolids)/((double)NProcessors);

	DecompositionLog << "Porosity: " << Porosity << endl;
	DecompositionLog << "Non-solids per CPU: " << NonSolidsPerCPU << endl << endl;

	LatticeZBlock* Plane = new LatticeZBlock[nZ];			//Divided into blocks in Z according to nZ

	int Remainder = NProcessors%nZ;
	int NPerBlock = (NProcessors - Remainder)/nZ;

	int i=0;
	while(i!=nZ){
		Plane[i].NSquares = NPerBlock;
		i++;
	}

	i=0;
	while(Remainder>0){
		Plane[i].NSquares++;
		Remainder--;

		i++;
		if(i==nZ){
			i=0;
		}
	}

	double NSolidsPerSlice = ((double)TotalNonSolids)/((double)nZ);

	DecompositionLog << "Dividing z into " << nZ << " partitions of average " << NSolidsPerSlice << " non-solids" << endl;

	int z0 = 0;
	double z0_ = 0;

	int z1 = 0;
	double z1_ = 0;

	double NSolidsPreceding = 0;

	i=0;
	while(i!=nZ){

		double NSolidsInSlice = ((double)(Plane[i].NSquares) * NonSolidsPerCPU);

		while(true){

			double N0 = 0;
			if(z1!=0){
				N0 = NonSolidsZCumulative[z1-1] - NSolidsPreceding;
			}
			double N1 = NonSolidsZCumulative[z1] - NSolidsPreceding;

			if(z1==NLattice_z){

				z1_ = 0;
				break;

			}else if(N0<NSolidsInSlice && N1>=NSolidsInSlice){

				z1_ = (NSolidsInSlice - N0)/NonSolidsZSlices[z1];
				break;

			}

			z1++;
		}

		NSolidsPreceding = NonSolidsZCumulative[z1-1] + z1_*NonSolidsZSlices[z1];

		Plane[i].z0 = z0;
		Plane[i].z0_ = z0_;
		Plane[i].z1 = z1;
		Plane[i].z1_ = z1_;
		Plane[i].nZ_d = ((double)z1 + z1_) - ((double)z0 + z0_);

		DecompositionLog << "Block range: " << ((double)z0 + z0_) << " -> " << ((double)z1 + z1_) << " with " << NSolidsInSlice << " non-solids (" << Plane[i].NSquares << " CPUs)" << endl;

		z0 = z1;
		z0_ = z1_;

		i++;
	}

	if(MPIMode1_DecomposeZ){

	//Obtain minimum error combination of rounded borders

	int nBorders = nZ - 1;

	unsigned long long FlagInt = (0xFFFFFFFFFFFFFFFF << nBorders);	//1111111111....1110000 (nZ rightmost 0s)

	unsigned long long MinErrorCombination;
	double MinError = DBL_MAX;

	while(FlagInt!=0){	//Wait until it overflows

		unsigned long long Bits = FlagInt;

		i=0;
		double Error = 0;
		while(i!=nBorders){

			double NSolidsInSlice = ((double)(Plane[i].NSquares) * NonSolidsPerCPU);

			bool Round = (bool)(Bits&1);	//1 or 0

			Plane[i].UpperRound = Round;
			Plane[i+1].LowerRound = Round;

			Bits >>= 1;

			if(i==0){

				int rIndex;
				
				if(Round){
					rIndex = (int)ceil(Plane[0].z1_);
				}else{
					rIndex = (int)floor(Plane[0].z1_);
				}

				int zIndex = Plane[0].z1 + rIndex;

				double delta = NonSolidsZCumulative[zIndex-1] - NSolidsInSlice;

				Error += (delta*delta);


				if(nBorders==1){	//2 Partitions

					NSolidsInSlice = ((double)(Plane[1].NSquares) * NonSolidsPerCPU);

					if(Plane[0].UpperRound){
						rIndex = (int)ceil(Plane[1].z0_);
					}else{
						rIndex = (int)floor(Plane[1].z0_);
					}

					zIndex = Plane[1].z0 + rIndex;

					delta = (NonSolidsZCumulative[NLattice_z-1] - NonSolidsZCumulative[zIndex-1]) - NSolidsInSlice;

					Error += (delta*delta);

				}

			}else{

				int rIndex0;
				if(Plane[i].LowerRound){
					rIndex0 = (int)ceil(Plane[i].z0_);
				}else{
					rIndex0 = (int)floor(Plane[i].z0_);
				}

				int rIndex1;
				if(Plane[i].UpperRound){
					rIndex1 = (int)ceil(Plane[i].z1_);
				}else{
					rIndex1 = (int)floor(Plane[i].z1_);
				}

				int zIndex0 = Plane[i].z0 + rIndex0;
				int zIndex1 = Plane[i].z1 + rIndex1;

				double delta = (NonSolidsZCumulative[zIndex1-1] - NonSolidsZCumulative[zIndex0-1]) - NSolidsInSlice;

				Error += (delta*delta);


				if(i==(nBorders-1)){

					NSolidsInSlice = ((double)(Plane[nBorders].NSquares) * NonSolidsPerCPU);

					if(Plane[i].LowerRound){
						rIndex0 = (int)ceil(Plane[nBorders].z0_);
					}else{
						rIndex0 = (int)floor(Plane[nBorders].z0_);
					}

					zIndex0 = Plane[nBorders].z0 + rIndex0;

					delta = (NonSolidsZCumulative[NLattice_z-1] - NonSolidsZCumulative[zIndex0-1]) - NSolidsInSlice;

					Error += (delta*delta);

				}

			}

			i++;
		}

		if(Error < MinError){

			MinErrorCombination = FlagInt;
			MinError = Error;

		}

		FlagInt++;
	}

	//Set lowest error combination

	DecompositionLog << endl;
	DecompositionLog << "Rounding borders to the minimum error combination: " << endl;

	unsigned long long Bits = MinErrorCombination;

	i=0;
	while(i!=nBorders){

		double NSolidsInSlice = ((double)(Plane[i].NSquares) * NonSolidsPerCPU);

		bool Round = (bool)(Bits&1);	//1 or 0

		Plane[i].UpperRound = Round;
		Plane[i+1].LowerRound = Round;

		Bits >>= 1;

		if(i==0){

			int rIndex;
				
			if(Round){
				rIndex = (int)ceil(Plane[0].z1_);
			}else{
				rIndex = (int)floor(Plane[0].z1_);
			}

			int zIndex = Plane[0].z1 + rIndex;

			double delta = NonSolidsZCumulative[zIndex-1] - NSolidsInSlice;
			double Err = fabs(delta)/NSolidsInSlice;

			Plane[0].z0Rounded = 0;
			Plane[0].z1Rounded = zIndex;

			DecompositionLog << "z from 0 to " << zIndex << ": Error = " << Err << endl;


			if(nBorders==1){		//2 Partitions

				NSolidsInSlice = ((double)(Plane[1].NSquares) * NonSolidsPerCPU);

				if(Plane[0].UpperRound){
					rIndex = (int)ceil(Plane[1].z0_);
				}else{
					rIndex = (int)floor(Plane[1].z0_);
				}

				zIndex = Plane[1].z0 + rIndex;

				delta = (NonSolidsZCumulative[NLattice_z-1] - NonSolidsZCumulative[zIndex-1]) - NSolidsInSlice;
				Err = fabs(delta)/NSolidsInSlice;

				Plane[nBorders].z0Rounded = zIndex;
				Plane[nBorders].z1Rounded = NLattice_z;

				DecompositionLog << "z from " << zIndex << " to " << NLattice_z << ": Error = " << Err << endl;

			}

		}else{

			int rIndex0;
			if(Plane[i].LowerRound){
				rIndex0 = (int)ceil(Plane[i].z0_);
			}else{
				rIndex0 = (int)floor(Plane[i].z0_);
			}

			int rIndex1;
			if(Plane[i].UpperRound){
				rIndex1 = (int)ceil(Plane[i].z1_);
			}else{
				rIndex1 = (int)floor(Plane[i].z1_);
			}

			int zIndex0 = Plane[i].z0 + rIndex0;
			int zIndex1 = Plane[i].z1 + rIndex1;

			double delta = (NonSolidsZCumulative[zIndex1-1] - NonSolidsZCumulative[zIndex0-1]) - NSolidsInSlice;
			double Err = fabs(delta)/NSolidsInSlice;

			Plane[i].z0Rounded = zIndex0;
			Plane[i].z1Rounded = zIndex1;

			DecompositionLog << "z from " << zIndex0 << " to " << zIndex1 << ": Error = " << Err << endl;

			if(i==(nBorders-1)){

				NSolidsInSlice = ((double)(Plane[nBorders].NSquares) * NonSolidsPerCPU);

				if(Plane[i].LowerRound){
					rIndex0 = (int)ceil(Plane[nBorders].z0_);
				}else{
					rIndex0 = (int)floor(Plane[nBorders].z0_);
				}

				zIndex0 = Plane[nBorders].z0 + rIndex0;

				delta = (NonSolidsZCumulative[NLattice_z-1] - NonSolidsZCumulative[zIndex0-1]) - NSolidsInSlice;
				Err = fabs(delta)/NSolidsInSlice;

				Plane[nBorders].z0Rounded = zIndex0;
				Plane[nBorders].z1Rounded = NLattice_z;

				DecompositionLog << "z from " << zIndex0 << " to " << NLattice_z << ": Error = " << Err << endl;

			}

		}

		i++;
	}

	}else{		//Not decomposing Z

		Plane[0].z0Rounded = 0;
		Plane[0].z1Rounded = NLattice_z;

	}

	//For each z block, divide into x,y cuboids

	int zBlock = 0;
	while(zBlock!=nZ){

		double ycoefficient = (double)NLattice_y / (double)NLattice_x;

		double N = sqrt( (double)Plane[zBlock].NSquares/ycoefficient );

		int nY = (int)ceil( ycoefficient*N );	//Number of blocks in y

		int* NonSolidsYPlane = new int[NLattice_y];			//Count number of non solids in each slice in y of this z block
		int* NonSolidsYCumulative = new int[NLattice_y];	//Cumulative number of non solids in z
		int TotalNonSolidsY = 0;

		Plane[zBlock].YBlocks = new LatticeYBlock[nY];
		Plane[zBlock].YBlockNum = nY;

		Remainder = (Plane[zBlock].NSquares)%nY;
		NPerBlock = (Plane[zBlock].NSquares - Remainder)/nY;

		i=0;
		while(i!=nY){
			Plane[zBlock].YBlocks[i].NSquares = NPerBlock;	//Number of processor blocks in y
			i++;
		}

		i=0;
		while(Remainder>0){
			Plane[zBlock].YBlocks[i].NSquares++;
			Remainder--;

			i++;
			if(i==nY){
				i=0;
			}
		}

		int y = 0;
		while(y!=NLattice_y){

			NonSolidsYPlane[y] = 0;

			int z = Plane[zBlock].z0Rounded;
			while(z!=Plane[zBlock].z1Rounded){

				int x = 0;
				while(x!=NLattice_x){

					if(!_DecompositionGetLattice(x,y,z)){
						NonSolidsYPlane[y]++;
						TotalNonSolidsY++;
					}

				x++;
				}

			z++;
			}

			NonSolidsYCumulative[y] = TotalNonSolidsY;

		y++;
		}

		DecompositionLog << endl;
		DecompositionLog << "z partition " << zBlock << " containing " << TotalNonSolidsY << " non-solids" << endl;

		double NonSolidsPerCPU = ((double)TotalNonSolidsY)/((double)Plane[zBlock].NSquares);

		DecompositionLog << "Non-solids per CPU: " << NonSolidsPerCPU << endl;

		int Remainder = Plane[zBlock].NSquares%nY;
		int NPerBlock = (Plane[zBlock].NSquares - Remainder)/nY;

		int i=0;
		while(i!=nY){
			Plane[zBlock].YBlocks[i].NSquares = NPerBlock;
			i++;
		}

		i=0;
		while(Remainder>0){
			Plane[zBlock].YBlocks[i].NSquares++;
			Remainder--;

			i++;
			if(i==nY){
				i=0;
			}
		}

		double NSolidsPerSlice = ((double)TotalNonSolidsY)/((double)nY);

		DecompositionLog << '\t' << "Dividing y into " << nY << " partitions of average " << NSolidsPerSlice << " non-solids" << endl;

		int y0 = 0;
		double y0_ = 0;

		int y1 = 0;
		double y1_ = 0;

		double NSolidsPreceding = 0;

		i=0;
		while(i!=nY){

			double NSolidsInSlice = (((double)(Plane[zBlock].YBlocks[i].NSquares)) * NonSolidsPerCPU);

			while(true){

				double N0 = 0;
				if(y1!=0){
					N0 = NonSolidsYCumulative[y1-1] - NSolidsPreceding;
				}
				double N1 = NonSolidsYCumulative[y1] - NSolidsPreceding;

				if(y1==NLattice_y){

					y1_ = 0;
					break;

				}else if(N0<NSolidsInSlice && N1>=NSolidsInSlice){

					y1_ = (NSolidsInSlice - N0)/NonSolidsYPlane[y1];
					break;

				}

				y1++;
			}

			NSolidsPreceding = NonSolidsYCumulative[y1-1] + y1_*NonSolidsYPlane[y1];

			Plane[zBlock].YBlocks[i].y0 = y0;
			Plane[zBlock].YBlocks[i].y0_ = y0_;
			Plane[zBlock].YBlocks[i].y1 = y1;
			Plane[zBlock].YBlocks[i].y1_ = y1_;
			Plane[zBlock].YBlocks[i].nY_d = ((double)y1 + y1_) - ((double)y0 + y0_);

			DecompositionLog << '\t' << "Block range: " << ((double)y0 + y0_) << " -> " << ((double)y1 + y1_) << " with " << NSolidsInSlice << " non-solids (" << Plane[zBlock].YBlocks[i].NSquares << " CPUs)" << endl;
						y0 = y1;
			y0_ = y1_;

			i++;
		}

		//Obtain minimum error combination of rounded borders

		int nBorders = nY - 1;

		unsigned long long FlagInt = (0xFFFFFFFFFFFFFFFF << nBorders);	//1111111111....1110000 (nZ rightmost 0s)

		unsigned long long MinErrorCombination;
		double MinError = DBL_MAX;

		while(FlagInt!=0){	//Wait until it overflows

			unsigned long long Bits = FlagInt;

			i=0;
			double Error = 0;
			while(i!=nBorders){

				double NSolidsInSlice = (((double)(Plane[zBlock].YBlocks[i].NSquares)) * NonSolidsPerCPU);

				bool Round = (bool)(Bits&1);	//1 or 0

				Plane[zBlock].YBlocks[i].UpperRound = Round;
				Plane[zBlock].YBlocks[i+1].LowerRound = Round;

				Bits >>= 1;

				if(i==0){

					int rIndex;
				
					if(Round){
						rIndex = (int)ceil(Plane[zBlock].YBlocks[0].y1_);
					}else{
						rIndex = (int)floor(Plane[zBlock].YBlocks[0].y1_);
					}

					int yIndex = Plane[zBlock].YBlocks[0].y1 + rIndex;

					double delta = NonSolidsYCumulative[yIndex-1] - NSolidsInSlice;

					Error += (delta*delta);


					if(nBorders==1){

						NSolidsInSlice = ((double)(Plane[zBlock].YBlocks[1].NSquares) * NonSolidsPerCPU);

						if(Plane[zBlock].YBlocks[i].UpperRound){
							rIndex = (int)ceil(Plane[zBlock].YBlocks[1].y0_);
						}else{
							rIndex = (int)floor(Plane[zBlock].YBlocks[1].y0_);
						}

						yIndex = Plane[zBlock].YBlocks[1].y0 + rIndex;

						delta = (NonSolidsYCumulative[NLattice_y-1] - NonSolidsYCumulative[yIndex-1]) - NSolidsInSlice;

						Error += (delta*delta);

					}

				}else{

					int rIndex0;
					if(Plane[zBlock].YBlocks[i].LowerRound){
						rIndex0 = (int)ceil(Plane[zBlock].YBlocks[i].y0_);
					}else{
						rIndex0 = (int)floor(Plane[zBlock].YBlocks[i].y0_);
					}

					int rIndex1;
					if(Plane[zBlock].YBlocks[i].UpperRound){
						rIndex1 = (int)ceil(Plane[zBlock].YBlocks[i].y1_);
					}else{
						rIndex1 = (int)floor(Plane[zBlock].YBlocks[i].y1_);
					}

					int yIndex0 = Plane[zBlock].YBlocks[i].y0 + rIndex0;
					int yIndex1 = Plane[zBlock].YBlocks[i].y1 + rIndex1;

					double delta = (NonSolidsYCumulative[yIndex1-1] - NonSolidsYCumulative[yIndex0-1]) - NSolidsInSlice;

					Error += (delta*delta);


					if(i==(nBorders-1)){

						NSolidsInSlice = ((double)(Plane[zBlock].YBlocks[nBorders].NSquares) * NonSolidsPerCPU);

						if(Plane[zBlock].YBlocks[i].LowerRound){
							rIndex0 = (int)ceil(Plane[zBlock].YBlocks[nBorders].y0_);
						}else{
							rIndex0 = (int)floor(Plane[zBlock].YBlocks[nBorders].y0_);
						}

						yIndex0 = Plane[zBlock].YBlocks[nBorders].y0 + rIndex0;

						delta = (NonSolidsYCumulative[NLattice_y-1] - NonSolidsYCumulative[yIndex0-1]) - NSolidsInSlice;

						Error += (delta*delta);

					}

				}

				i++;
			}

			if(Error < MinError){

				MinErrorCombination = FlagInt;
				MinError = Error;

			}

			FlagInt++;
		}

		//Set lowest error combination

		DecompositionLog << '\t' << endl;
		DecompositionLog << '\t' << "Rounding borders to the minimum error combination: " << endl;

		unsigned long long Bits = MinErrorCombination;

		i=0;
		while(i!=nBorders){

			double NSolidsInSlice = ((double)(Plane[zBlock].YBlocks[i].NSquares) * NonSolidsPerCPU);

			bool Round = (bool)(Bits&1);	//1 or 0

			Plane[zBlock].YBlocks[i].UpperRound = Round;
			Plane[zBlock].YBlocks[i+1].LowerRound = Round;

			Bits >>= 1;

			if(i==0){

				int rIndex;
				
				if(Round){
					rIndex = (int)ceil(Plane[zBlock].YBlocks[0].y1_);
				}else{
					rIndex = (int)floor(Plane[zBlock].YBlocks[0].y1_);
				}

				int yIndex = Plane[zBlock].YBlocks[0].y1 + rIndex;

				double delta = NonSolidsYCumulative[yIndex-1] - NSolidsInSlice;
				double Err = fabs(delta)/NSolidsInSlice;

				Plane[zBlock].YBlocks[0].y0Rounded = 0;
				Plane[zBlock].YBlocks[0].y1Rounded = yIndex;

				DecompositionLog << '\t' << "y from 0 to " << yIndex << ": Error = " << Err << endl;


				if(nBorders==1){	//2 Partitions

					NSolidsInSlice = ((double)(Plane[zBlock].YBlocks[1].NSquares) * NonSolidsPerCPU);

					if(Plane[zBlock].YBlocks[0].UpperRound){
						rIndex = (int)ceil(Plane[zBlock].YBlocks[1].y0_);
					}else{
						rIndex = (int)floor(Plane[zBlock].YBlocks[1].y0_);
					}

					yIndex = Plane[zBlock].YBlocks[1].y0 + rIndex;

					delta = (NonSolidsYCumulative[NLattice_y-1] - NonSolidsYCumulative[yIndex-1]) - NSolidsInSlice;
					Err = fabs(delta)/NSolidsInSlice;

					Plane[zBlock].YBlocks[1].y0Rounded = yIndex;
					Plane[zBlock].YBlocks[1].y1Rounded = NLattice_y;

					DecompositionLog << '\t' << "y from " << yIndex << " to " << NLattice_y << ": Error = " << Err << endl;

				}

			}else{

				int rIndex0;
				if(Plane[zBlock].YBlocks[i].LowerRound){
					rIndex0 = (int)ceil(Plane[zBlock].YBlocks[i].y0_);
				}else{
					rIndex0 = (int)floor(Plane[zBlock].YBlocks[i].y0_);
				}

				int rIndex1;
				if(Plane[zBlock].YBlocks[i].UpperRound){
					rIndex1 = (int)ceil(Plane[zBlock].YBlocks[i].y1_);
				}else{
					rIndex1 = (int)floor(Plane[zBlock].YBlocks[i].y1_);
				}

				int yIndex0 = Plane[zBlock].YBlocks[i].y0 + rIndex0;
				int yIndex1 = Plane[zBlock].YBlocks[i].y1 + rIndex1;

				double delta = (NonSolidsYCumulative[yIndex1-1] - NonSolidsYCumulative[yIndex0-1]) - NSolidsInSlice;
				double Err = fabs(delta)/NSolidsInSlice;

				Plane[zBlock].YBlocks[i].y0Rounded = yIndex0;
				Plane[zBlock].YBlocks[i].y1Rounded = yIndex1;

				DecompositionLog << '\t' << "y from " << yIndex0 << " to " << yIndex1 << ": Error = " << Err << endl;

				if(i==(nBorders-1)){

					NSolidsInSlice = ((double)(Plane[zBlock].YBlocks[nBorders].NSquares) * NonSolidsPerCPU);

					if(Plane[zBlock].YBlocks[i].LowerRound){
						rIndex0 = (int)ceil(Plane[zBlock].YBlocks[nBorders].y0_);
					}else{
						rIndex0 = (int)floor(Plane[zBlock].YBlocks[nBorders].y0_);
					}

					yIndex0 = Plane[zBlock].YBlocks[nBorders].y0 + rIndex0;

					delta = (NonSolidsYCumulative[NLattice_y-1] - NonSolidsYCumulative[yIndex0-1]) - NSolidsInSlice;
					Err = fabs(delta)/NSolidsInSlice;

					Plane[zBlock].YBlocks[nBorders].y0Rounded = yIndex0;
					Plane[zBlock].YBlocks[nBorders].y1Rounded = NLattice_y;

					DecompositionLog << '\t' << "y from " << yIndex0 << " to " << NLattice_y << ": Error = " << Err << endl;

				}

			}

			i++;
		}

		//For each y block, divide in x to obtain decomposed domain

		int yBlock = 0;
		while(yBlock!=nY){

			LatticeYBlock* Block = &(Plane[zBlock].YBlocks[yBlock]);

			int nX = Block->NSquares;	//Number of blocks in x

			int* NonSolidsXPlane = new int[NLattice_x];			//Count number of non solids in each slice in x of this y block
			int* NonSolidsXCumulative = new int[NLattice_x];	//Cumulative number of non solids in x
			int TotalNonSolidsX = 0;

			Block->XBlocks = new LatticeXBlock[nX];
			Block->XBlockNum = nX;

			int x = 0;
			while(x!=NLattice_x){

				NonSolidsXPlane[x] = 0;

				int z = Plane[zBlock].z0Rounded;
				while(z!=Plane[zBlock].z1Rounded){

					int y = Block->y0Rounded;
					while(y!=Block->y1Rounded){

						if(!_DecompositionGetLattice(x,y,z)){
							NonSolidsXPlane[x]++;
							TotalNonSolidsX++;
						}

					y++;
					}

				z++;
				}

				NonSolidsXCumulative[x] = TotalNonSolidsX;

			x++;
			}

			DecompositionLog << endl;
			DecompositionLog << '\t' << "y partition " << yBlock << " containing " << TotalNonSolidsX << " non-solids" << endl;

			double NonSolidsPerCPU = ((double)TotalNonSolidsX)/((double)Block->NSquares);

			DecompositionLog << "\tNon-solids per CPU: " << NonSolidsPerCPU << endl;

			double NSolidsPerSlice = ((double)TotalNonSolidsX)/((double)nX);

			DecompositionLog << '\t' << "Dividing x into " << nX << " partitions of average " << NSolidsPerSlice << " non-solids" << endl;

			int x0 = 0;
			double x0_ = 0;

			int x1 = 0;
			double x1_ = 0;

			double NSolidsPreceding = 0;

			i=0;
			while(i!=nX){

				while(true){

					double N0 = 0;
					if(x1!=0){
						N0 = NonSolidsXCumulative[x1-1] - NSolidsPreceding;
					}
					double N1 = NonSolidsXCumulative[x1] - NSolidsPreceding;

					if(x1==NLattice_x){

						x1_ = 0;
						break;

					}else if(N0<NonSolidsPerCPU && N1>=NonSolidsPerCPU){

						x1_ = (NonSolidsPerCPU - N0)/NonSolidsXPlane[x1];
						break;

					}

					x1++;
				}

				NSolidsPreceding = NonSolidsXCumulative[x1-1] + x1_*NonSolidsXPlane[x1];

				Block->XBlocks[i].x0 = x0;
				Block->XBlocks[i].x0_ = x0_;
				Block->XBlocks[i].x1 = x1;
				Block->XBlocks[i].x1_ = x1_;
				Block->XBlocks[i].nX_d = ((double)x1 + x1_) - ((double)x0 + x0_);

				DecompositionLog << "\t\t" << "Block range: " << ((double)x0 + x0_) << " -> " << ((double)x1 + x1_) << " with " << NonSolidsPerCPU << " non-solids" << endl;

				x0 = x1;
				x0_ = x1_;

				i++;
			}

			//Obtain minimum error combination of rounded borders

			int nBorders = nX - 1;

			unsigned long long FlagInt = (0xFFFFFFFFFFFFFFFF << nBorders);	//1111111111....1110000 (nZ rightmost 0s)

			unsigned long long MinErrorCombination;
			double MinError = DBL_MAX;

			while(FlagInt!=0){	//Wait until it overflows

				unsigned long long Bits = FlagInt;

				i=0;
				double Error = 0;
				while(i!=nBorders){

					bool Round = (bool)(Bits&1);	//1 or 0

					Block->XBlocks[i].UpperRound = Round;
					Block->XBlocks[i+1].LowerRound = Round;

					Bits >>= 1;

					if(i==0){

						int rIndex;
				
						if(Round){
							rIndex = (int)ceil(Block->XBlocks[0].x1_);
						}else{
							rIndex = (int)floor(Block->XBlocks[0].x1_);
						}

						int xIndex = Block->XBlocks[0].x1 + rIndex;

						double delta = NonSolidsXCumulative[xIndex-1] - NonSolidsPerCPU;

						Error += (delta*delta);


						if(nBorders==1){

							if(Block->XBlocks[0].UpperRound){
								rIndex = (int)ceil(Block->XBlocks[1].x0_);
							}else{
								rIndex = (int)floor(Block->XBlocks[1].x0_);
							}

							xIndex = Block->XBlocks[1].x0 + rIndex;

							delta = (NonSolidsXCumulative[NLattice_x-1] - NonSolidsXCumulative[xIndex-1]) - NonSolidsPerCPU;

							Error += (delta*delta);

						}

					}else{

						int rIndex0;
						if(Block->XBlocks[i].LowerRound){
							rIndex0 = (int)ceil(Block->XBlocks[i].x0_);
						}else{
							rIndex0 = (int)floor(Block->XBlocks[i].x0_);
						}

						int rIndex1;
						if(Block->XBlocks[i].UpperRound){
							rIndex1 = (int)ceil(Block->XBlocks[i].x1_);
						}else{
							rIndex1 = (int)floor(Block->XBlocks[i].x1_);
						}

						int xIndex0 = Block->XBlocks[i].x0 + rIndex0;
						int xIndex1 = Block->XBlocks[i].x1 + rIndex1;

						double delta = (NonSolidsXCumulative[xIndex1-1] - NonSolidsXCumulative[xIndex0-1]) - NonSolidsPerCPU;

						Error += (delta*delta);


						if(i==(nBorders-1)){

							if(Block->XBlocks[i].LowerRound){
								rIndex0 = (int)ceil(Block->XBlocks[nBorders].x0_);
							}else{
								rIndex0 = (int)floor(Block->XBlocks[nBorders].x0_);
							}

							xIndex0 = Block->XBlocks[nBorders].x0 + rIndex0;

							delta = (NonSolidsXCumulative[NLattice_x-1] - NonSolidsXCumulative[xIndex0-1]) - NonSolidsPerCPU;

							Error += (delta*delta);

						}

					}

					i++;
				}

				if(Error < MinError){

					MinErrorCombination = FlagInt;
					MinError = Error;

				}

				FlagInt++;
			}

			//Set lowest error combination

			DecompositionLog << "\t\t" << endl;
			DecompositionLog << "\t\t" << "Rounding borders to the minimum error combination: " << endl;

			unsigned long long Bits = MinErrorCombination;

			i=0;
			while(i!=nBorders){

				bool Round = (bool)(Bits&1);	//1 or 0

				Block->XBlocks[i].UpperRound = Round;
				Block->XBlocks[i+1].LowerRound = Round;

				Bits >>= 1;

				if(i==0){

					int rIndex;
				
					if(Round){
						rIndex = (int)ceil(Block->XBlocks[0].x1_);
					}else{
						rIndex = (int)floor(Block->XBlocks[0].x1_);
					}

					int xIndex = Block->XBlocks[0].x1 + rIndex;

					double delta = NonSolidsXCumulative[xIndex-1] - NonSolidsPerCPU;
					double Err = fabs(delta)/NonSolidsPerCPU;

					Block->XBlocks[0].x0Rounded = 0;
					Block->XBlocks[0].x1Rounded = xIndex;

					DecompositionLog << "\t\t" << "x from 0 to " << xIndex << ": Error = " << Err << endl;


					if(nBorders==1){

						if(Block->XBlocks[0].UpperRound){
							rIndex = (int)ceil(Block->XBlocks[1].x0_);
						}else{
							rIndex = (int)floor(Block->XBlocks[1].x0_);
						}

						xIndex = Block->XBlocks[1].x0 + rIndex;

						delta = (NonSolidsXCumulative[NLattice_x-1] - NonSolidsXCumulative[xIndex-1]) - NonSolidsPerCPU;
						Err = fabs(delta)/NonSolidsPerCPU;

						Block->XBlocks[1].x0Rounded = xIndex;
						Block->XBlocks[1].x1Rounded = NLattice_x;

						DecompositionLog << "\t\t" << "x from " << xIndex << " to " << NLattice_x << ": Error = " << Err << endl;

					}

				}else{

					int rIndex0;
					if(Block->XBlocks[i].LowerRound){
						rIndex0 = (int)ceil(Block->XBlocks[i].x0_);
					}else{
						rIndex0 = (int)floor(Block->XBlocks[i].x0_);
					}

					int rIndex1;
					if(Block->XBlocks[i].UpperRound){
						rIndex1 = (int)ceil(Block->XBlocks[i].x1_);
					}else{
						rIndex1 = (int)floor(Block->XBlocks[i].x1_);
					}

					int xIndex0 = Block->XBlocks[i].x0 + rIndex0;
					int xIndex1 = Block->XBlocks[i].x1 + rIndex1;

					double delta = (NonSolidsXCumulative[xIndex1-1] - NonSolidsXCumulative[xIndex0-1]) - NonSolidsPerCPU;
					double Err = fabs(delta)/NonSolidsPerCPU;

					Block->XBlocks[i].x0Rounded = xIndex0;
					Block->XBlocks[i].x1Rounded = xIndex1;

					DecompositionLog << "\t\t" << "x from " << xIndex0 << " to " << xIndex1 << ": Error = " << Err << endl;


					if(i==(nBorders-1)){

						if(Block->XBlocks[i].LowerRound){
							rIndex0 = (int)ceil(Block->XBlocks[nBorders].x0_);
						}else{
							rIndex0 = (int)floor(Block->XBlocks[nBorders].x0_);
						}

						xIndex0 = Block->XBlocks[nBorders].x0 + rIndex0;

						delta = (NonSolidsXCumulative[NLattice_x-1] - NonSolidsXCumulative[xIndex0-1]) - NonSolidsPerCPU;
						Err = fabs(delta)/NonSolidsPerCPU;

						Block->XBlocks[nBorders].x0Rounded = xIndex0;
						Block->XBlocks[nBorders].x1Rounded = NLattice_x;

						DecompositionLog << "\t\t" << "x from " << xIndex0 << " to " << NLattice_x << ": Error = " << Err << endl;

					}

				}

				i++;
			}

			delete[] NonSolidsXPlane;
			delete[] NonSolidsXCumulative;

			yBlock++;
		}

		delete[] NonSolidsYPlane;
		delete[] NonSolidsYCumulative;

		zBlock++;
	}

	delete[] NonSolidsZSlices;
	delete[] NonSolidsZCumulative;

	//Organise CPU domains

	Domains = new CPUDomain[NProcessors];		//Thread 0 domains

	i=0;
	z = 0;
	while(z!=nZ){

		LatticeZBlock* ZBlock = &(Plane[z]);

		int y=0;
		while(y!=ZBlock->YBlockNum){

			LatticeYBlock* YBlock = &(ZBlock->YBlocks[y]);

			int x=0;
			while(x!=YBlock->XBlockNum){

				LatticeXBlock* XBlock = &(YBlock->XBlocks[x]);

				Domains[i].x0 = XBlock->x0Rounded;
				Domains[i].x1 = XBlock->x1Rounded-1;
				Domains[i].y0 = YBlock->y0Rounded;
				Domains[i].y1 = YBlock->y1Rounded-1;
				Domains[i].z0 = ZBlock->z0Rounded;
				Domains[i].z1 = ZBlock->z1Rounded-1;

				/*
				Coords C;
				GetCoordinates(Domains[i].x0,Domains[i].y0,Domains[i].z0,&C);	//Gets coordinates of adjoining voxels

				Domains[i].xn = C.xn;
				Domains[i].yn = C.yn;
				Domains[i].zn = C.zn;

				GetCoordinates(Domains[i].x1,Domains[i].y1,Domains[i].z1,&C);

				Domains[i].xp = C.xp;
				Domains[i].yp = C.yp;
				Domains[i].zp = C.zp;
				
				Domains[i].DomainWidthX = Domains[i].x1 - Domains[i].x0 + 3;
				Domains[i].DomainWidthY = Domains[i].y1 - Domains[i].y0 + 3;
				Domains[i].DomainWidthZ = Domains[i].z1 - Domains[i].z0 + 3;
				*/

				Domains[i].DomainWidthX = Domains[i].x1 - Domains[i].x0 + 1;
				Domains[i].DomainWidthY = Domains[i].y1 - Domains[i].y0 + 1;
				Domains[i].DomainWidthZ = Domains[i].z1 - Domains[i].z0 + 1;

				DecompositionLog << Domains[i].x0 << " -> " << Domains[i].x1 << ", " << Domains[i].y0 << " -> " << Domains[i].y1 << ", " << Domains[i].z0 << " -> " << Domains[i].z1 << endl;
				
				i++;

				x++;
			}

			y++;
		}

		z++;
	}

	i=0;
	while(i!=NProcessors){

		Domains[i].NeighbourDomainsXn = new int[NProcessors];
		Domains[i].nNeighbourDomainsXn = 0;
		Domains[i].NeighbourDomainsXp = new int[NProcessors];
		Domains[i].nNeighbourDomainsXp = 0;
		Domains[i].NeighbourDomainsYn = new int[NProcessors];
		Domains[i].nNeighbourDomainsYn = 0;
		Domains[i].NeighbourDomainsYp = new int[NProcessors];
		Domains[i].nNeighbourDomainsYp = 0;
		Domains[i].NeighbourDomainsZn = new int[NProcessors];
		Domains[i].nNeighbourDomainsZn = 0;
		Domains[i].NeighbourDomainsZp = new int[NProcessors];
		Domains[i].nNeighbourDomainsZp = 0;

		DecompositionLog << "[" << i << "]: " << Domains[i].x0 << " -> " << Domains[i].x1 << ", " << Domains[i].y0 << " -> " << Domains[i].y1 << ", " << Domains[i].z0 << " -> " << Domains[i].z1 << endl;

		Coords Cn, Cp;
		GetCoordinates(Domains[i].x0,Domains[i].y0,Domains[i].z0,&Cn);	//Gets coordinates of adjoining voxels in negative direction
		GetCoordinates(Domains[i].x1,Domains[i].y1,Domains[i].z1,&Cp);	//Gets coordinates of adjoining voxels in positive direction

		int i2=0;
		while(i2!=NProcessors){

			//if(Domains[i2].x1 == Domains[i].xn){
			if(Domains[i2].x1 == Cn.xn){
				if(((Domains[i2].y0>=Domains[i].y0&&Domains[i2].y0<=Domains[i].y1)||(Domains[i2].y1>=Domains[i].y0&&Domains[i2].y1<=Domains[i].y1)||(Domains[i2].y0<Domains[i].y0&&Domains[i2].y1>Domains[i].y1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsXn[Domains[i].nNeighbourDomainsXn] = i2;
					Domains[i].nNeighbourDomainsXn++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in xn" << endl;

				}
			}
			//if(Domains[i2].x0 == Domains[i].xp){
			if(Domains[i2].x0 == Cp.xp){
				if(((Domains[i2].y0>=Domains[i].y0&&Domains[i2].y0<=Domains[i].y1)||(Domains[i2].y1>=Domains[i].y0&&Domains[i2].y1<=Domains[i].y1)||(Domains[i2].y0<Domains[i].y0&&Domains[i2].y1>Domains[i].y1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsXp[Domains[i].nNeighbourDomainsXp] = i2;
					Domains[i].nNeighbourDomainsXp++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in xp" << endl;
				}
			}
			//if(Domains[i2].y1 == Domains[i].yn){
			if(Domains[i2].y1 == Cn.yn){
				if(((Domains[i2].x0>=Domains[i].x0&&Domains[i2].x0<=Domains[i].x1)||(Domains[i2].x1>=Domains[i].x0&&Domains[i2].x1<=Domains[i].x1)||(Domains[i2].x0<Domains[i].x0&&Domains[i2].x1>Domains[i].x1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsYn[Domains[i].nNeighbourDomainsYn] = i2;
					Domains[i].nNeighbourDomainsYn++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in yn" << endl;
				}
			}
			//if(Domains[i2].y0 == Domains[i].yp){
			if(Domains[i2].y0 == Cp.yp){
				if(((Domains[i2].x0>=Domains[i].x0&&Domains[i2].x0<=Domains[i].x1)||(Domains[i2].x1>=Domains[i].x0&&Domains[i2].x1<=Domains[i].x1)||(Domains[i2].x0<Domains[i].x0&&Domains[i2].x1>Domains[i].x1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsYp[Domains[i].nNeighbourDomainsYp] = i2;
					Domains[i].nNeighbourDomainsYp++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in yp" << endl;
				}
			}
			//if(Domains[i2].z1 == Domains[i].zn){
			if(Domains[i2].z1 == Cn.zn){
				if(((Domains[i2].x0>=Domains[i].x0&&Domains[i2].x0<=Domains[i].x1)||(Domains[i2].x1>=Domains[i].x0&&Domains[i2].x1<=Domains[i].x1)||(Domains[i2].x0<Domains[i].x0&&Domains[i2].x1>Domains[i].x1))
					&& ((Domains[i2].y0>=Domains[i].y0&&Domains[i2].y0<=Domains[i].y1)||(Domains[i2].y1>=Domains[i].y0&&Domains[i2].y1<=Domains[i].y1)||(Domains[i2].y0<Domains[i].y0&&Domains[i2].y1>Domains[i].y1))){
					
					Domains[i].NeighbourDomainsZn[Domains[i].nNeighbourDomainsZn] = i2;
					Domains[i].nNeighbourDomainsZn++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in zn" << endl;
				}
			}
			//if(Domains[i2].z0 == Domains[i].zp){
			if(Domains[i2].z0 == Cp.zp){
				if(((Domains[i2].x0>=Domains[i].x0&&Domains[i2].x0<=Domains[i].x1)||(Domains[i2].x1>=Domains[i].x0&&Domains[i2].x1<=Domains[i].x1)||(Domains[i2].x0<Domains[i].x0&&Domains[i2].x1>Domains[i].x1))
					&& ((Domains[i2].y0>=Domains[i].y0&&Domains[i2].y0<=Domains[i].y1)||(Domains[i2].y1>=Domains[i].y0&&Domains[i2].y1<=Domains[i].y1)||(Domains[i2].y0<Domains[i].y0&&Domains[i2].y1>Domains[i].y1))){
				
					Domains[i].NeighbourDomainsZp[Domains[i].nNeighbourDomainsZp] = i2;
					Domains[i].nNeighbourDomainsZp++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in zp" << endl;
				}
			}

			i2++;
		}
		DecompositionLog << endl;

		i++;
	}

	//Open vector field file and create array of pointers to z,y coordinates in file

	ifstream fin(InputFilePath(VelocitiesFileName),ios::binary);

	if(fin!=0){
		cout << "Reading in velocity field from '" << VelocitiesFileName << "'" << endl;
		DecompositionLog << "Reading in velocity field from '" << InputFilePath(VelocitiesFileName) << "'" << endl;
	}else{
		cout << "Could not open velocity field file" << endl;
		DecompositionLog << "Could not open velocity field file" << endl;
		return 2;
	}

	__int64 FileLength;
	__int64* FilePointers = new __int64[NLattice_z*NLattice_y];		//Location in file stream of each z,y coordinate values
	_DecompositionGetStreamPtrArray(fin, &FileLength, FilePointers);

	//Send structs and geometry to other threads

	i=1;
	while(i!=NProcessors){

		MPI_Send(Domains, NProcessors*sizeof(CPUDomain), MPI_BYTE, i, 0, MPI_COMM_WORLD);		//Send CPUDomain structs

		int ci = 0;
		while(ci!=NProcessors){
			MPI_Send(Domains[ci].NeighbourDomainsXn, sizeof(int)*Domains[ci].nNeighbourDomainsXn, MPI_BYTE, i, 0, MPI_COMM_WORLD);
			MPI_Send(Domains[ci].NeighbourDomainsXp, sizeof(int)*Domains[ci].nNeighbourDomainsXp, MPI_BYTE, i, 0, MPI_COMM_WORLD);
			MPI_Send(Domains[ci].NeighbourDomainsYn, sizeof(int)*Domains[ci].nNeighbourDomainsYn, MPI_BYTE, i, 0, MPI_COMM_WORLD);
			MPI_Send(Domains[ci].NeighbourDomainsYp, sizeof(int)*Domains[ci].nNeighbourDomainsYp, MPI_BYTE, i, 0, MPI_COMM_WORLD);
			MPI_Send(Domains[ci].NeighbourDomainsZn, sizeof(int)*Domains[ci].nNeighbourDomainsZn, MPI_BYTE, i, 0, MPI_COMM_WORLD);
			MPI_Send(Domains[ci].NeighbourDomainsZp, sizeof(int)*Domains[ci].nNeighbourDomainsZp, MPI_BYTE, i, 0, MPI_COMM_WORLD);
		ci++;
		}

		//Send geometry

		int nX = Domains[i].x1 - Domains[i].x0 + 1;
		int nY = Domains[i].y1 - Domains[i].y0 + 1;
		int nZ = Domains[i].z1 - Domains[i].z0 + 1;

		//int nVoxels = (nX+2)*(nY+2)*(nZ+2);			//Domain overlap
		int nVoxels = nX*nY*nZ;

		unsigned char* Geometry = new unsigned char[nVoxels];

		if(Geometry==0){	//Memory allocation failed
			cout << "Memory allocation failed at (0): stopping execution" << endl;
			DecompositionLog << "Memory allocation failed at (0)" << endl;
			delete[] FilePointers;
			fin.close();
			return 1;
		}

		int nNonSolids = 0;

		/*
		int z = -1;
		while(z!=(nZ+1)){
		int y = -1;
		while(y!=(nY+1)){
		int x = -1;
		while(x!=(nX+1)){
		*/
		int z = 0;
		while(z!=nZ){
		int y = 0;
		while(y!=nY){
		int x = 0;
		while(x!=nX){
					
			int VoxelX = Domains[i].x0 + x;
			int VoxelY = Domains[i].y0 + y;
			int VoxelZ = Domains[i].z0 + z;

			//CycleVoxelCoordinates(&VoxelX,&VoxelY,&VoxelZ);

			//int Index = (nX+2)*((nY+2)*(z+1) + (y+1)) + (x+1);
			int Index = nX*(nY*z + y) + x;

			if(_DecompositionGetLattice(VoxelX,VoxelY,VoxelZ)){		//Solid
				Geometry[Index] = VoxelSolid;
			}else{
				Geometry[Index] =  _DecompositionGetVoxelSolid(VoxelX,VoxelY,VoxelZ);
				nNonSolids++;
			}

		x++;
		}
		y++;
		}
		z++;
		}

		MPI_Send(&nNonSolids, 1, MPI_INT, i, 0, MPI_COMM_WORLD);	//Send number of non solids

		char ThreadStatusFlag;	//Check thread initialised memory properly
		MPI_Status Stat;
		MPI_Recv(&ThreadStatusFlag, 1, MPI_CHAR, i, 0, MPI_COMM_WORLD, &Stat);

		if(ThreadStatusFlag==1){
			DecompositionLog << "Memory allocation failed at thread [" << i << "]: stopping execution" << endl;
			cout << "Memory allocation failed at thread [" << i << "]: stopping execution" << endl;
			delete[] Geometry;
			delete[] FilePointers;
			fin.close();
			return 1;
		}

		int dCount = 0;
		int dLength = 1048576;		//1MB
		while(dCount < nVoxels){
			int l = dLength;
			if(nVoxels-dCount < l){
				l = nVoxels-dCount;
			}

			MPI_Send(&Geometry[dCount], l, MPI_BYTE, i, 0, MPI_COMM_WORLD);			//Send geometry to thread

			dCount += l;
		}

		DecompositionLog << "Sent geometry to thread " << i << endl;
		cout << "Sent geometry to thread " << i << endl;

		delete[] Geometry;

		//Vector field

		double* VectorField = new double[nNonSolids*3];

		if(VectorField==0){	//Memory allocation failed
			cout << "Memory allocation failed at (1): stopping execution" << endl;
			DecompositionLog << "Memory allocation failed at (1)" << endl;
			delete[] FilePointers;
			fin.close();
			return 1;
		}

		_DecompositionGetVectorField(fin, FileLength, FilePointers, VectorField, &(Domains[i]));	

		dCount = 0;
		dLength = 1048576;		//3MB * 8 (double)
		while(dCount < nNonSolids){
			int l = dLength;
			if(nNonSolids-dCount < l){
				l = nNonSolids-dCount;
			}

			MPI_Send(&VectorField[dCount*3], l*3, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);			//Send vector field to thread

			dCount += l;
		}

		DecompositionLog << "Sent vector field to thread " << i << endl;
		cout << "Sent vector field to thread " << i << endl;

		delete[] VectorField;

		i++;
	}

	//Thread 0 domain

	if(ThreadID==0){

		memcpy(&ProcessorDomain,&(Domains[0]),sizeof(CPUDomain));

		int nX = Domains[0].x1 - Domains[0].x0 + 1;
		int nY = Domains[0].y1 - Domains[0].y0 + 1;
		int nZ = Domains[0].z1 - Domains[0].z0 + 1;

		//int nVoxels = (nX+2)*(nY+2)*(nZ+2);			//Domain overlap
		int nVoxels = nX*nY*nZ;

		LatticeBase = new unsigned int[nVoxels];

		//_CPUDomainWidthX = nX+2;
		//_CPUDomainWidthY = nY+2;
		_CPUDomainWidthX = nX;
		_CPUDomainWidthY = nY;
		
		//Geometry

		int nNonSolids = 0;

		/*
		int z = -1;
		while(z!=(nZ+1)){
		int y = -1;
		while(y!=(nY+1)){
		int x = -1;
		while(x!=(nX+1)){
		*/
		int z = 0;
		while(z!=nZ){
		int y = 0;
		while(y!=nY){
		int x = 0;
		while(x!=nX){
					
			int VoxelX = Domains[0].x0 + x;
			int VoxelY = Domains[0].y0 + y;
			int VoxelZ = Domains[0].z0 + z;

			CycleVoxelCoordinates(&VoxelX,&VoxelY,&VoxelZ);

			if(!_DecompositionGetLattice(VoxelX,VoxelY,VoxelZ)){		//Non-Solid
				nNonSolids++;
			}

		x++;
		}
		y++;
		}
		z++;
		}

		Lattice = new Voxel[nNonSolids+1];

		if(LatticeBase==0 || Lattice==0){			//Memory allocation failed
			DecompositionLog << "Memory allocation failed at thread [0]: stopping execution" << endl;
			cout << "Memory allocation failed at thread [0]: stopping execution" << endl;
			delete[] FilePointers;
			fin.close();
			return 1;
		}

		Lattice[0].Solid = VoxelSolid;
		Lattice[0].Vx = 0;
		Lattice[0].Vy = 0;
		Lattice[0].Vz = 0;

		int LatticeIndex = 1;

		/*
		z = -1;
		while(z!=(nZ+1)){
		int y = -1;
		while(y!=(nY+1)){
		int x = -1;
		while(x!=(nX+1)){
		*/
		z = 0;
		while(z!=nZ){
		int y = 0;
		while(y!=nY){
		int x = 0;
		while(x!=nX){
					
			int VoxelX = Domains[0].x0 + x;
			int VoxelY = Domains[0].y0 + y;
			int VoxelZ = Domains[0].z0 + z;

			//CycleVoxelCoordinates(&VoxelX,&VoxelY,&VoxelZ);
			//int Index = (nX+2)*((nY+2)*(z+1) + (y+1)) + (x+1);

			int Index = nX*(nY*z + y) + x;

			if(_DecompositionGetLattice(VoxelX,VoxelY,VoxelZ)){		//Solid

				LatticeBase[Index] = 0;

			}else{

				LatticeBase[Index] = LatticeIndex;
				Lattice[LatticeIndex].Solid = _DecompositionGetVoxelSolid(VoxelX,VoxelY,VoxelZ);
				LatticeIndex++;

			}

		x++;
		}
		y++;
		}
		z++;
		}

		DecompositionLog << "Sent geometry to thread 0" << endl;
		cout << "Sent geometry to thread 0" << endl;

		//Vectorfield

		_DecompositionGetVectorField(fin, FileLength, FilePointers, &Lattice[1], &(Domains[0]));

		DecompositionLog << "Sent vector field to thread 0" << endl;
		cout << "Sent vector field to thread 0" << endl;

	}

	//

	fin.close();

	//Free memory

	delete[] FilePointers;

	i=0;
	while(i!=nZ){

		int yc=0;
		while(yc!=Plane[i].YBlockNum){
			delete[] Plane[i].YBlocks[yc].XBlocks;
			yc++;
		}

		delete[] Plane[i].YBlocks;
		i++;
	}
	delete[] Plane;

	delete[] _DecompositionLattice;

	cout << "Lattice decomposition completed" << endl;

	DecompositionInitialiseLattice(NProcessors, ThreadID);	//Initialise lattice velocities

	return 0;
}

void _DecompositionGetStreamPtrArray(ifstream& InFile, __int64* FileLength, __int64* StreamPtr){		//Create array of filestream pointers to z,y coordinates in file

	StreamPtr[0] = 0;
	int streamc = 1;

	InFile.seekg(0,ios::end);
	__int64 flength = InFile.tellg();	//File length
	InFile.seekg(0,ios::beg);

	*FileLength = flength;

	const __int64 buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	int component=0;

	bool state = 1;					//Read char (true) or whitespace/return (false)
	
	__int64 nRead = 0;				//Number of chars read from file in total
	__int64 readl = 0;				//Number of chars read from file into buffer
	__int64 readpos = 0;			//Read position in buffer

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				delete[] buf;
				return;
			}
			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]==' '||buf[readpos]=='\n'||buf[readpos]=='\r'){	//Between values
			state = 0;
		}else if(state==0){							//Beginning of next value

			state = 1;

			component++;
			if(component==3){
				
			component=0;
			x++;
			if(x==NLattice_x){

				StreamPtr[streamc] = nRead - readl + readpos;
				streamc++;

				x=0;
				y++;
				if(y==NLattice_y){
					y=0;
					z++;
					if(z==NLattice_z){
						delete[] buf;
						return;
					}
				}
			}

			}

		}

		readpos++;
	}

}

void _DecompositionReadVectorFieldValues(ifstream& InFile, __int64 FileLength, int n, bool StartOnValue){

	if(n==0 && StartOnValue){
		return;
	}

	const __int64 buffersize = 10240;		//10KB

	char* buf = new char[buffersize];

	__int64 nRead = InFile.tellg();		//Number of chars read from file in total
	__int64 readl = 0;					//Number of chars read from file into buffer
	__int64 readpos = 0;				//Read position in buffer

	int state = 0;
	int count = -1;

	if(StartOnValue){
		state = 1;
		count = 0;
	}

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (FileLength-nRead)){
				readl = FileLength-nRead;
			}
			if(readl == 0){
				break;
			}
			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]==' '||buf[readpos]=='\n'||buf[readpos]=='\r'){	//Between values
			state = 0;
		}else if(state==0){							//Beginning of next value
			state = 1;
			count++;

			if(count == n){
				InFile.seekg(nRead-readl+readpos);
				break;
			}

		}

		readpos++;
	}

	delete[] buf;

}

void _DecompositionReadVectorFieldValues(ifstream& InFile, __int64 FileLength, double* Array, int ArrOffset, int n, bool StartOnValue){

	const int buffersize = 10240;		//10KB

	char* buf = new char[buffersize];

	__int64 nRead = InFile.tellg();		//Number of chars read from file in total
	__int64 readl = 0;					//Number of chars read from file into buffer
	__int64 readpos = 0;				//Read position in buffer

	int state = 0;
	int count = -1;

	if(StartOnValue){
		state = 1;
		count = 0;
	}

	char VBuffer[24];				//Value buffer
	int VBufl = 0;

	while(true){

		if(readpos == readl){

			readl = buffersize;

			if(buffersize > (FileLength-nRead)){
				readl = FileLength-nRead;
			}
			if(readl == 0){
				break;
			}

			InFile.read(buf,readl);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]==' '||buf[readpos]=='\n'||buf[readpos]=='\r'){	//Between values

			if(VBufl!=0){

				VBuffer[VBufl]='\0';
				double v = atof(VBuffer);
						
				if(fabs(v) < MIN_VECTOR_RESOLUTION){
					v = 0;
				}

				Array[ArrOffset + count] = v;

				VBufl = 0;
			}

			state = 0;

		}else if(state==0){							//Beginning of next value

			state = 1;
			count++;

			if(count == n){
				InFile.seekg(nRead-readl+readpos);
				break;
			}

			VBuffer[0] = buf[readpos];
			VBufl = 1;

		}else{		//Number 

			VBuffer[VBufl] = buf[readpos];
			VBufl++;

		}

		readpos++;
	}

	delete[] buf;

}

void _DecompositionGetVectorField(ifstream& InFile, __int64 FileLength, __int64* StreamPtr, double* VectorField, CPUDomain* Domain){

	/*
	int nX = Domain->x1 - Domain->x0 + 3;
	int nY = Domain->y1 - Domain->y0 + 3;
	int nZ = Domain->z1 - Domain->z0 + 3;
	*/

	int nX = Domain->x1 - Domain->x0 + 1;
	int nY = Domain->y1 - Domain->y0 + 1;
	int nZ = Domain->z1 - Domain->z0 + 1;

	double* xValues = new double[nX*3];
	int LatticeIndex = 0;

	int zi = 0;
	while(zi!=nZ){
		/*
		int z = Domain->zn + zi;
		while(z>=NLattice_z){ z-=NLattice_z; }
		*/
		int z = Domain->z0 + zi;

		int yi = 0;
		while(yi!=nY){
			/*
			int y = Domain->yn + yi;
			while(y>=NLattice_y){ y-=NLattice_y; }
			*/
			int y = Domain->y0 + yi;

			InFile.seekg(StreamPtr[z*NLattice_y + y]);

			/*
			int x = Domain->xn;
			if(x!=0){
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,0,3,true);			//Read in xn values

			x = Domain->x0;
			if(x != Domain->xn+1){
				InFile.seekg(StreamPtr[z*NLattice_y + y]);
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,3,(nX-2)*3,true);		//Read in x values

			x = Domain->xp;
			if(x != Domain->x1+1){
				InFile.seekg(StreamPtr[z*NLattice_y + y]);
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,(nX-1)*3,3,true);		//Read in xp values

			int i=0;
			while(i!=nX){
				x = Domain->xn + i;
				while(x>=NLattice_x){ x-=NLattice_x; }

				if(!_DecompositionGetLattice(x,y,z)){	//If non solid
					VectorField[LatticeIndex*3] = xValues[i*3];
					VectorField[LatticeIndex*3+1] = xValues[i*3+1];
					VectorField[LatticeIndex*3+2] = xValues[i*3+2];
					LatticeIndex++;
				}
				i++;
			}
			*/

			int x = Domain->x0;
			if(x != 0){
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,0,nX*3,true);		//Read in x values

			int i=0;
			while(i!=nX){
				x = Domain->x0 + i;

				if(!_DecompositionGetLattice(x,y,z)){	//If non solid
					VectorField[LatticeIndex*3] = xValues[i*3];
					VectorField[LatticeIndex*3+1] = xValues[i*3+1];
					VectorField[LatticeIndex*3+2] = xValues[i*3+2];
					LatticeIndex++;
				}
				i++;
			}

			yi++;
		}

		zi++;
	}

	delete[] xValues;

}

void _DecompositionGetVectorField(ifstream& InFile, __int64 FileLength, __int64* StreamPtr, Voxel* VectorField, CPUDomain* Domain){

	/*
	int nX = Domain->x1 - Domain->x0 + 3;
	int nY = Domain->y1 - Domain->y0 + 3;
	int nZ = Domain->z1 - Domain->z0 + 3;
	*/

	int nX = Domain->x1 - Domain->x0 + 1;
	int nY = Domain->y1 - Domain->y0 + 1;
	int nZ = Domain->z1 - Domain->z0 + 1;

	double* xValues = new double[nX*3];
	int LatticeIndex = 0;

	int zi = 0;
	while(zi!=nZ){
		/*
		int z = Domain->zn + zi;
		while(z>=NLattice_z){ z-=NLattice_z; }
		*/
		int z = Domain->z0 + zi;

		int yi = 0;
		while(yi!=nY){
			/*
			int y = Domain->yn + yi;
			while(y>=NLattice_y){ y-=NLattice_y; }
			*/
			int y = Domain->y0 + yi;

			InFile.seekg(StreamPtr[z*NLattice_y + y]);

			/*
			int x = Domain->xn;
			if(x!=0){
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,0,3,true);			//Read in xn values

			x = Domain->x0;
			if(x != Domain->xn+1){
				InFile.seekg(StreamPtr[z*NLattice_y + y]);
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,3,(nX-2)*3,true);		//Read in x values

			x = Domain->xp;
			if(x != Domain->x1+1){
				InFile.seekg(StreamPtr[z*NLattice_y + y]);
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,(nX-1)*3,3,true);		//Read in xp values

			int i=0;
			while(i!=nX){
				x = Domain->xn + i;
				while(x>=NLattice_x){ x-=NLattice_x; }

				if(!_DecompositionGetLattice(x,y,z)){	//If non solid
					VectorField[LatticeIndex*3] = xValues[i*3];
					VectorField[LatticeIndex*3+1] = xValues[i*3+1];
					VectorField[LatticeIndex*3+2] = xValues[i*3+2];
					LatticeIndex++;
				}
				i++;
			}
			*/

			int x = Domain->x0;
			if(x != 0){
				_DecompositionReadVectorFieldValues(InFile,FileLength,x*3,true);
			}

			_DecompositionReadVectorFieldValues(InFile,FileLength,xValues,0,nX*3,true);		//Read in x values

			int i=0;
			while(i!=nX){
				x = Domain->x0 + i;

				if(!_DecompositionGetLattice(x,y,z)){	//If non solid
					VectorField[LatticeIndex].Vx = xValues[i*3];
					VectorField[LatticeIndex].Vy = xValues[i*3+1];
					VectorField[LatticeIndex].Vz = xValues[i*3+2];
					LatticeIndex++;
				}
				i++;
			}

			yi++;
		}

		zi++;
	}

	delete[] xValues;

}

unsigned char _DecompositionGetVoxelSolid(int x, int y, int z){

	//0 = loop boundary conditions, 1 = solid boundary at domain edge

	unsigned char SolidFlags = 0x00;

	if(_DecompositionGetLattice(x,y,z)){
		SolidFlags |= VoxelSolid;
	}

	Coords C;
	GetCoordinates(x,y,z,&C);

	if(C.zp==0 && BoundaryConditionZ==1){					//Positive z neighbour over solid boundary
		SolidFlags |= PositiveZSolid;
	}else{
		if(_DecompositionGetLattice(x,y,C.zp)){
			SolidFlags |= PositiveZSolid;
		}
	}

	if(C.zn==(NLattice_z-1) && BoundaryConditionZ==1){		//Negative z neighbour over solid boundary
		SolidFlags |= NegativeZSolid;
	}else{
		if(_DecompositionGetLattice(x,y,C.zn)){
			SolidFlags |= NegativeZSolid;
		}
	}

	if(C.yp==0 && BoundaryConditionY==1){					//Positive y neighbour over solid boundary
		SolidFlags |= PositiveYSolid;
	}else{
		if(_DecompositionGetLattice(x,C.yp,z)){
			SolidFlags |= PositiveYSolid;
		}
	}

	if(C.yn==(NLattice_y-1) && BoundaryConditionY==1){		//Negative y neighbour over solid boundary
		SolidFlags |= NegativeYSolid;
	}else{
		if(_DecompositionGetLattice(x,C.yn,z)){
			SolidFlags |= NegativeYSolid;
		}
	}

	if(C.xp==0 && BoundaryConditionX==1){					//Positive x neighbour over solid boundary
		SolidFlags |= PositiveXSolid;
	}else{
		if(_DecompositionGetLattice(C.xp,y,z)){
			SolidFlags |= PositiveXSolid;
		}
	}

	if(C.xn==(NLattice_x-1) && BoundaryConditionX==1){		//Negative x neighbour over solid boundary
		SolidFlags |= NegativeXSolid;
	}else{
		if(_DecompositionGetLattice(C.xn,y,z)){
			SolidFlags |= NegativeXSolid;
		}
	}

	return SolidFlags;
}

int ReceiveLatticeDecompositionData(int NProcessors, int ThreadID){

	MPI_Status Stat;

	Domains = new CPUDomain[NProcessors];
	MPI_Recv(Domains, NProcessors*sizeof(CPUDomain), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);	//Receive CPUDomain struct

	int i=0;
	while(i!=NProcessors){
		Domains[i].NeighbourDomainsXn = new int[Domains[i].nNeighbourDomainsXn];
		Domains[i].NeighbourDomainsXp = new int[Domains[i].nNeighbourDomainsXp];
		Domains[i].NeighbourDomainsYn = new int[Domains[i].nNeighbourDomainsYn];
		Domains[i].NeighbourDomainsYp = new int[Domains[i].nNeighbourDomainsYp];
		Domains[i].NeighbourDomainsZn = new int[Domains[i].nNeighbourDomainsZn];
		Domains[i].NeighbourDomainsZp = new int[Domains[i].nNeighbourDomainsZp];

		//Receive neighbouring domain info
		MPI_Recv(Domains[i].NeighbourDomainsXn, sizeof(int)*Domains[i].nNeighbourDomainsXn, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
		MPI_Recv(Domains[i].NeighbourDomainsXp, sizeof(int)*Domains[i].nNeighbourDomainsXp, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
		MPI_Recv(Domains[i].NeighbourDomainsYn, sizeof(int)*Domains[i].nNeighbourDomainsYn, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
		MPI_Recv(Domains[i].NeighbourDomainsYp, sizeof(int)*Domains[i].nNeighbourDomainsYp, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
		MPI_Recv(Domains[i].NeighbourDomainsZn, sizeof(int)*Domains[i].nNeighbourDomainsZn, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
		MPI_Recv(Domains[i].NeighbourDomainsZp, sizeof(int)*Domains[i].nNeighbourDomainsZp, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
	i++;
	}

	memcpy(&ProcessorDomain,&Domains[ThreadID],sizeof(CPUDomain));

	int nX = ProcessorDomain.x1 - ProcessorDomain.x0 + 1;
	int nY = ProcessorDomain.y1 - ProcessorDomain.y0 + 1;
	int nZ = ProcessorDomain.z1 - ProcessorDomain.z0 + 1;
	
	/*
	int nVoxels = (nX+2)*(nY+2)*(nZ+2);
	_CPUDomainWidthX = nX+2;
	_CPUDomainWidthY = nY+2;
	*/

	int nVoxels = nX*nY*nZ;

	_CPUDomainWidthX = nX;
	_CPUDomainWidthY = nY;

	int nNonSolids;	//Obtain number of non-solids in domain
	MPI_Recv(&nNonSolids, 1, MPI_INT, 0, 0,MPI_COMM_WORLD, &Stat);

	LatticeBase = new unsigned int[nVoxels];	//Sparse storage pointer array
	Lattice = new Voxel[nNonSolids+1];			//Non solids voxels

	char MemSuccess = 0;

	if(LatticeBase==0 || Lattice==0){			//Memory allocation failed
		MemSuccess = 1;
		MPI_Send(&MemSuccess, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		return 1;
	}else{
		MPI_Send(&MemSuccess, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	}

	Lattice[0].Solid = VoxelSolid;
	Lattice[0].Vx = 0;
	Lattice[0].Vy = 0;
	Lattice[0].Vz = 0;

	//Receive Geometry
	
	const int dLength = 1048576;		//1MB
	char* DataBuf = new char[dLength];
	int dCount = 0;

	int LatticeIndex = 1;

	while(dCount < nVoxels){
		int l = dLength;
		if(nVoxels-dCount < l){
			l = nVoxels-dCount;
		}

		MPI_Recv(DataBuf, l, MPI_BYTE, 0, 0,MPI_COMM_WORLD, &Stat);

		int ci = 0;
		while(ci!=l){

			if(DataBuf[ci] == VoxelSolid){

				LatticeBase[dCount + ci] = 0;

			}else{

				LatticeBase[dCount + ci] = LatticeIndex;
				Lattice[LatticeIndex].Solid = DataBuf[ci];
				LatticeIndex++;

			}

			ci++;
		}

		dCount += l;
	}

	delete[] DataBuf;

	//Receive Vectorfield

	double* VfBuf = new double[dLength*3];	//24MB

	dCount = 0;

	while(dCount < nNonSolids){
		int l = dLength;
		if(nNonSolids-dCount < l){
			l = nNonSolids - dCount;
		}

		MPI_Recv(VfBuf, l*3, MPI_DOUBLE, 0, 0,MPI_COMM_WORLD, &Stat);

		int ci = 0;
		while(ci!=l){
			Lattice[dCount + ci + 1].Vx = VfBuf[ci*3];
			Lattice[dCount + ci + 1].Vy = VfBuf[ci*3 + 1];
			Lattice[dCount + ci + 1].Vz = VfBuf[ci*3 + 2];
			ci++;
		}

		dCount += l;
	}

	delete[] VfBuf;

	DecompositionInitialiseLattice(NProcessors, ThreadID);

	return 0;
}

void DecompositionInitialiseLattice(int NProcessors, int ThreadID){

	MPI_Status Stat;

	//Initialise lattice velocities

	short z=0;
	while(z!=ProcessorDomain.DomainWidthZ){
	short y=0;
	while(y!=ProcessorDomain.DomainWidthY){
	short x=0;
	while(x!=ProcessorDomain.DomainWidthX){

		if(!IsElementSolid(x,y,z)){

			Voxel* LatticeElement = GetLattice(x,y,z);

			short xL = x;
			short yL = y;
			short zL = z;
			DomainToLattice(xL,yL,zL);

			Coords CL;
			GetCoordinates(xL,yL,zL,&CL);

			if(LatticeElement->Solid&PositiveXSolid || LatticeElement->Solid&VoxelSolid || (x==(ProcessorDomain.DomainWidthX-1) && CL.xp!=ProcessorDomain.x0)){
				LatticeElement->u2 = 0;
			}else{
				int _xp = x + 1;
				if(_xp == ProcessorDomain.DomainWidthX){ _xp = ProcessorDomain.x0; }
				LatticeElement->u2 = (LatticeElement->Vx + GetLattice(_xp,y,z)->Vx)/2;
			}

			if(LatticeElement->Solid&NegativeXSolid || LatticeElement->Solid&VoxelSolid || (x==0 && CL.xn!=ProcessorDomain.x1)){
				LatticeElement->u1 = 0;
			}else{
				int _xn = x - 1;
				if(_xn == -1){ _xn = ProcessorDomain.x1; }
				LatticeElement->u1 = (LatticeElement->Vx + GetLattice(_xn,y,z)->Vx)/2;
			}

			if(LatticeElement->Solid&PositiveYSolid || LatticeElement->Solid&VoxelSolid || (y==(ProcessorDomain.DomainWidthY-1) && CL.yp!=ProcessorDomain.y0)){
				LatticeElement->v2 = 0;
			}else{
				int _yp = y + 1;
				if(_yp == ProcessorDomain.DomainWidthY){ _yp = ProcessorDomain.y0; }
				LatticeElement->v2 = (LatticeElement->Vy + GetLattice(x,_yp,z)->Vy)/2;
			}

			if(LatticeElement->Solid&NegativeYSolid || LatticeElement->Solid&VoxelSolid || (y==0 && CL.yn!=ProcessorDomain.y1)){
				LatticeElement->v1 = 0;
			}else{
				int _yn = y - 1;
				if(_yn == -1){ _yn = ProcessorDomain.y1; }
				LatticeElement->v1 = (LatticeElement->Vy + GetLattice(x,_yn,z)->Vy)/2;
			}

			if(LatticeElement->Solid&PositiveZSolid || LatticeElement->Solid&VoxelSolid || (z==(ProcessorDomain.DomainWidthZ-1) && CL.zp!=ProcessorDomain.z0)){
				LatticeElement->w2 = 0;
			}else{
				int _zp = z + 1;
				if(_zp == ProcessorDomain.DomainWidthZ){ _zp = ProcessorDomain.z0; }
				LatticeElement->w2 = (LatticeElement->Vz + GetLattice(x,y,_zp)->Vz)/2;
			}

			if(LatticeElement->Solid&NegativeZSolid || LatticeElement->Solid&VoxelSolid || (z==0 && CL.zn!=ProcessorDomain.z1)){
				LatticeElement->w1 = 0;
			}else{
				int _zn = z - 1;
				if(_zn == -1){ _zn = ProcessorDomain.z1; }
				LatticeElement->w1 = (LatticeElement->Vz + GetLattice(x,y,_zn)->Vz)/2;
			}

		}

	x++;
	}
	y++;
	}
	z++;
	}

	//Recieve edge velocities from neighbours

	int i=0;
	while(i!=NProcessors){

		int DomainID;

		if(ThreadID==0){
			DomainID = i;

			cout << endl << "Turn of thread " << i << endl;

			int _i = 1;
			while(_i!=NProcessors){
				MPI_Send(&DomainID, 1, MPI_INT, _i, 0, MPI_COMM_WORLD);
				_i++;
			}
		}else{
			MPI_Recv(&DomainID, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);	//Domain which will request neighbours' values [from thread 0]
		}

		if(DomainID == ThreadID){	//Current thread's turn to be master

			short Rect[7];			//Range on lattice thread wants velocities from, {x0,x1,y0,y1,z0,z1,VelocityComponent}

			// X

			int c = 0;
			while(c!=ProcessorDomain.nNeighbourDomainsXn){
				if(ProcessorDomain.NeighbourDomainsXn[c]==ThreadID){
					c++;
					continue;
				}
				CPUDomain* Dmn = &Domains[ProcessorDomain.NeighbourDomainsXn[c]];

				short x0 = ProcessorDomain.x0 - 1;				//x0 = x1 = Domain.x0-1
				if(x0<0){ x0 += NLattice_x; }
				short x1 = x0;

				short y0 = ProcessorDomain.y0;
				short y1 = ProcessorDomain.y1;

				if(Dmn->y0 > y0){
					y0 = Dmn->y0;
				}
				if(Dmn->y1 < y1){
					y1 = Dmn->y1;
				}

				short z0 = ProcessorDomain.z0;
				short z1 = ProcessorDomain.z1;

				if(Dmn->z0 > z0){
					z0 = Dmn->z0;
				}
				if(Dmn->z1 < z1){
					z1 = Dmn->z1;
				}

				Rect[0] = x0;
				Rect[1] = x1;
				Rect[2] = y0;
				Rect[3] = y1;
				Rect[4] = z0;
				Rect[5] = z1;

				cout << "Requesting data from thread " << ProcessorDomain.NeighbourDomainsXn[c] << " of range " << x0 << " -> "  << x1 << ", " << y0 << " -> "  << y1 << ", " << z0 << " -> "  << z1 << endl;

				Rect[6] = 0;	//X component

				MPI_Send(Rect, 7, MPI_SHORT, ProcessorDomain.NeighbourDomainsXn[c], 0, MPI_COMM_WORLD);

				int DatLength = ((int)x1-(int)x0+1)*((int)y1-(int)y0+1)*((int)z1-(int)z0+1);
				double* VDat = new double[DatLength];	//Recieve thread's data

				MPI_Recv(VDat, DatLength, MPI_DOUBLE, ProcessorDomain.NeighbourDomainsXn[c], 0, MPI_COMM_WORLD, &Stat);

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				short _z = z0;
				while(_z<=z1){
				short _y = y0;
				while(_y<=y1){
					int Index = ((int)y1-(int)y0+1)*((int)_z-(int)z0) + ((int)_y-(int)y0);
					if(!IsElementSolid(0,_y,_z)){
						Voxel* LatticeElement = GetLattice(0,_y,_z);
						if(!(LatticeElement->Solid&NegativeXSolid)){
							LatticeElement->u1 = (LatticeElement->Vx + VDat[Index])/2;
						}
					}
				_y++;
				}
				_z++;
				}

				delete[] VDat;

				c++;
			}

			c = 0;
			while(c!=ProcessorDomain.nNeighbourDomainsXp){
				if(ProcessorDomain.NeighbourDomainsXp[c]==ThreadID){
					c++;
					continue;
				}
				CPUDomain* Dmn = &Domains[ProcessorDomain.NeighbourDomainsXp[c]];

				short x0 = ProcessorDomain.x1 + 1;				//x0 = x1 = Domain.x1 + 1
				if(x0>=NLattice_x){ x0 -= NLattice_x; }
				short x1 = x0;

				short y0 = ProcessorDomain.y0;
				short y1 = ProcessorDomain.y1;

				if(Dmn->y0 > y0){
					y0 = Dmn->y0;
				}
				if(Dmn->y1 < y1){
					y1 = Dmn->y1;
				}

				short z0 = ProcessorDomain.z0;
				short z1 = ProcessorDomain.z1;

				if(Dmn->z0 > z0){
					z0 = Dmn->z0;
				}
				if(Dmn->z1 < z1){
					z1 = Dmn->z1;
				}

				Rect[0] = x0;
				Rect[1] = x1;
				Rect[2] = y0;
				Rect[3] = y1;
				Rect[4] = z0;
				Rect[5] = z1;

				cout << "Requesting data from thread " << ProcessorDomain.NeighbourDomainsXp[c] << " of range " << x0 << " -> "  << x1 << ", " << y0 << " -> "  << y1 << ", " << z0 << " -> "  << z1 << endl;

				Rect[6] = 0;	//X component

				MPI_Send(Rect, 7, MPI_SHORT, ProcessorDomain.NeighbourDomainsXp[c], 0, MPI_COMM_WORLD);

				int DatLength = ((int)x1-(int)x0+1)*((int)y1-(int)y0+1)*((int)z1-(int)z0+1);
				double* VDat = new double[DatLength];	//Recieve thread's data

				MPI_Recv(VDat, DatLength, MPI_DOUBLE, ProcessorDomain.NeighbourDomainsXp[c], 0, MPI_COMM_WORLD, &Stat);

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				short _z = z0;
				while(_z<=z1){
				short _y = y0;
				while(_y<=y1){
					int Index = ((int)y1-(int)y0+1)*((int)_z-(int)z0) + ((int)_y-(int)y0);
					if(!IsElementSolid(ProcessorDomain.DomainWidthX-1,_y,_z)){
						Voxel* LatticeElement = GetLattice(ProcessorDomain.DomainWidthX-1,_y,_z);
						if(!(LatticeElement->Solid&PositiveXSolid)){
							LatticeElement->u2 = (LatticeElement->Vx + VDat[Index])/2;
						}
					}
				_y++;
				}
				_z++;
				}

				delete[] VDat;

				c++;
			}

			// Y

			c = 0;
			while(c!=ProcessorDomain.nNeighbourDomainsYn){
				if(ProcessorDomain.NeighbourDomainsYn[c]==ThreadID){
					c++;
					continue;
				}
				CPUDomain* Dmn = &Domains[ProcessorDomain.NeighbourDomainsYn[c]];

				short x0 = ProcessorDomain.x0;
				short x1 = ProcessorDomain.x1;

				if(Dmn->x0 > x0){
					x0 = Dmn->x0;
				}
				if(Dmn->x1 < x1){
					x1 = Dmn->x1;
				}

				short y0 = ProcessorDomain.y0 - 1;				//y0 = y1 = Domain.y0 - 1
				if(y0<0){ y0 += NLattice_y; }
				short y1 = y0;

				short z0 = ProcessorDomain.z0;
				short z1 = ProcessorDomain.z1;

				if(Dmn->z0 > z0){
					z0 = Dmn->z0;
				}
				if(Dmn->z1 < z1){
					z1 = Dmn->z1;
				}

				Rect[0] = x0;
				Rect[1] = x1;
				Rect[2] = y0;
				Rect[3] = y1;
				Rect[4] = z0;
				Rect[5] = z1;

				cout << "Requesting data from thread " << ProcessorDomain.NeighbourDomainsYn[c] << " of range " << x0 << " -> "  << x1 << ", " << y0 << " -> "  << y1 << ", " << z0 << " -> "  << z1 << endl;

				Rect[6] = 1;	//Y component

				MPI_Send(Rect, 7, MPI_SHORT, ProcessorDomain.NeighbourDomainsYn[c], 0, MPI_COMM_WORLD);

				int DatLength = ((int)x1-(int)x0+1)*((int)y1-(int)y0+1)*((int)z1-(int)z0+1);
				double* VDat = new double[DatLength];	//Recieve thread's data

				MPI_Recv(VDat, DatLength, MPI_DOUBLE, ProcessorDomain.NeighbourDomainsYn[c], 0, MPI_COMM_WORLD, &Stat);

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				short _z = z0;
				while(_z<=z1){
				short _x = x0;
				while(_x<=x1){
					int Index = ((int)x1-(int)x0+1)*((int)_z-(int)z0) + ((int)_x-(int)x0);
					if(!IsElementSolid(_x,0,_z)){
						Voxel* LatticeElement = GetLattice(_x,0,_z);
						if(!(LatticeElement->Solid&NegativeYSolid)){
							LatticeElement->v1 = (LatticeElement->Vy + VDat[Index])/2;
						}
					}
				_x++;
				}
				_z++;
				}

				delete[] VDat;

				c++;
			}

			c = 0;
			while(c!=ProcessorDomain.nNeighbourDomainsYp){
				if(ProcessorDomain.NeighbourDomainsYp[c]==ThreadID){
					c++;
					continue;
				}
				CPUDomain* Dmn = &Domains[ProcessorDomain.NeighbourDomainsYp[c]];

				short x0 = ProcessorDomain.x0;
				short x1 = ProcessorDomain.x1;

				if(Dmn->x0 > x0){
					x0 = Dmn->x0;
				}
				if(Dmn->x1 < x1){
					x1 = Dmn->x1;
				}

				short y0 = ProcessorDomain.y1 + 1;				//y0 = y1 = Domain.y0 - 1
				if(y0>=NLattice_y){ y0 -= NLattice_y; }
				short y1 = y0;

				short z0 = ProcessorDomain.z0;
				short z1 = ProcessorDomain.z1;

				if(Dmn->z0 > z0){
					z0 = Dmn->z0;
				}
				if(Dmn->z1 < z1){
					z1 = Dmn->z1;
				}

				Rect[0] = x0;
				Rect[1] = x1;
				Rect[2] = y0;
				Rect[3] = y1;
				Rect[4] = z0;
				Rect[5] = z1;

				cout << "Requesting data from thread " << ProcessorDomain.NeighbourDomainsYp[c] << " of range " << x0 << " -> "  << x1 << ", " << y0 << " -> "  << y1 << ", " << z0 << " -> "  << z1 << endl;

				Rect[6] = 1;	//Y component

				MPI_Send(Rect, 7, MPI_SHORT, ProcessorDomain.NeighbourDomainsYp[c], 0, MPI_COMM_WORLD);

				int DatLength = ((int)x1-(int)x0+1)*((int)y1-(int)y0+1)*((int)z1-(int)z0+1);
				double* VDat = new double[DatLength];	//Recieve thread's data

				MPI_Recv(VDat, DatLength, MPI_DOUBLE, ProcessorDomain.NeighbourDomainsYp[c], 0, MPI_COMM_WORLD, &Stat);

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				short _z = z0;
				while(_z<=z1){
				short _x = x0;
				while(_x<=x1){
					int Index = ((int)x1-(int)x0+1)*((int)_z-(int)z0) + ((int)_x-(int)x0);
					if(!IsElementSolid(_x,ProcessorDomain.DomainWidthY-1,_z)){
						Voxel* LatticeElement = GetLattice(_x,ProcessorDomain.DomainWidthY-1,_z);
						if(!(LatticeElement->Solid&PositiveYSolid)){
							LatticeElement->v2 = (LatticeElement->Vy + VDat[Index])/2;
						}
					}
				_x++;
				}
				_z++;
				}

				delete[] VDat;

				c++;
			}

			// Z

			c = 0;
			while(c!=ProcessorDomain.nNeighbourDomainsZn){
				if(ProcessorDomain.NeighbourDomainsZn[c]==ThreadID){
					c++;
					continue;
				}
				CPUDomain* Dmn = &Domains[ProcessorDomain.NeighbourDomainsZn[c]];

				short x0 = ProcessorDomain.x0;
				short x1 = ProcessorDomain.x1;

				if(Dmn->x0 > x0){
					x0 = Dmn->x0;
				}
				if(Dmn->x1 < x1){
					x1 = Dmn->x1;
				}

				short y0 = ProcessorDomain.y0;
				short y1 = ProcessorDomain.y1;

				if(Dmn->y0 > y0){
					y0 = Dmn->y0;
				}
				if(Dmn->y1 < y1){
					y1 = Dmn->y1;
				}

				short z0 = ProcessorDomain.z0 - 1;				//z0 = z1 = Domain.z0-1
				if(z0<0){ z0 += NLattice_z; }
				short z1 = z0;

				Rect[0] = x0;
				Rect[1] = x1;
				Rect[2] = y0;
				Rect[3] = y1;
				Rect[4] = z0;
				Rect[5] = z1;

				cout << "Requesting data from thread " << ProcessorDomain.NeighbourDomainsZn[c] << " of range " << x0 << " -> "  << x1 << ", " << y0 << " -> "  << y1 << ", " << z0 << " -> "  << z1 << endl;

				Rect[6] = 2;	//Z component

				MPI_Send(Rect, 7, MPI_SHORT, ProcessorDomain.NeighbourDomainsZn[c], 0, MPI_COMM_WORLD);

				int DatLength = ((int)x1-(int)x0+1)*((int)y1-(int)y0+1)*((int)z1-(int)z0+1);
				double* VDat = new double[DatLength];	//Recieve thread's data

				MPI_Recv(VDat, DatLength, MPI_DOUBLE, ProcessorDomain.NeighbourDomainsZn[c], 0, MPI_COMM_WORLD, &Stat);

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				short _y = y0;
				while(_y<=y1){
				short _x = x0;
				while(_x<=x1){
					int Index = ((int)x1-(int)x0+1)*((int)_y-(int)y0) + ((int)_x-(int)x0);
					if(!IsElementSolid(_x,_y,0)){
						Voxel* LatticeElement = GetLattice(_x,_y,0);
						if(!(LatticeElement->Solid&NegativeZSolid)){
							LatticeElement->w1 = (LatticeElement->Vz + VDat[Index])/2;
						}
					}
				_x++;
				}
				_y++;
				}

				delete[] VDat;

				c++;
			}

			c = 0;
			while(c!=ProcessorDomain.nNeighbourDomainsZp){
				if(ProcessorDomain.NeighbourDomainsZp[c]==ThreadID){
					c++;
					continue;
				}
				CPUDomain* Dmn = &Domains[ProcessorDomain.NeighbourDomainsZp[c]];

				short x0 = ProcessorDomain.x0;
				short x1 = ProcessorDomain.x1;

				if(Dmn->x0 > x0){
					x0 = Dmn->x0;
				}
				if(Dmn->x1 < x1){
					x1 = Dmn->x1;
				}

				short y0 = ProcessorDomain.y0;
				short y1 = ProcessorDomain.y1;

				if(Dmn->y0 > y0){
					y0 = Dmn->y0;
				}
				if(Dmn->y1 < y1){
					y1 = Dmn->y1;
				}

				short z0 = ProcessorDomain.z1 + 1;				//z0 = z1 = Domain.z1+1
				if(z0>=NLattice_z){ z0 -= NLattice_z; }
				short z1 = z0;

				Rect[0] = x0;
				Rect[1] = x1;
				Rect[2] = y0;
				Rect[3] = y1;
				Rect[4] = z0;
				Rect[5] = z1;

				cout << "Requesting data from thread " << ProcessorDomain.NeighbourDomainsZp[c] << " of range " << x0 << " -> "  << x1 << ", " << y0 << " -> "  << y1 << ", " << z0 << " -> "  << z1 << endl;

				Rect[6] = 2;	//Z component

				MPI_Send(Rect, 7, MPI_SHORT, ProcessorDomain.NeighbourDomainsZp[c], 0, MPI_COMM_WORLD);

				int DatLength = ((int)x1-(int)x0+1)*((int)y1-(int)y0+1)*((int)z1-(int)z0+1);
				double* VDat = new double[DatLength];	//Recieve thread's data

				MPI_Recv(VDat, DatLength, MPI_DOUBLE, ProcessorDomain.NeighbourDomainsZp[c], 0, MPI_COMM_WORLD, &Stat);

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				short _y = y0;
				while(_y<=y1){
				short _x = x0;
				while(_x<=x1){
					int Index = ((int)x1-(int)x0+1)*((int)_y-(int)y0) + ((int)_x-(int)x0);
					if(!IsElementSolid(_x,_y,ProcessorDomain.DomainWidthZ-1)){
						Voxel* LatticeElement = GetLattice(_x,_y,ProcessorDomain.DomainWidthZ-1);
						if(!(LatticeElement->Solid&PositiveZSolid)){
							LatticeElement->w2 = (LatticeElement->Vz + VDat[Index])/2;
						}
					}
				_x++;
				}
				_y++;
				}

				delete[] VDat;

				c++;
			}

			Rect[6] = -1;			//Inform all threads that the current thread is done
			int _i = 0;
			while(_i!=NProcessors){
				if(_i!=ThreadID){
				MPI_Send(Rect, 7, MPI_SHORT, _i, 0, MPI_COMM_WORLD);
				}
				_i++;
			}

		}else{						//Answer any requests

			while(true){

				short Rect[7];		//Range on lattice thread wants velocities from, {x0,x1,y0,y1,z0,z1,VelocityComponent}
				MPI_Recv(Rect, 7, MPI_SHORT, DomainID, 0, MPI_COMM_WORLD, &Stat);

				if(Rect[6]==-1){	//Finished or no data needed from current thread
					break;
				}

				short x0 = Rect[0];		//Index of first x in range
				short x1 = Rect[1];		//Index of last (inclusive) x in range
				short y0 = Rect[2];
				short y1 = Rect[3];
				short z0 = Rect[4];
				short z1 = Rect[5];

				int DatLength = ((int)x1-(int)x0+1)*((int)y1-(int)y0+1)*((int)z1-(int)z0+1);
				double* VDat = new double[DatLength];	//Data to send back

				LatticeToDomain(x0,y0,z0);
				LatticeToDomain(x1,y1,z1);

				cout << "Thread " << ThreadID << " recieved request of domain coordinates: "  << x0 << " -> "  << x1 << ", " << y0 << " -> "  << y1 << ", " << z0 << " -> "  << z1 << endl;

				short _z = z0;
				while(_z<=z1){
				short _y = y0;
				while(_y<=y1){
				short _x = x0;
				while(_x<=x1){
					int Index = ((int)x1-(int)x0+1)*(((int)y1-(int)y0+1)*((int)_z-(int)z0) + ((int)_y-(int)y0)) + ((int)_x-(int)x0);
					if(IsElementSolid(_x,_y,_z)){
						VDat[Index] = 0;
					}else{
						if(Rect[6]==0){
							VDat[Index] = GetLattice(_x,_y,_z)->Vx;
						}else if(Rect[6]==1){
							VDat[Index] = GetLattice(_x,_y,_z)->Vy;
						}else{
							VDat[Index] = GetLattice(_x,_y,_z)->Vz;
						}
					}
				_x++;
				}
				_y++;
				}
				_z++;
				}

				MPI_Send(VDat, DatLength, MPI_DOUBLE, DomainID, 0, MPI_COMM_WORLD);

				delete[] VDat;

			}

		}

		i++;
	}

	//Exchange lattice info (porosity etc)

	LatticeInformation DomainLatticeInfo;
	
	double VMax = 0;
	Vector V = {0,0,0};
	int ncount=0;
	
	z=0;
	while(z!=ProcessorDomain.DomainWidthZ){
	short y=0;
	while(y!=ProcessorDomain.DomainWidthY){
	short x=0;
	while(x!=ProcessorDomain.DomainWidthX){

		Voxel* LatticeElement = GetLattice(x,y,z);

		if(!(LatticeElement->Solid&VoxelSolid)){
			V.X += LatticeElement->Vx;
			V.Y += LatticeElement->Vy;
			V.Z += LatticeElement->Vz;
			ncount++;

			if(fabs(LatticeElement->Vx) > VMax){
				VMax = fabs(LatticeElement->Vx);
			}
			if(fabs(LatticeElement->Vy) > VMax){
				VMax = fabs(LatticeElement->Vy);
			}
			if(fabs(LatticeElement->Vz) > VMax){
				VMax = fabs(LatticeElement->Vz);
			}
		}

	x++;
	}
	y++;
	}
	z++;
	}

	V.X /= ((double)ncount);
	V.Y /= ((double)ncount);
	V.Z /= ((double)ncount);

	DomainLatticeInfo.AverageFlowVelocity.X = V.X;
	DomainLatticeInfo.AverageFlowVelocity.Y = V.Y;
	DomainLatticeInfo.AverageFlowVelocity.Z = V.Z;

	DomainLatticeInfo.MaxFlowVelocityComponent = VMax;

	DomainLatticeInfo.NSolids = (ProcessorDomain.DomainWidthX*ProcessorDomain.DomainWidthY*ProcessorDomain.DomainWidthZ) - ncount;
	DomainLatticeInfo.NNonSolids = ncount;

	DomainLatticeInfo.Porosity = ((double)ncount)/((double)(ProcessorDomain.DomainWidthX*ProcessorDomain.DomainWidthY*ProcessorDomain.DomainWidthZ));

	ProcessorDomain.LatticeInfo = DomainLatticeInfo;					//Processor domain's lattice info

	if(ThreadID==0){

		LatticeInformation* DomainLatticeInfos = new LatticeInformation[NProcessors];

		memcpy(&DomainLatticeInfos[0], &DomainLatticeInfo, sizeof(LatticeInformation));

		Domains[0].LatticeInfo = DomainLatticeInfos[0];
		ProcessorDomain.LatticeInfo = DomainLatticeInfos[0];

		double VMax_Lattice = DomainLatticeInfos[0].MaxFlowVelocityComponent;		//Max flow velocity
		Vector V_Lattice;															//Average flow velocity
		V_Lattice.X = DomainLatticeInfos[0].AverageFlowVelocity.X;
		V_Lattice.Y = DomainLatticeInfos[0].AverageFlowVelocity.Y;
		V_Lattice.Z = DomainLatticeInfos[0].AverageFlowVelocity.Z;
		int ncount_Lattice = DomainLatticeInfos[0].NNonSolids;						//Number of non solids
		int scount_Lattice = DomainLatticeInfos[0].NSolids;							//Number of solids

		i = 1;
		while(i!=NProcessors){
			MPI_Recv(&DomainLatticeInfos[i], sizeof(LatticeInformation), MPI_BYTE, i, 0, MPI_COMM_WORLD, &Stat);

			Domains[i].LatticeInfo = DomainLatticeInfos[i];

			if(DomainLatticeInfos[i].MaxFlowVelocityComponent > VMax_Lattice){
				VMax_Lattice = DomainLatticeInfos[i].MaxFlowVelocityComponent;
			}

			V_Lattice.X += DomainLatticeInfos[i].AverageFlowVelocity.X;
			V_Lattice.Y += DomainLatticeInfos[i].AverageFlowVelocity.Y;
			V_Lattice.Z += DomainLatticeInfos[i].AverageFlowVelocity.Z;

			ncount_Lattice += DomainLatticeInfos[i].NNonSolids;
			scount_Lattice += DomainLatticeInfos[i].NSolids;		

			i++;
		}

		//Global LatticeInfo struct

		LatticeInfo.AverageFlowVelocity.X = V_Lattice.X / ((double)NProcessors);
		LatticeInfo.AverageFlowVelocity.Y = V_Lattice.Y / ((double)NProcessors);
		LatticeInfo.AverageFlowVelocity.Z = V_Lattice.Z / ((double)NProcessors);

		LatticeInfo.MaxFlowVelocityComponent = VMax_Lattice;

		LatticeInfo.NNonSolids = ncount_Lattice;
		LatticeInfo.NSolids = scount_Lattice;

		i = 1;
		while(i!=NProcessors){
			MPI_Send(&LatticeInfo, sizeof(LatticeInformation), MPI_BYTE, i, 0, MPI_COMM_WORLD);

			int _i = 0;
			while(_i!=NProcessors){
				MPI_Send(&DomainLatticeInfos[_i], sizeof(LatticeInformation), MPI_BYTE, i, 0, MPI_COMM_WORLD);
				_i++;
			}

			i++;
		}

		delete[] DomainLatticeInfos;

		cout << endl << "Exchanged local domain information" << endl;

	}else{

		//Send local lattice information
		MPI_Send(&DomainLatticeInfo, sizeof(LatticeInformation), MPI_BYTE, 0, 0, MPI_COMM_WORLD);

		//Recieve global lattice information
		MPI_Recv(&LatticeInfo, sizeof(LatticeInformation), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);

		//Recieve domains' lattice infos
		int _i = 0;
		while(_i!=NProcessors){
			MPI_Recv(&(Domains[_i].LatticeInfo), sizeof(LatticeInformation), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &Stat);
			_i++;
		}

	}

}

void DomainInitialiseParticlesUniform(int NProcessors, int ThreadID, int* NParticlesThread){

	if(ParticleNum == 0){
		NParticles = LatticeInfo.NNonSolids * ParticlesPerVoxel;

		int i=0;
		while(i!=NProcessors){
			Domains[i].NParticlesDomain = Domains[i].LatticeInfo.NNonSolids * ParticlesPerVoxel;
			i++;
		}
	}else{
		NParticles = ParticleNum;
		
		int PCount = 0;
		int i=0;
		while(i!=NProcessors){
			Domains[i].NParticlesDomain = (int)floor( ((double)Domains[i].LatticeInfo.NNonSolids / (double)LatticeInfo.NNonSolids) * (double)NParticles );
			PCount += Domains[i].NParticlesDomain;
			i++;
		}

		i=0;
		while(PCount!=0){
			Domains[i].NParticlesDomain++;
			PCount--;
			i++;
			if(i==NProcessors){
				i = 0;
			}
		}
	}

	ProcessorDomain.NParticlesDomain = Domains[ThreadID].NParticlesDomain;
	*NParticlesThread = Domains[ThreadID].NParticlesDomain;

	if(ThreadID == 0){
		cout << "Number of particles: " << NParticles << endl;

		int i=0;
		while(i!=NProcessors){
			cout << "Thread " << i << " begins with " << Domains[i].NParticlesDomain << " particles" << endl;
			i++;
		}

		cout << endl;
	}

	ParticleArrayLength = (int)ceil( (double)(*NParticlesThread) * 1.10 );		//Allocate particle array with 10% extra space

	Particles = new Particle[ParticleArrayLength];

	if(Particles==0){
		cout << "Memory allocation failure: Thread [" << ThreadID << "] particle array" << endl;
		return;
	}

	double pSpacing = pow( (double)(ProcessorDomain.LatticeInfo.NNonSolids) , 1.0/3.0 ) / pow( (double)(*NParticlesThread) , 1.0/3.0 );	//Distance between particles
	
	double Px = pSpacing/2;
	double Py = pSpacing/2;
	double Pz = pSpacing/2;

	int i=0;
	while(i!=(*NParticlesThread)){

		Voxel* LatticeElement = GetLattice((int)floor(Px),(int)floor(Py),(int)floor(Pz));
				
		if(!(LatticeElement->Solid&VoxelSolid)){

			Particles[i].X = Px;
			Particles[i].Y = Py;
			Particles[i].Z = Pz;

			Particles[i].VoxelX = (short)floor(Particles[i].X);
			Particles[i].VoxelY = (short)floor(Particles[i].Y);
			Particles[i].VoxelZ = (short)floor(Particles[i].Z);

			double InitX = Particles[i].X;
			double InitY = Particles[i].Y;
			double InitZ = Particles[i].Z;

			DomainToLattice(InitX,InitY,InitZ);

			Particles[i].InitialX = InitX;			//Set initial positions on lattice rather than domain
			Particles[i].InitialY = InitY;
			Particles[i].InitialZ = InitZ;

			Particles[i].LastVoxelX = Particles[i].VoxelX;
			Particles[i].LastVoxelY = Particles[i].VoxelY;
			Particles[i].LastVoxelZ = Particles[i].VoxelZ;

			Particles[i].XCycles = 0;
			Particles[i].YCycles = 0;
			Particles[i].ZCycles = 0;

			i++;
		}

		Px += pSpacing;

		if(Px >= (double)ProcessorDomain.DomainWidthX){

			Px = pSpacing/2;
			Py += pSpacing;

			if(Py >= (double)ProcessorDomain.DomainWidthY){

				Py = pSpacing/2;
				Pz += pSpacing;

				if(Pz >= (double)ProcessorDomain.DomainWidthZ){

					Pz = pSpacing/2;

				}
			}
		}

	}

}

#endif