//Compile options
#define CompileMPI					1	//Compile with MPI functions

//#define AdvectionMode				1	//0 = Pollock algorithm, 1 = Predictor-Corrector

#ifdef WIN32
#define WindowsLinux				0	//0 for Windows, 1 for Linux
#else
#define WindowsLinux				1	//0 for Windows, 1 for Linux
#endif 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Shared Mem Multithreading options
int MultiThreadingSharedMem_TNum = 7;		//Number of threads to run in shared memory mode

//Lattice
short NLattice_x = 100;						//Lattice of NLattice_x x NLattice_y x NLattice_z
short NLattice_y = 100;
short NLattice_z = 100;

int BoundaryConditionX = 0;					//Boundary conditions in X,Y,Z. Loop boundary = 0, Solid boundary = 1
int BoundaryConditionY = 0;
int BoundaryConditionZ = 0;

double LatticeResolutionMicron = 4.9;		//Resolution of lattice in micrometres

int ParticleNum = 0;						//Number of tracer particles. Set to 0 to use particles per voxel...
int ParticlesPerVoxel = 1;					//For uniform initialisation, number of particles per voxel. Use a cube number eg 1, 8, 27, 64...

double DiffusionCoefficient = 0.00001;		//Diffusion coefficient for random walk length

double SimulationTimestep = 0;				//Timestep for the simulation. Set to 0 for automatic (optimum) value
double SimulationTimeScale = 1000000;		//Value of dt corresponding physically to 1 second

double OutputTimeInterval = 10000;		//Time interval for writing output data
double SimulationTimeMax =  0;			//End time of simulation. Set to 0 for no limit

//Input folder
char InputFilesFolder[1024]		= "H:\\2012 ChemEng PhD\\VC++ Projects\\Input Files\\BeadElements\\";

bool SolidsFileBin = false;			//Is solids file binary format
bool VelocitiesFileBin = false;		//Is velocities file binary format
bool VelocitiesFileSparse = false;	//Sparse velocity files

//Input files
char SolidsFileName[1024]		= "BeadElementSolids_R=15.txt";			//File specifying solid voxels
char VelocitiesFileName[1024]	= "BeadElementVelocity_R=15.dat";		//File containing vector field

//Output folder
char OutputFilesFolder[1024]	= "G:\\Output Files\\";

//Output files
char GeometryOutputFile[1024]	= "Geometry t=%ST.vtk";
char DataOutputFile[1024]		= "Data.txt";
char AdvectionDebugFile[1024]	= "Advection Debug.txt";		//Advection debug output. Use %ID for thread ID

//LB-Disp Working Directory
char WorkingFolder[1024]		= "G:\\Working Directory\\";

//Working Directory Files
char SimulationOutputFile[1024]	= "SimulationOutput.dat";	//Output simulation variables and particles
char GeometryUpdateFile[1024]	= "GeometryUpdate.txt";
char WorkingGeometry[1024]		= "Geometry.dat";
char WorkingVelocities[1024]	= "Velocity.dat";

int NAlterationsEnd = 0;			//Number of voxels changed before stopping for velocity update step

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SimulationContinue = true;


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
#include <ctype.h>

using namespace std;

//////////////////////////////

//Advection Debug
//Debug mode:	0 = off
//				1 = minimal (catches infinite loops only - fastest)
//				2 = detailed (captures infinite loops and infs and NaNs and generates debug data - fast but rhobust)
//				3 = highest (captures all of above + infs and NaNs in tau values - fast but not necessary)

#define DebugAdvectionMode	2		//Debug mode
#define OutputDebugData		false	//Whether to write to output file

bool DebugParticle = false;
bool TauError = false;

ofstream AdvectDebug;

inline bool isInfNaN(double x){
	return (!(x <= DBL_MAX && x >= -DBL_MAX));
}

//Console output
#define ConsoleSpacer "-----------------------------------------------------------"

//////////////////////////////

#if WindowsLinux==0		//Windows headers

#include <Windows.h>
#include <process.h>

#else					//Linux headers

#define _snprintf snprintf
#define Sleep usleep

#endif

#define PI 3.14159265358979323846264338328
#define Min_dV (1.0E-9)							//Minimum difference between opposing face velocities for which equations hold (otherwise treat as special cases)
#define MIN_VECTOR_RESOLUTION Min_dV			//Vectors read in smaller than this are set to 0

int NParticles;								//Number of particles

struct Voxel{				//sizeof(Voxel) = 88
double Vx;					//Voxel vector x
double Vy;					//Voxel vector y
double Vz;					//Voxel vector z
double u1;					//Voxel face velocity x1
double u2;					//Voxel face velocity x2
double v1;					//Voxel face velocity y1
double v2;					//Voxel face velocity y2
double w1;					//Voxel face velocity z1
double w2;					//Voxel face velocity z2 
double DissolutionFraction;	//0 = whole, 1 = completely dissolved
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

struct Particle{			//sizeof(Particle) = 80 (73 unpadded)
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
char Type;					//Particle type (index in ParticleTypes array); 0 = unreactive
};

struct ParticleType{
	double RandomWalkLength;		//Random walk length for timestep
	double DiffusionCoefficient;	//Particle type diffusion coefficient
	double DissolutionFraction;		//Fraction of a solid voxel the particle dissolves
	double Mol;						//Number of mols the particle represents
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

struct Vector{
double X;
double Y;
double Z;
};

struct LatticeInformation{
	Vector AverageFlowVelocity;
	double MaxFlowVelocityComponent;
	int NNonSolids;
	int NSolids;
	int NVoxels;
	double Porosity;
};

//////////////////////////////////////////

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
		
		long double Tfr2 = ((long double)UINT_MAX)*(ThreadID + 1.0 + TFraction)/(double)(ThreadNum + 1.1);

		m_z = (unsigned int)floor(Tfr2);

		if(m_z==0x9068FFFF){		//z = 0x9068FFFF fixed point
			m_z = 362436069;
		}

		GetUniform();

	}

};

Voxel* GetLattice(int x, int y, int z);
char* WorkingFilePath(char* FileName);

class GeometryAlterationLog{

private:
	struct VoxelAlteration{
		short Alteration;
		short X;
		short Y;
		short Z;
	};

	static const int InitialArraySize = 1000;

	VoxelAlteration** Alterations;
	int* AlterationsLength;
	int* AlterationsIndex;

	int T_Num;

	int NAlterations;

public:

	GeometryAlterationLog(int ThreadNum){

		NAlterations = 0;

		Alterations = new VoxelAlteration*[ThreadNum];
		AlterationsLength = new int[ThreadNum];
		AlterationsIndex = new int[ThreadNum];

		int i=0;
		while(i!=ThreadNum){

			Alterations[i] = new VoxelAlteration[InitialArraySize];

			AlterationsLength[i] = InitialArraySize;
			AlterationsIndex[i] = 0;

			i++;
		}

		T_Num = ThreadNum;
	}

	void LogAlteration(short x, short y, short z, short Alteration, int ThreadID){ 

		NAlterations++;

		int Index = AlterationsIndex[ThreadID];

		AlterationsIndex[ThreadID]++;

		if(Index == AlterationsLength[ThreadID]){
			ExtendArray(ThreadID);
		}

		Alterations[ThreadID][Index].Alteration = Alteration;
		Alterations[ThreadID][Index].X = x;
		Alterations[ThreadID][Index].Y = y;
		Alterations[ThreadID][Index].Z = z;

		if(NAlterations >= NAlterationsEnd){
			SimulationContinue = false;
		}

	}

	bool ExtendArray(int ID){

		int newSize = (int)( (double)(AlterationsLength[ID]) * 1.2 );		//Increase by 20%

		VoxelAlteration* newArr = new VoxelAlteration[newSize];

		if(newArr==0){
			return false;
		}

		memcpy( newArr, Alterations[ID], sizeof(VoxelAlteration)*AlterationsLength[ID] );	//Copy old array into new

		delete[] Alterations[ID];

		Alterations[ID] = newArr;

		AlterationsLength[ID] = newSize;

		return true;
	}

	void OutputAlterations(){

		ofstream outfile;
		outfile.open(WorkingFilePath(GeometryUpdateFile),ios_base::trunc);

		cout << "Outputting to " << WorkingFilePath(GeometryUpdateFile) << endl;

		int i=0;
		while(i!=T_Num){

			int c=0;
			while(c!=AlterationsIndex[i]){

				Voxel* Vx = GetLattice(Alterations[i][c].X, Alterations[i][c].Y, Alterations[i][c].Z);

				outfile << Alterations[i][c].Alteration << ' ';
				outfile << Alterations[i][c].X << ' ' << Alterations[i][c].Y << ' ' << Alterations[i][c].Z << ' ';
				outfile << Vx->Vx << ' ' << Vx->Vy << ' ' << Vx->Vz << endl;

				c++;
			}

			i++;
		}

		outfile.close();

	}

	~GeometryAlterationLog(){

		int i=0;
		while(i!=T_Num){

			delete[] Alterations[i];

			i++;
		}

		delete[] Alterations;
		delete[] AlterationsLength;
		delete[] AlterationsIndex;

	}

};

RandomNumber Random;								//Initialise class called Random. Get next number with: Random.GetUniform();

//////////////////////////////////////////

char* InputFilePath(char* FileName);
char* OutputFilePath(char* FileName);
char* OutputFilePath(char* FileName, double t);

void CreateLattice();													//Sets lattice entirely non solid, all flow vectors to 0
void CreateWorkingLattice();
void InitialiseLattice();												//Calculates and sets voxel face velocities
void InterpolateVelocity(int SolidX, int SolidY, int SolidZ, Voxel* Vxl);	//Interpolates velocity of a voxel from its neighbours and sets face velocities
void ReadSolids(unsigned int* NumberOfNonSolids, char* FilePath);
void ReadSolids(unsigned int* NumberOfNonSolids);						//Reads in the solids file for non-decomposed lattice
void ReadSolidsBin(unsigned int* NumberOfNonSolids, char* FilePath);
void ReadSolidsBin(unsigned int* NumberOfNonSolids);					//Reads in a binary solids file for non-decomposed lattice
void ReadVelocities(char* FilePath);
void ReadVelocities();													//Sets lattice velocity field from file
void ReadVelocitiesSparse(char* FilePath);
void ReadVelocitiesSparse();
void ReadVelocitiesBin(char* FilePath);
void ReadVelocitiesBin();
void ReadVelocitiesBinSparse(char* FilePath);
void ReadVelocitiesBinSparse();											//Sets lattice velocity field from binary format file
double GetTimestep();													//Calculate optimum timestep
void GetCoordinates(short x, short y, short z, Coords* C);
void GetCoordinates(int i, Coords* C);
void CycleCoordinates(double* Px, double* Py, double* Pz);
void OutputGeometry(int n);												//Output geometry file
void InitialiseParticles();												//Called from main to set initial particle positions
void InitialiseParticlesUniform();										//Disitrbutes particles evenly throughout the lattice
void SetLatticeInfo();													//Calculates average velocities, porosity etc in LatticeInfo struct
void WriteOutParameters();												//Writes out average velocities, porosity, Peclet number to the console
void RunTimeStep(double dt, int ThreadID);								//Carries out an advection and diffusion step on all particles
void RandomWalkParticle(int i, double dt, double RandomWalkLength, int ThreadID);
void VoxelAfterDisplacement(Particle* P, Vector* V, Coords* C);
void RandomWalkVector(double L, Vector* V);								//Returns a random walk vector of length L
void ParticleSolidReaction(Particle* P, int SolidX, int SolidY, int SolidZ, bool* ParticleRemoved, int ThreadID);		//Carry out solid-particle reaction
void AdvectParticle(Particle* P, double dt, int ThreadID);
void CaseNoSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt);
void CaseOneSolid(Particle* P, Coords* C, Voxel* LatticeElement, double* dt);
void CaseTwoSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt);
void CaseThreeSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt);
void CaseFourSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt);

//////////////////////////////////////////

Voxel* Lattice;							//Array of lattice elements. Element 0 is a solid voxel

int _CPUDomainWidthX;					//Number of X elements in array
int _CPUDomainWidthY;					//Number of Y elements in array

inline bool IsElementSolid(short x, short y, short z){
	return (Lattice[_CPUDomainWidthX*(_CPUDomainWidthY*z + y) + x].Solid&VoxelSolid);
}

inline Voxel* GetLattice(int x, int y, int z){
	return &(Lattice[_CPUDomainWidthX*(_CPUDomainWidthY*z + y) + x]);
}

ParticleType* ParticleTypes;		//Types of particle
Particle* Particles;				//Dynamic array to allocate number of particles assigned to each thread

LatticeInformation LatticeInfo;		//Struct to be filled with (full) lattice information

double AdvectionThreshold;			//Minimum velocity to carry out advection calculation

GeometryAlterationLog* AlterationLog;	//Record alterations to the geometry

//////////////////////////////////////////

char* WorkingFilePath(char* FileName){

	char* Folder = WorkingFolder;

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

struct SimulationStateStruct{
	int ThreadNum;
	int NParticles;
	double SimulationTime;			//t
	double SimulationTimeStep;		//dt
	double SimulationTimeMax;		//tmax
	int LogInterval;				//LogInterval
	int TimeSteps;					//TimeSteps
	int GeometryOutputIndex;		//n
};

struct PartialVoxel{
	double Fraction;
	short X;
	short Y;
	short Z;
};

void OutputSimulationState(SimulationStateStruct* State){

	FILE* OutFile = fopen(WorkingFilePath(SimulationOutputFile), "wb");

	fwrite(State, sizeof(SimulationStateStruct), 1, OutFile);		//Write state variables to file

	fwrite(Particles, sizeof(Particle), NParticles, OutFile);		//Write out particle array

	//Find number of voxels partially dissolved

	int nPartialVoxels = 0;

	short x=0;
	while(x!=NLattice_x){
	short y=0;
	while(y!=NLattice_y){
	short z=0;
	while(z!=NLattice_z){

		Voxel* LatticeElement = GetLattice(x,y,z);

		if((LatticeElement->Solid&VoxelSolid) && LatticeElement->DissolutionFraction != 0){
			nPartialVoxels++;
		}

	z++;
	}
	y++;
	}
	x++;
	}
	
	fwrite(&nPartialVoxels, sizeof(int), 1, OutFile);

	PartialVoxel* PartialVx = new PartialVoxel[nPartialVoxels];

	int c = 0;
	x=0;
	while(x!=NLattice_x){
	short y=0;
	while(y!=NLattice_y){
	short z=0;
	while(z!=NLattice_z){

		Voxel* LatticeElement = GetLattice(x,y,z);

		if((LatticeElement->Solid&VoxelSolid) && LatticeElement->DissolutionFraction != 0){
			PartialVx[c].X = x;
			PartialVx[c].Y = y;
			PartialVx[c].Z = z;
			PartialVx[c].Fraction = LatticeElement->DissolutionFraction;
			c++;
		}

	z++;
	}
	y++;
	}
	x++;
	}

	fwrite(PartialVx, sizeof(PartialVoxel), nPartialVoxels, OutFile);
	
	delete[] PartialVx;
	fclose(OutFile);
}

bool ReadSimulationState(SimulationStateStruct* State){

	FILE* InFile = fopen(WorkingFilePath(SimulationOutputFile), "rb");

	if(InFile == 0){	//Failed to open file
		return false;
	}

	//Read state variables

	fread(State, sizeof(SimulationStateStruct), 1, InFile);		
	 
	NParticles = State->NParticles;

	//Create particles and read in from file

	Particles = new Particle[NParticles];				

	fread(Particles, sizeof(Particle), State->NParticles, InFile);

	//Create lattice from working directory

	CreateWorkingLattice();		//Create lattice from geometry and velocity in working directory

	InitialiseLattice();		//Calculate border velocities

	SetLatticeInfo();			//Fills in global LatticeInfo struct

	//Read in partially dissolved voxels

	int nPartialVoxels;
	fread(&nPartialVoxels, sizeof(int), 1, InFile);

	PartialVoxel* PartialVx = new PartialVoxel[nPartialVoxels];

	fread(PartialVx, sizeof(PartialVoxel), nPartialVoxels, InFile);

	int i=0;
	while(i!=nPartialVoxels){
		GetLattice(PartialVx[i].X, PartialVx[i].Y, PartialVx[i].Z)->DissolutionFraction = PartialVx[i].Fraction;
		i++;
	}

	delete[] PartialVx;

	//

	fclose(InFile);

	//Initialise particle types
	
	ParticleTypes = new ParticleType[2];

	//Non-reactive particle
	ParticleTypes[0].DissolutionFraction = 0;
	ParticleTypes[0].Mol = 0;

	//Reactive particle
	ParticleTypes[1].DissolutionFraction = 0.5;
	ParticleTypes[1].Mol = 0.1;

	return true;
}

void ReadSolidsBin(unsigned int* NumberOfNonSolids){

	ReadSolidsBin(NumberOfNonSolids, InputFilePath(SolidsFileName));

}

void ReadSolidsBin(unsigned int* NumberOfNonSolids, char* FilePath){

	FILE* InFile = fopen(FilePath, "rb");
	
	if(InFile!=0){
		cout << "Reading in solids from '" << FilePath << "'" << endl;
	}else{
		cout << "Could not open solids file" << endl;
		return;
	}

	__int64 flength;

#if WindowsLinux==0
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif
	
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
			fread(buf, 1, readl, InFile);
			nRead += readl;
			readpos = 0;
		}

		//Voxel* LatticeElement = GetLattice(x,y,z);

		if(buf[readpos]){		//Solid

		//	SetLatticeBase(x, y, z, 0);
			GetLattice(x,y,z)->Solid |= VoxelSolid;

		}else{	//Non solid

		//	SetLatticeBase(x, y, z, nNonSolids+1);
			GetLattice(x,y,z)->Solid &= (~VoxelSolid);

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

		readpos++;
	}

ReadSolidsReturnRJMP:
	*NumberOfNonSolids = nNonSolids;
	delete[] buf;
	fclose(InFile);
	return;

}

void ReadSolids(unsigned int* NumberOfNonSolids){

	ReadSolids(NumberOfNonSolids, InputFilePath(SolidsFileName));

}

void ReadSolids(unsigned int* NumberOfNonSolids, char* FilePath){

	FILE* InFile = fopen(FilePath, "rb");
	
	if(InFile!=0){
		cout << "Reading in solids from '" << FilePath << "'" << endl;
	}else{
		cout << "Could not open solids file" << endl;
		return;
	}

	__int64 flength;

#if WindowsLinux==0
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif
	
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
			fread(buf, 1, readl, InFile);
			nRead += readl;
			readpos = 0;
		}

		if(buf[readpos]=='1'||buf[readpos]=='0'){

			//Voxel* LatticeElement = GetLattice(x,y,z);

			if(buf[readpos]=='1'){	//Solid

			//	SetLatticeBase(x, y, z, 0);
				GetLattice(x,y,z)->Solid |= VoxelSolid;

			}else{	//Non solid

			//	SetLatticeBase(x, y, z, nNonSolids+1);
				GetLattice(x,y,z)->Solid &= (~VoxelSolid);

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
	fclose(InFile);
	return;
}

void ReadVelocitiesBinSparse(){

	ReadVelocitiesBinSparse(InputFilePath(VelocitiesFileName));

}

void ReadVelocitiesBinSparse(char* FilePath){

	FILE* InFile = fopen(FilePath, "rb");

	if(InFile!=0){
		cout << "Reading in velocity field from '" << FilePath << "'" << endl;
	}else{
		cout << "Could not open velocity field file" << endl;
		return;
	}

	__int64 flength;

#if WindowsLinux==0
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif

	char initChar;
	fread(&initChar, 1, 1, InFile);

	rewind(InFile);

	const __int64 buffersize = 24*1048576;		//24MB

	char* buf = new char[buffersize];

	double v[3];					//Vector

	short x = 0;
	short y = 0;
	short z = 0;

	while(true){
		if(!IsElementSolid(x,y,z)){ break; }
		x++;
		if(x==NLattice_x){
			x=0;
			y++;
			if(y==NLattice_y){
				y=0;
				z++;
				if(z==NLattice_z){
					goto ReadVelocitiesReturnRJMP;
				}
			}
		}
	}
	
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
				goto ReadVelocitiesReturnRJMP;
			}
			fread(buf, 1, readl, InFile);
			nRead += readl;
			readpos = 0;
		}

		if((readl-readpos) < 24){
			cout << "File read failed..." << endl;
			goto ReadVelocitiesReturnRJMP;
		}

		Voxel* LatticeElement = GetLattice(x,y,z);

		memcpy(v, &buf[readpos], 24);

		LatticeElement->Vx = v[0];
		LatticeElement->Vy = v[1];
		LatticeElement->Vz = v[2];
		
		while(true){
			x++;
			if(x==NLattice_x){
				x=0;
				y++;
				if(y==NLattice_y){
					y=0;
					z++;
					if(z==NLattice_z){
						goto ReadVelocitiesReturnRJMP;
					}
				}
			}

			if(!IsElementSolid(x,y,z)){ break; }
		}

		readpos += 24;
	}

ReadVelocitiesReturnRJMP:
	delete[] buf;
	fclose(InFile);
	return;
}

void ReadVelocitiesBin(){

	ReadVelocitiesBin(InputFilePath(VelocitiesFileName));

}

void ReadVelocitiesBin(char* FilePath){

	FILE* InFile = fopen(FilePath, "rb");

	if(InFile!=0){
		cout << "Reading in velocity field from '" << FilePath << "'" << endl;
	}else{
		cout << "Could not open velocity field file" << endl;
		return;
	}

	__int64 flength;

#if WindowsLinux==0
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif

	char initChar;
	fread(&initChar, 1, 1, InFile);

	rewind(InFile);

	const __int64 buffersize = 24*1048576;		//24MB

	char* buf = new char[buffersize];

	double v[3];					//Vector

	short x = 0;
	short y = 0;
	short z = 0;
	
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
				goto ReadVelocitiesReturnRJMP;
			}
			fread(buf, 1, readl, InFile);
			nRead += readl;
			readpos = 0;
		}

		if((readl-readpos) < 24){
			cout << "File read failed..." << endl;
			goto ReadVelocitiesReturnRJMP;
		}

		Voxel* LatticeElement = GetLattice(x,y,z);

		memcpy(v, &buf[readpos], 24);

		LatticeElement->Vx = v[0];
		LatticeElement->Vy = v[1];
		LatticeElement->Vz = v[2];
		
		x++;
		if(x==NLattice_x){
			x=0;
			y++;
			if(y==NLattice_y){
				y=0;
				z++;
				if(z==NLattice_z){
					goto ReadVelocitiesReturnRJMP;
				}
			}
		}

		readpos += 24;
	}

ReadVelocitiesReturnRJMP:
	delete[] buf;
	fclose(InFile);
	return;
}

void ReadVelocitiesSparse(){

	ReadVelocitiesSparse(InputFilePath(VelocitiesFileName));

}

void ReadVelocitiesSparse(char* FilePath){

	FILE* InFile = fopen(FilePath, "rb");

	if(InFile!=0){
		cout << "Reading in velocity field from '" << FilePath << "'" << endl;
	}else{
		cout << "Could not open velocity field file" << endl;
		return;
	}

	__int64 flength;

#if WindowsLinux==0
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif

	char initChar;
	fread(&initChar, 1, 1, InFile);

	rewind(InFile);

	const __int64 buffersize = 1048576;		//1MB

	char* buf = new char[buffersize];

	int z=0;
	int y=0;
	int x=0;
	int component=0;

	while(true){
		if(!IsElementSolid(x,y,z)){ break; }
		x++;
		if(x==NLattice_x){
			x=0;
			y++;
			if(y==NLattice_y){
				y=0;
				z++;
				if(z==NLattice_z){
					goto ReadVelocitiesReturnRJMP;
				}
			}
		}
	}

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
				goto ReadVelocitiesReturnRJMP;
			}
			fread(buf, 1, readl, InFile);
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

			while(true){
				x++;
				if(x==NLattice_x){
					x=0;
					y++;
					if(y==NLattice_y){
						y=0;
						z++;
						if(z==NLattice_z){
							goto ReadVelocitiesReturnRJMP;
						}
					}
				}

				if(!IsElementSolid(x,y,z)){ break; }
			}

			}

		}else{		//Number 

			VBuffer[VBufl] = buf[readpos];
			VBufl++;

		}

		readpos++;
	}

ReadVelocitiesReturnRJMP:
	delete[] buf;
	fclose(InFile);
	return;
}

void ReadVelocities(){

	ReadVelocities(InputFilePath(VelocitiesFileName));

}

void ReadVelocities(char* FilePath){

	FILE* InFile = fopen(FilePath, "rb");

	if(InFile!=0){
		cout << "Reading in velocity field from '" << FilePath << "'" << endl;
	}else{
		cout << "Could not open velocity field file" << endl;
		return;
	}

	__int64 flength;

#if WindowsLinux==0
	_fseeki64(InFile, 0, SEEK_END);	//Seek end of file
	flength = _ftelli64(InFile);
	rewind(InFile);					//Seek beginning of file
#else
	fseeko64(InFile, 0, SEEK_END);	//Seek end of file
	flength = ftello64(InFile);
	rewind(InFile);					//Seek beginning of file
#endif

	char initChar;
	fread(&initChar, 1, 1, InFile);

	rewind(InFile);

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
				goto ReadVelocitiesReturnRJMP;
			}
			fread(buf, 1, readl, InFile);
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
						goto ReadVelocitiesReturnRJMP;
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

ReadVelocitiesReturnRJMP:
	delete[] buf;
	fclose(InFile);
	return;
}

void CreateWorkingLattice(){

	_CPUDomainWidthX = NLattice_x;
	_CPUDomainWidthY = NLattice_y;

	Lattice = new Voxel[((int)NLattice_x)*((int)NLattice_y)*((int)NLattice_z)];	//Array of non solid voxels + element 0 (solid)

	if(Lattice==0){
		cout << "Lattice memory allocation failed." << endl;
		return;
	}

	memset(Lattice, 0, ((int)NLattice_x)*((int)NLattice_y)*((int)NLattice_z)*sizeof(Voxel));

	unsigned int n;

	if(SolidsFileBin){
		ReadSolidsBin(&n, WorkingFilePath(WorkingGeometry));
	}else{
		ReadSolids(&n, WorkingFilePath(WorkingGeometry));
	}
	
	if(VelocitiesFileBin){
		if(VelocitiesFileSparse){
			ReadVelocitiesBinSparse(WorkingFilePath(WorkingVelocities));
		}else{
			ReadVelocitiesBin(WorkingFilePath(WorkingVelocities));
		}
	}else{
		if(VelocitiesFileSparse){
			ReadVelocitiesSparse(WorkingFilePath(WorkingVelocities));
		}else{
			ReadVelocities(WorkingFilePath(WorkingVelocities));
		}
	}

	cout << ConsoleSpacer << endl;

}

void CreateLattice(){		//Allocate lattice and read in solids and velocities

	_CPUDomainWidthX = NLattice_x;
	_CPUDomainWidthY = NLattice_y;

	Lattice = new Voxel[((int)NLattice_x)*((int)NLattice_y)*((int)NLattice_z)];	//Array of non solid voxels + element 0 (solid)

	if(Lattice==0){
		cout << "Lattice memory allocation failed." << endl;
		return;
	}

	//Zero Lattice

	int x=0;
	while(x!=NLattice_x){
	int y=0;
	while(y!=NLattice_y){
	int z=0;
	while(z!=NLattice_z){

		Voxel* LatticeElement = GetLattice(x,y,z);

		LatticeElement->Solid = 0;

		LatticeElement->Vx = 0.0;
		LatticeElement->Vy = 0.0;
		LatticeElement->Vz = 0.0;

		LatticeElement->DissolutionFraction = 0;

	z++;
	}
	y++;
	}
	x++;
	}

	unsigned int n;

	if(SolidsFileBin){
		ReadSolidsBin(&n);
	}else{
		ReadSolids(&n);
	}
	
	if(VelocitiesFileBin){
		if(VelocitiesFileSparse){
			ReadVelocitiesBinSparse();
		}else{
			ReadVelocitiesBin();
		}
	}else{
		if(VelocitiesFileSparse){
			ReadVelocitiesSparse();
		}else{
			ReadVelocities();
		}
	}

	cout << ConsoleSpacer << endl;
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

		bool ESolid = IsElementSolid(x,y,z);

		Voxel* LatticeElement = GetLattice(x,y,z);

		//Calculate velocities at the faces
		Voxel* LatticeElement_xn = GetLattice(xn,y,z);

		if(!ESolid){
			LatticeElement->u1 = 0.5*(LatticeElement->Vx + LatticeElement_xn->Vx);	//Velocity of negative x face
		}

			if(LatticeElement_xn->Solid&VoxelSolid  || (BoundaryConditionX==1 && x==0)){
				LatticeElement->u1=0;
				LatticeElement->Solid |= NegativeXSolid;
			}

		Voxel* LatticeElement_xp = GetLattice(xp,y,z);

		if(!ESolid){
			LatticeElement->u2 = 0.5*(LatticeElement->Vx + LatticeElement_xp->Vx);	//Velocity of positive x face
		}

			if(LatticeElement_xp->Solid&VoxelSolid || (BoundaryConditionX==1 && x==(NLattice_x-1))){
				LatticeElement->u2=0;
				LatticeElement->Solid |= PositiveXSolid;
			}


		Voxel* LatticeElement_yn = GetLattice(x,yn,z);

		if(!ESolid){
			LatticeElement->v1 = 0.5*(LatticeElement->Vy + LatticeElement_yn->Vy);	//Velocity of negative y face
		}

			if(LatticeElement_yn->Solid&VoxelSolid || (BoundaryConditionY==1 && y==0)){
				LatticeElement->v1=0;
				LatticeElement->Solid |= NegativeYSolid;
			}

		Voxel* LatticeElement_yp = GetLattice(x,yp,z);

		if(!ESolid){
			LatticeElement->v2 = 0.5*(LatticeElement->Vy + LatticeElement_yp->Vy);	//Velocity of positive y face
		}

			if(LatticeElement_yp->Solid&VoxelSolid || (BoundaryConditionY==1 && y==(NLattice_y-1))){
				LatticeElement->v2=0;
				LatticeElement->Solid |= PositiveYSolid;
			}

		Voxel* LatticeElement_zn = GetLattice(x,y,zn);

		if(!ESolid){
			LatticeElement->w1 = 0.5*(LatticeElement->Vz + LatticeElement_zn->Vz);	//Velocity of negative z face
		}

			if((LatticeElement_zn->Solid&VoxelSolid) || (BoundaryConditionZ==1 && z==0)){
				LatticeElement->w1=0;
				LatticeElement->Solid |= NegativeZSolid;
			}

		Voxel* LatticeElement_zp = GetLattice(x,y,zp);

		if(!ESolid){
			LatticeElement->w2 = 0.5*(LatticeElement->Vz + LatticeElement_zp->Vz);	//Velocity of positive z face
		}

			if(LatticeElement_zp->Solid&VoxelSolid || (BoundaryConditionZ==1 && z==(NLattice_z-1))){
				LatticeElement->w2=0;
				LatticeElement->Solid |= PositiveZSolid;
			}

	z++;
	}
	y++;
	}
	x++;
	}

}

bool OutputGeometryBin(){

	int Sz = ((int)NLattice_x)*((int)NLattice_y)*((int)NLattice_z);

	bool* Geometry = new bool[Sz];

	if(Geometry == 0){
		return false;
	}

	int i=0;
	while(i!=Sz){
		Geometry[i] = (Lattice[i].Solid&VoxelSolid);
		i++;
	}

	FILE* OutFile = fopen(WorkingFilePath(WorkingGeometry), "wb");

	fwrite(Geometry, 1, Sz, OutFile);		//Write state variables to file

	fclose(OutFile);

	delete[] Geometry;
	
}

void OutputGeometry(int n){

	ofstream outfile;
	outfile.open(OutputFilePath(GeometryOutputFile, (double)n),ios_base::trunc);
	
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Geometry Output" << endl;
	outfile << "ASCII" << endl;
	outfile << "DATASET STRUCTURED_POINTS" << endl;
	outfile << "DIMENSIONS " << NLattice_x << ' ' << NLattice_y << ' ' << NLattice_z << endl;
	outfile << "ORIGIN 0 0 0" << endl;
	outfile << "SPACING 1 1 1" << endl;
	outfile << "POINT_DATA " << ((int)NLattice_x)*((int)NLattice_y)*((int)NLattice_z) << endl;
	outfile << "SCALARS sample_scalars float" << endl;
	outfile << "LOOKUP_TABLE default" << endl;

	int z=0;
	while(z!=NLattice_z){
	int y=0;
	while(y!=NLattice_y){
	int x=0;
	while(x!=NLattice_x){

		if(IsElementSolid(x,y,z)){
			outfile << "1 ";
		}else{
			outfile << "0 ";
		}

	x++;
	}
	y++;
	}
	z++;
	}
	
	/*
	//Output particles
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Geometry Output" << endl;
	outfile << "ASCII" << endl;
	outfile << "DATASET POLYDATA" << endl;
	outfile << "POINTS " << NParticles << " float" << endl;

	int i=0;
	while(i!=NParticles){

		outfile << (float)Particles[i].X << " " << (float)Particles[i].Y << " " << (float)Particles[i].Z << " ";

		i++;
	}

	outfile << endl;
	outfile << "VERTICES " << NParticles << " " << 2*NParticles << endl;

	i=0;
	while(i!=NParticles){
		outfile << "1 " << i << " ";
		i++;
	}*/

	outfile.close();
	
}

void InitialiseParticles(){			//Initialise particles

	ParticleTypes = new ParticleType[2];

	//Non-reactive particle
	ParticleTypes[0].DissolutionFraction = 0;
	ParticleTypes[0].Mol = 0;

	//Reactive particle
	ParticleTypes[1].DissolutionFraction = 0.5;
	ParticleTypes[1].Mol = 0.1;
	
	if(ParticleNum == 0){	//Number of particles set per voxel

		NParticles = LatticeInfo.NNonSolids * ParticlesPerVoxel;
		
	}else{

		NParticles = ParticleNum;

	}

	Particles = new Particle[NParticles];

	if(Particles == 0){
		cout << "Memory allocation failure for particle array";
	}

	cout << "Number of particles: " << NParticles << endl;

	InitialiseParticlesUniform();
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

			Particles[i].Type = 1;

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

void WriteOutParameters(){
	//cout << sizeof(Particle) << ", " << sizeof(Voxel) << ", " << sizeof(unsigned short) << endl << endl;
	
	double AverageSpeed = sqrt((LatticeInfo.AverageFlowVelocity.X*LatticeInfo.AverageFlowVelocity.X) + (LatticeInfo.AverageFlowVelocity.Y*LatticeInfo.AverageFlowVelocity.Y) + (LatticeInfo.AverageFlowVelocity.Z*LatticeInfo.AverageFlowVelocity.Z));
	double Peclet = AverageSpeed / DiffusionCoefficient;

	cout << "Lattice dimensions (voxels) = (" << NLattice_x << ", " << NLattice_y << ", " << NLattice_z << ")" << endl;
	cout << "Lattice dimensions (micrometres) = (" << ((double)NLattice_x * LatticeResolutionMicron) << ", " << ((double)NLattice_y * LatticeResolutionMicron) << ", " << ((double)NLattice_z * LatticeResolutionMicron) << ")" << endl;
	cout << "Lattice resolution = " << LatticeResolutionMicron << " microns" << endl;
	cout << "Porosity = " << LatticeInfo.Porosity*100 << "%" << endl;
	cout << "Average velocity = (" << LatticeInfo.AverageFlowVelocity.X << ", " << LatticeInfo.AverageFlowVelocity.Y << ", " << LatticeInfo.AverageFlowVelocity.Z << ") Lu/s" << endl;
	cout << "Average velocity = (" << LatticeInfo.AverageFlowVelocity.X * LatticeResolutionMicron << ", " << LatticeInfo.AverageFlowVelocity.Y* LatticeResolutionMicron << ", " << LatticeInfo.AverageFlowVelocity.Z * LatticeResolutionMicron << ") microns/s" << endl;
	cout << "Diffusion Coefficient Dm = " << DiffusionCoefficient << endl;
	cout << "Peclet number = " << Peclet << " x Characteristic Length" << endl;

}

double GetTimestep(){			//Calculates ideal value of dt for simulation

	double dt;

	if(LatticeInfo.MaxFlowVelocityComponent >= Min_dV){
		double x = ( ( -sqrt(6*DiffusionCoefficient) + sqrt(6*DiffusionCoefficient + 4*LatticeInfo.MaxFlowVelocityComponent) )/(2*LatticeInfo.MaxFlowVelocityComponent) );
		dt = x*x;
	}else{
		dt = 1 / (24*DiffusionCoefficient);
	}

	return dt;
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

void CycleCoordinates(double* Px, double* Py, double* Pz){	//Applies boundary conditions to coordinates outside of lattice range

	double LSizeX = (double)NLattice_x;
	double LSizeY = (double)NLattice_y;
	double LSizeZ = (double)NLattice_z;

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

void RandomWalkVector(double L, Vector* V){		//Returns a random walk vector of length L

	double zp = 2.0*Random.GetUniform() - 1.0;

	double phi = acos(zp);
	double theta = 2*PI*Random.GetUniform();

	double LSinPhi = L*sin(phi);

	V->X = LSinPhi*cos(theta);
	V->Y = LSinPhi*sin(theta);
	V->Z = L*zp;
}

void RandomWalkParticle(int i, double dt, double RandomWalkLength, int ThreadID){

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

	bool VSolid;					//Is bounding voxel solid
	int SolidX, SolidY, SolidZ;		//Solid coordinate

	if(Ztime<=Ytime && Ztime<=Xtime){		//Reaches z boundary first

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

			SolidX = CurrentVoxelX;
			SolidY = CurrentVoxelY;
			SolidZ = C.zn;

		}else if(V.Z>0 && !VSolid){				//if moving in +z and neighbour is not solid

			P->Z = (double)C.zp;			//Move particle to voxel boundary and into upper voxel
			P->VoxelZ = C.zp;

		}else{									//if moving in +z and neighbour is solid

			P->Z = (double)(CurrentVoxelZ+1);			//Move particle to voxel boundary but inside current voxel

			SolidX = CurrentVoxelX;
			SolidY = CurrentVoxelY;
			SolidZ = C.zp;

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

			SolidX = CurrentVoxelX;
			SolidY = C.yn;
			SolidZ = CurrentVoxelZ;

		}else if(V.Y>0 && !VSolid){				//if moving in +y and neighbour is not solid

			P->Y = (double)C.yp;			//Move particle to voxel boundary and into right voxel
			P->VoxelY = C.yp;

		}else{									//if moving in +y and neighbour is solid

			P->Y = (double)(CurrentVoxelY+1);					//Move particle to voxel boundary but inside current voxel

			SolidX = CurrentVoxelX;
			SolidY = C.yp;
			SolidZ = CurrentVoxelZ;

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

			SolidX = C.xn;
			SolidY = CurrentVoxelY;
			SolidZ = CurrentVoxelZ;

		}else if(V.X>0 && !VSolid){				//if moving in +x and neighbour is not solid

			P->X = (double)C.xp;					//Move particle to voxel boundary and into right voxel
			P->VoxelX = C.xp;

		}else{									//if moving in +x and neighbour is solid

			P->X = (double)(CurrentVoxelX+1);		//Move particle to voxel boundary but inside current voxel

			SolidX = C.xp;
			SolidY = CurrentVoxelY;
			SolidZ = CurrentVoxelZ;

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

	if(VSolid && P->Type!=0){		//Reaction

		bool ParticleRemoved;

		ParticleSolidReaction(P, SolidX, SolidY, SolidZ, &ParticleRemoved, ThreadID);

		if(ParticleRemoved){
			return;
		}

	}

	}	//end of while loop

}

void ParticleSolidReaction(Particle* P, int SolidX, int SolidY, int SolidZ, bool* ParticleRemoved, int ThreadID){

	Voxel* V = GetLattice(SolidX,SolidY,SolidZ);

	//Fraction of solid particle dissolves
	double ParticleDissolutionFraction = ParticleTypes[P->Type].DissolutionFraction;

	double DisFractionBefore = V->DissolutionFraction;
	double DisFractionAfter = DisFractionBefore + ParticleDissolutionFraction;

	if(DisFractionAfter >= 1){		//Remove voxel

		V->Solid &= (~VoxelSolid);				//Make non-solid

		Coords C;
		GetCoordinates(SolidX,SolidY,SolidZ,&C);

		//Unset neighbour solid flags
		GetLattice(C.xn,C.y,C.z)->Solid &= (~PositiveXSolid);
		GetLattice(C.xp,C.y,C.z)->Solid &= (~NegativeXSolid);
		GetLattice(C.x,C.yn,C.z)->Solid &= (~PositiveYSolid);
		GetLattice(C.x,C.yp,C.z)->Solid &= (~NegativeYSolid);
		GetLattice(C.x,C.y,C.zn)->Solid &= (~PositiveZSolid);
		GetLattice(C.x,C.y,C.zp)->Solid &= (~NegativeZSolid);

		InterpolateVelocity(SolidX,SolidY,SolidZ, V);

		//Adjust lattice info
		LatticeInfo.NNonSolids++;
		LatticeInfo.NSolids--;
		LatticeInfo.Porosity = ((double)LatticeInfo.NNonSolids)/((double)(NLattice_x*NLattice_y*NLattice_z));

		//Probability of removing particle
		double PRemoved = (1 - DisFractionBefore) / ParticleDissolutionFraction;

		if(Random.GetUniform() <= PRemoved){	//Remove
			*ParticleRemoved = true;
		}else{
			*ParticleRemoved = false;
		}

		//Log alteration for output
		AlterationLog->LogAlteration((short)SolidX, (short)SolidY, (short)SolidZ, 0, ThreadID);

	}else{
		V->DissolutionFraction = DisFractionAfter;

		*ParticleRemoved = true;
	}

	if(!(*ParticleRemoved)){
		return;
	}

	//Remove particle and reassign randomly in medium

	while(true){	//Generate random position then check if non-solid

		double newX = (Random.GetUniform() * (double)(NLattice_x-1)) + 0.5;
		double newY = (Random.GetUniform() * (double)(NLattice_y-1)) + 0.5;
		double newZ = (Random.GetUniform() * (double)(NLattice_z-1)) + 0.5;

		int Vx = (int)newX;
		int Vy = (int)newY;
		int Vz = (int)newZ;

		if(!IsElementSolid(Vx,Vy,Vz)){

			P->VoxelX = Vx;
			P->VoxelY = Vy;
			P->VoxelZ = Vz;

			P->X = newX;
			P->Y = newY;
			P->Z = newZ;

			break;
		}

	}
	
	return;

}

void InterpolateVelocity(int SolidX, int SolidY, int SolidZ, Voxel* Vxl){

	Coords C;
	GetCoordinates(SolidX,SolidY,SolidZ,&C);

	Voxel* VxlXn = GetLattice(C.xn,C.y,C.z);
	Voxel* VxlXp = GetLattice(C.xp,C.y,C.z);
	Voxel* VxlYn = GetLattice(C.x,C.yn,C.z);
	Voxel* VxlYp = GetLattice(C.x,C.yp,C.z);
	Voxel* VxlZn = GetLattice(C.x,C.y,C.zn);
	Voxel* VxlZp = GetLattice(C.x,C.y,C.zp);

	double VXn = VxlXn->Vx;
	double VXp = VxlXp->Vx;
	double VYn = VxlYn->Vy;
	double VYp = VxlYp->Vy;
	double VZn = VxlZn->Vz;
	double VZp = VxlZp->Vz;

	if(VxlXn->Solid&VoxelSolid){ VXn = 0; }
	if(VxlXp->Solid&VoxelSolid){ VXp = 0; }
	if(VxlYn->Solid&VoxelSolid){ VYn = 0; }
	if(VxlYp->Solid&VoxelSolid){ VYp = 0; }
	if(VxlZn->Solid&VoxelSolid){ VZn = 0; }
	if(VxlZp->Solid&VoxelSolid){ VZp = 0; }

	//Simple average of neighbour vectors

	Vxl->Vx = 0.5 * (VXn + VXp);
	Vxl->Vy = 0.5 * (VYn + VYp);
	Vxl->Vz = 0.5 * (VZn + VZp);

	//Calculate new face velocities

	if(VxlXn->Solid&VoxelSolid){
		Vxl->u1 = 0;
	}else{
		Vxl->u1 = 0.5 * (Vxl->Vx + VxlXn->Vx);
		VxlXn->u2 = Vxl->u1;
	}
	if(VxlXp->Solid&VoxelSolid){ 
		Vxl->u2 = 0; 
	}else{
		Vxl->u2 = 0.5 * (Vxl->Vx + VxlXp->Vx);
		VxlXp->u1 = Vxl->u2;
	}
	if(VxlYn->Solid&VoxelSolid){ 
		Vxl->v1 = 0; 
	}else{
		Vxl->v1 = 0.5 * (Vxl->Vy + VxlYn->Vy);
		VxlYn->v2 = Vxl->v1;
	}
	if(VxlYp->Solid&VoxelSolid){ 
		Vxl->v2 = 0; 
	}else{
		Vxl->v2 = 0.5 * (Vxl->Vy + VxlYp->Vy);
		VxlYp->v1 = Vxl->v2;
	}
	if(VxlZn->Solid&VoxelSolid){ 
		Vxl->w1 = 0; 
	}else{
		Vxl->w1 = 0.5 * (Vxl->Vz + VxlZn->Vz);
		VxlZn->w2 = Vxl->w1;
	}
	if(VxlZp->Solid&VoxelSolid){ 
		Vxl->w2 = 0;
	}else{
		Vxl->w2 = 0.5 * (Vxl->Vz + VxlZp->Vz);
		VxlZp->w1 = Vxl->w2;
	}

}

void RunTimeStep(double dt, int ParticleIndex0, int ParticleIndex1, int ThreadID){

	double RandomWalkLength = sqrt(6*DiffusionCoefficient*dt);	//Computationally intensive, precalculate where possible

	int i = ParticleIndex0;
	while(i <= ParticleIndex1){

		Particle* P = &(Particles[i]);

		AdvectParticle(P, dt, ThreadID);
		RandomWalkParticle(i, dt, RandomWalkLength, ThreadID);	//Moves particle along a random path

	i++;
	}

}

void RunTimeStep(double dt, int ThreadID){

	double RandomWalkLength = sqrt(6*DiffusionCoefficient*dt);	//Computationally intensive, precalculate where possible

	int i=0;
	while(i!=NParticles){

		Particle* P = &(Particles[i]);

		AdvectParticle(P, dt, ThreadID);
		RandomWalkParticle(i, dt, RandomWalkLength, ThreadID);	//Moves particle along a random path

	i++;
	}

}

struct ThreadStruct{
	int ID;						//Thread index
	int ParticlesIndex0;		//Index of first particle under thread's auspices
	int ParticlesIndex1;		//Index of last particle under thread's auspices
	int Status;					//Running or waiting

	double t;					//Simulation time
	double tmax;				//Max simulation time
	double dt;					//Step time
	int TimeSteps;				//Number of time steps dt carried out
	int LogInterval;			//Number of timesteps between recalculating dispersion coefficient or propagator graph
	double RelativeBinWidth;	//Proportion of average flow distance for z bins in propagator distribution
	int GeometryOutputIndex;
};

void SyncThreads(ThreadStruct* Threads, int ThreadNum){

	while(true){
		bool Proceed = true;
		int i=1;
		while(i!=ThreadNum){
			if(Threads[i].Status != 1){ Proceed = false; }
			i++;
		}
		if(Proceed){ break; }
		Sleep(1);
	}

	return;
}

#if WindowsLinux==0				//Windows multi-threading

unsigned int __stdcall ThreadMain(void* data){

#else							//Linux multi-threading
void* ThreadMain(void* data){
#endif

	ThreadStruct* Thread = (ThreadStruct*)data;

	int ThreadID = Thread->ID;

	ThreadStruct* Threads;					//Array of all threads
	if(ThreadID == 0){
		Threads = (ThreadStruct*)data;
	}

	int ThreadNum = MultiThreadingSharedMem_TNum;
	int NParticlesThread = Thread->ParticlesIndex1 - Thread->ParticlesIndex0 + 1;

	//printf("Greetings from thread %i, dealing with %i particles from index %i\n" , Thread->ID, NParticlesThread, Thread->ParticlesIndex0);

	double t = Thread->t;									//Simulation time
	double tmax = Thread->tmax;								//Max simulation time
	double dt = Thread->dt;									//Step time
	int TimeSteps = Thread->TimeSteps;						//Number of time steps dt carried out
	int LogInterval = Thread->LogInterval;					//Number of timesteps between recalculating dispersion coefficient or propagator graph
	double RelativeBinWidth = Thread->RelativeBinWidth;		//Proportion of average flow distance for z bins in propagator distribution

	int ParticlesIndex0 = Thread->ParticlesIndex0;
	int ParticlesIndex1 = Thread->ParticlesIndex1;

	int GeometryOutputIndex = Thread->GeometryOutputIndex;

	ofstream outfile;

	if(ThreadID == 0){

		outfile.open(OutputFilePath(DataOutputFile));
		outfile << "Index\tTime\tPorosity\tNo. Solids" << endl;
		outfile << "0\t" << t << "\t" << LatticeInfo.Porosity << '\t' << LatticeInfo.NSolids << endl;
		outfile.close();

		OutputGeometry(0);

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

		RunTimeStep(dt, ParticlesIndex0, ParticlesIndex1, ThreadID);		//Simulate thread's particles

		if(ThreadID==0){
			cout << "Porosity = " << LatticeInfo.Porosity*100 << "%" << endl;
			cout << "NNonSolids = " << LatticeInfo.NNonSolids << endl;
		}

		TimeSteps++;
		t+=dt;

		if(!SimulationContinue){

			Thread->Status = 1;

			if(ThreadID==0){
				SyncThreads(Threads, ThreadNum);	//Wait for other threads to complete

				AlterationLog->OutputAlterations();

				SimulationStateStruct State;		//Record simulation state
				State.GeometryOutputIndex = GeometryOutputIndex;
				State.LogInterval = LogInterval;
				State.NParticles = NParticles;
				State.SimulationTime = t;
				State.SimulationTimeMax = tmax;
				State.SimulationTimeStep = dt;
				State.ThreadNum = ThreadNum;
				State.TimeSteps = TimeSteps;

				OutputSimulationState(&State);

				cout << "Run velocity update..." << endl;
			}
			return 0;
		}

		if(TimeSteps == LogInterval){

			double dt_remainder = (double)OutputTimeInterval - (double)TimeSteps*dt;		//Carry out any remaining time
			if(dt_remainder > (10E-3)){
				RunTimeStep(dt_remainder, ParticlesIndex0, ParticlesIndex1, ThreadID);
				t += dt_remainder;
			}

			if(ThreadID == 0){

				SyncThreads(Threads, ThreadNum);	//Wait for other threads to complete

				//Write outputs

				cout << "Output to file... " << endl;
				OutputGeometry(GeometryOutputIndex);

				outfile.open(OutputFilePath(DataOutputFile),ios::app);
				outfile << GeometryOutputIndex << "\t" << t << "\t" << LatticeInfo.Porosity << '\t' << LatticeInfo.NSolids << endl;
				outfile.close();

				GeometryOutputIndex++;

				//

				int i=1;
				while(i!=ThreadNum){
					Threads[i].Status = 0;		//Allow threads to continue
					i++;
				}

				if(LatticeInfo.NNonSolids >= (NLattice_x*NLattice_y*NLattice_z)){
					cout << "Complete dissolution." << endl;
					return 0;
				}	

			}else{

				Thread->t = t;
				Thread->TimeSteps = TimeSteps;

				Thread->Status = 1;

				while(true){		//Wait for thread 0 to output before continuing
					if(Thread->Status == 0){ break; }
					Sleep(1);
				}

				if(!SimulationContinue){
					return 0;
				}

				if(LatticeInfo.NNonSolids >= (NLattice_x*NLattice_y*NLattice_z)){
					return 0;
				}

			}
			 
			TimeSteps = 0;

		}

	}

	return 0;
}

int MainThreading(int argc, char *argv[]){

	if(MultiThreadingSharedMem_TNum==1){
		cout << "Running mode: Single thread mode" << endl;
	}else{
		cout << "Running mode: Multithread, shared memory with " << MultiThreadingSharedMem_TNum << " threads" << endl;
	}

	if(DebugAdvectionMode>=2){
		cout << "[Running in debug mode " << DebugAdvectionMode << "]" << endl;
		AdvectDebug.open(OutputFilePath(AdvectionDebugFile,0),ios::trunc);
	}

	//Test to see if during a simulation

	SimulationStateStruct SimulationState;
	bool ContinuingSimulation = false;

	ifstream WorkingGeometryOpen(WorkingFilePath(WorkingGeometry));
	if(!WorkingGeometryOpen){

		//Start from the beginning

		//Initialise Lattice	////////////////////////////////////

		CreateLattice();			//Read in from goemetry/velocity files

		InitialiseLattice();		//Calculate border velocities

		SetLatticeInfo();			//Fills in global LatticeInfo struct

		Random.TimeSeed();			//Seed random number generator based on local time

		InitialiseParticles();		//Initialise particles positions

		WriteOutParameters();		//Write out diffusion coefficient and average velocities to console

		OutputGeometry(0);			//Output geometry

	}else{
		WorkingGeometryOpen.close();
		ContinuingSimulation = true;
		cout << "Continuing simulation" << endl;

		//Continue from mid-simulation

		ReadSimulationState(&SimulationState);	//Read in simulation variables and particles

		Random.TimeSeed();			//Seed random number generator based on local time

		WriteOutParameters();		//Write out diffusion coefficient and average velocities to console

	}


	//Simulation	////////////////////////////////////////////

	double t = 0;						//Simulation time
	double dt = SimulationTimestep;		//Timestep
	double tmax = SimulationTimeMax;	//Simulation end time
	int TimeSteps = 0;					//Number of time steps dt carried out between log intervals
	int GeometryOutputIndex = 1;		//Index of geometry output files
	int LogInterval;					//Number of timesteps between outputting geometry

	if(!ContinuingSimulation){			//Starting from beginning of simulation

		if(tmax == 0){
			tmax = DBL_MAX;					//No limit
		}

		if(SimulationTimestep == 0){
			dt = GetTimestep();				//Optimum timestep
			dt = OutputTimeInterval / ceil(OutputTimeInterval/dt);	//Fit evenly into OutputTimeInterval
		}

		LogInterval = (int)floor(OutputTimeInterval / dt);	//Number of timesteps between outputting geometry

	}else{

		//Read in from simulation data file

		t = SimulationState.SimulationTime;
		dt = SimulationState.SimulationTimeStep;
		tmax = SimulationState.SimulationTimeMax;
		TimeSteps = SimulationState.TimeSteps;
		GeometryOutputIndex = SimulationState.GeometryOutputIndex;
		LogInterval = SimulationState.LogInterval;

	}

	cout << "Timestep dt = " << dt << endl;
	cout << "Porosity = " << LatticeInfo.Porosity*100 << "%" << endl;
	cout << "NNonSolids = " << LatticeInfo.NNonSolids << endl;

	cout << ConsoleSpacer << endl;

	int ThreadNum = MultiThreadingSharedMem_TNum;

	ThreadStruct* Threads = new ThreadStruct[ThreadNum];

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

	AlterationLog = new GeometryAlterationLog(ThreadNum);

	Threads[0].ID = 0;
	Threads[0].ParticlesIndex0 = 0;
	Threads[0].ParticlesIndex1 = ThreadParticleCount[0] - 1;

	Threads[0].t = t;
	Threads[0].dt = dt;
	Threads[0].tmax = tmax;
	Threads[0].LogInterval = LogInterval;
	Threads[0].TimeSteps = TimeSteps;
	Threads[0].GeometryOutputIndex = GeometryOutputIndex;

	i=1;
	while(i!=ThreadNum){
		
		Threads[i].ID = i;
		Threads[i].ParticlesIndex0 = Threads[i-1].ParticlesIndex1 + 1;
		Threads[i].ParticlesIndex1 = Threads[i].ParticlesIndex0 + ThreadParticleCount[i] - 1;

		Threads[i].t = t;
		Threads[i].dt = dt;
		Threads[i].tmax = tmax;
		Threads[i].LogInterval = LogInterval;
		Threads[i].TimeSteps = TimeSteps;
		Threads[i].GeometryOutputIndex = GeometryOutputIndex;

		Threads[i].Status = 1;	//Paused state
		
#if WindowsLinux==0				//Windows multi-threading
		_beginthreadex(0, 0, ThreadMain, &Threads[i], 0, 0);	//Threads
#else
		pthread_t tid;
		pthread_create(&tid, 0, &ThreadMain, (void*)(&Threads[i]));
#endif

		i++;
	}

	ThreadMain((void*)(&Threads[0]));	//Thread 0

	delete AlterationLog;
	delete[] ThreadParticleCount;
	delete[] Threads;

	delete[] Lattice;
	delete[] Particles;

	return 0;
}



int main(int argc, char *argv[]){


	MainThreading(argc,argv);

	
}

void AdvectParticle(Particle* P, double dt, int ThreadID){

#if DebugAdvectionMode>=2
	DebugParticle = false;
	double Restoredt;
	Particle RestoreParticle;
	
	//Store original particle state in case need to debug later
	Restoredt = dt;
	RestoreParticle.VoxelX = P->VoxelX;
	RestoreParticle.VoxelY = P->VoxelY;
	RestoreParticle.VoxelZ = P->VoxelZ;
	RestoreParticle.X = P->X;
	RestoreParticle.Y = P->Y;
	RestoreParticle.Z = P->Z;

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

	if(LatticeElement->Vx < AdvectionThreshold && LatticeElement->Vy < AdvectionThreshold && LatticeElement->Vz < AdvectionThreshold){
		return;
	}

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
	if(DebugParticle && OutputDebugData){
		AdvectDebug << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ") with " << NSolids << " solids" << endl;
		AdvectDebug << "Voxel Information:" << endl << "(Vx, Vy, Vz) = " << "(" << LatticeElement->Vx << ", " << LatticeElement->Vy << ", " << LatticeElement->Vz << ")" << endl;
		AdvectDebug << "u1 = " << LatticeElement->u1 << endl << "u2 = " << LatticeElement->u2 << endl << "v1 = " << LatticeElement->v1 << endl << "v2 = " << LatticeElement->v2 << endl << "w1 = " << LatticeElement->w1 << endl << "w2 = " << LatticeElement->w2 << endl;
	}
#endif

	switch(NSolids){

		case 0:
			CaseNoSolids(P, &C, LatticeElement, &dt);
			break;

		case 1:
			CaseOneSolid(P, &C, LatticeElement, &dt);
			break;

		case 2:
			CaseTwoSolids(P, &C, LatticeElement, &dt);
			break;
			
		case 3:
			CaseThreeSolids(P, &C, LatticeElement, &dt);
			break;

		case 4:
			CaseFourSolids(P, &C, LatticeElement, &dt);
			break;

		case 5:
		case 6:
			dt = 0;
			break;

		}

#if DebugAdvectionMode==1		//Just catch infinite loops

	if(cyclecount==8){
		cout << "Particle with " << NSolids << " solids exceeded 8 cycles, terminating." << endl;
		cout << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;
		AdvectDebug << "Particle with " << NSolids << " solids exceeded 8 cycles, terminating." << endl;
		AdvectDebug << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;

		dt = 0;
		break;
	}

#elif DebugAdvectionMode>=2		//Catch infinite loops, 1.#IND errors and particles in solid voxels
	
	if(cyclecount == 8){
		if(!OutputDebugData){
			cout << "Particle with " << NSolids << " solids exceeded 8 cycles." << endl;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;
		}else if(!DebugParticle){
			cout << "Trouble with particle with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Trouble with particle with " << NSolids << " solids. Debugging:" << endl;

			dt = Restoredt;							//Restore and debug
			P->VoxelX = RestoreParticle.VoxelX;
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;
			DebugParticle = true;
		
			goto AdvectDebugRestart;
		}else{
			AdvectDebug << "Exceeded 8 cycles. Terminating." << endl << endl;
			DebugParticle = false;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;
		}
	}

	//Inf and Nans

	if(isInfNaN(P->X) || isInfNaN(P->Y) || isInfNaN(P->Z)){
		if(!OutputDebugData){
			cout << "Indef. fault after advection: Particle with " << NSolids << " solids." << endl;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;

		}else if(!DebugParticle){
			cout << "Indef. fault after advection: Particle with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Indef. fault after advection: Particle with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Error in position: (" << P->X << ", " << P->Y << ", " << P->Z << ")" << endl;

			dt = Restoredt;
			P->VoxelX = RestoreParticle.VoxelX;
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;
			DebugParticle = true;
		
			goto AdvectDebugRestart;

		}else{
			AdvectDebug << endl << ConsoleSpacer << endl;
			DebugParticle = false;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;
		}

	}

	if(IsElementSolid(P->VoxelX, P->VoxelY, P->VoxelZ)){
		if(!OutputDebugData){
			cout << "Particle in solid after advection with " << NSolids << " solids." << endl;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;

		}else if(!DebugParticle){
			cout << "Particle in solid after advection with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Particle in solid after advection with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Particle found in solid at: (" << P->X << ", " << P->Y << ", " << P->Z << ")" << endl;
			AdvectDebug << "Voxel (solid): (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;

			dt = Restoredt;
			P->VoxelX = RestoreParticle.VoxelX;
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;
			DebugParticle = true;
		
			goto AdvectDebugRestart;
		}else{
			AdvectDebug << endl << ConsoleSpacer << endl;
			DebugParticle = false;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;
		}
	}

	if(P->X < (P->VoxelX-1) || P->X > (P->VoxelX+2) || P->Y < (P->VoxelY-1) || P->Y > (P->VoxelY+2) || P->Z < (P->VoxelZ-1) || P->Z > (P->VoxelZ+2)){
		if(!OutputDebugData){
			cout << "Particle position invalid after " << NSolids << " solids." << endl;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;

		}else if(!DebugParticle){
			cout << "Particle position invalid after " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Particle position invalid after " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Position: (" << P->X << ", " << P->Y << ", " << P->Z << "), Voxel: (" << P->VoxelX << ", " << P->VoxelY << ", " << P->VoxelZ << ")" << endl;
			
			dt = Restoredt;
			P->VoxelX = RestoreParticle.VoxelX;
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;
			DebugParticle = true;
		
			goto AdvectDebugRestart;
		}else{
			AdvectDebug << endl << ConsoleSpacer << endl;
			DebugParticle = false;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;
		}
	}

#endif
#if DebugAdvectionMode >= 3		//Catch everything including indef taus

	if(TauError){		//Indef taus
		
		TauError = false;

		if(!DebugParticle){
			cout << "Tau Indef on particle with " << NSolids << " solids. Debugging:" << endl;
			AdvectDebug << "Tau Indef on particle with " << NSolids << " solids. Debugging:" << endl;

			dt = Restoredt;							//Restore and debug
			P->VoxelX = RestoreParticle.VoxelX;
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;
			DebugParticle = true;
		
			goto AdvectDebugRestart;
		}else{
			AdvectDebug << endl << ConsoleSpacer << endl;
			DebugParticle = false;

			P->VoxelX = RestoreParticle.VoxelX;		//Restore for randomwalk
			P->VoxelY = RestoreParticle.VoxelY;
			P->VoxelZ = RestoreParticle.VoxelZ;
			P->X = RestoreParticle.X;
			P->Y = RestoreParticle.Y;
			P->Z = RestoreParticle.Z;	
			dt = 0;

			break;
		}
	}

#endif

	}

#if DebugAdvection
	DebugParticle = false;
#endif
	
}

void CaseNoSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt){

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

	double Px = P->X;					//X coordinate of particle
	double Py = P->Y;					//Y coordinate of particle
	double Pz = P->Z;					//Z coordinate of particle

	//cout << "Position: (" << Px << ", " << Py << ", " << Pz << "), Voxel: (" << x << ", " << y << ", " << z << ")" << endl;

	double Vx = (Px-x1)*(u2-u1) + u1;		//Velocity in x = ((Px-x1)/dx)*(u2-u1) + u1;
	double Vy = (Py-y1)*(v2-v1) + v1;		//Velocity in y = ((Py-y1)/dy)*(v2-v1) + v1;
	double Vz = (Pz-z1)*(w2-w1) + w1;		//Velocity in z = ((Pz-z1)/dz)*(w2-w1) + w1;

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
	if(Vx==0 || (u2<0 && u1>0) || (Vx>0 && u2<=0) || (Vx<0 && u1>=0)){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 1" << endl; }
			tau_x = DBL_MAX;
	}else if(fabs(u2-u1)>Min_dV){		
		if(Vx>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 2" << endl; }
			tau_x = log(u2/(u1 + (u2-u1)*(Px-x1)))/(u2-u1);		// = (dx/(u2-u1))*log((u2*dx)/(u1*dx+(u2-u1)*(Px-x1)));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 3" << endl; }
			tau_x = log(u1/(u2 + (u1-u2)*(x2-Px)))/(u2-u1);		// = (dx/(u2-u1))*log((u1*dx)/(u2*dx+(u1-u2)*(x2-Px)));
		}
	}else{
		if(Vx>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 4" << endl; }
			tau_x =	(x2-Px)/Vx;									// = 2*(x2-Px)/(u1+u2);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 5" << endl; }
			tau_x = (x1-Px)/Vx;									// = 2*(x1-Px)/(u1+u2);
		}
	}

	double tau_y;

	if(Vy==0 || (v2<0 && v1>0) || (Vy>0 && v2<=0) || (Vy<0 && v1>=0)){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 1" << endl; }
			tau_y = DBL_MAX;
	}else if(fabs(v2-v1)>Min_dV){
		if(Vy>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2" << endl; }
			tau_y = log(v2/(v1 + (v2-v1)*(Py-y1)))/(v2-v1);	// = (dy/(v2-v1))*log((v2*dy)/(v1*dy+(v2-v1)*(Py-y1)));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 3" << endl; }
			tau_y = log(v1/(v2 + (v1-v2)*(y2-Py)))/(v2-v1);	// = (dy/(v2-v1))*log((v1*dy)/(v2*dy+(v1-v2)*(y2-Py)));
		}
	}else{
		if(Vy>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 4" << endl; }
			tau_y = (y2-Py)/Vy;								// = 2*(y2-Py)/(v1+v2);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 5" << endl; }
			tau_y = (y1-Py)/Vy;								// = 2*(y1-Py)/(v1+v2);
		}
	}
		
	double tau_z;
	if(Vz==0 || (w2<0 && w1>0) || (Vz>0 && w2<=0) || (Vz<0 && w1>=0)){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
			tau_z = DBL_MAX;
	}else if(fabs(w2-w1)>Min_dV){
		if(Vz>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
			tau_z = log(w2/(w1 + (w2-w1)*(Pz-z1)))/(w2-w1);	// = (dz/(w2-w1))*log((w2*dz)/(w1*dz+(w2-w1)*(Pz-z1)));
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 3" << endl; }
			tau_z = log(w1/(w2 + (w1-w2)*(z2-Pz)))/(w2-w1);	// = dz/(w2-w1))*log((w1*dz)/(w2*dz+(w1-w2)*(z2-Pz)));
		}
	}else{
		if(Vz>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 4" << endl; }
			tau_z = (z2-Pz)/Vz;								// = 2*(z2-Pz)/(w1+w2);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 5" << endl; }
			tau_z = (z1-Pz)/Vz;								// = 2*(z1-Pz)/(w1+w2);
		}
	}

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << tau_x << ", " << tau_y << ", " << tau_z << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(tau_x)||isInfNaN(tau_y)||isInfNaN(tau_z)){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(tau_x)||tau_x<0){
		tau_x = DBL_MAX;
	}
	if(isInfNaN(tau_y)||tau_y<0){
		tau_y = DBL_MAX;
	}
	if(isInfNaN(tau_z)||tau_z<0){
		tau_z = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 1" << endl; }
			double u1_u2u1 = u1/(u2-u1);
			xe = x1 - u1_u2u1 + (u1_u2u1 + (Px-x1)) * exp((u2-u1)*time);	// = x1 - (u1*dx/(u2-u1)) + ((u1*dx/(u2-u1))+(Px-x1)) * exp((u2-u1)*time/dx);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 2" << endl; }
			xe = Px + time*Vx;
		}

		P->X = xe;
	}

	if(!ReachesBorderY){		//If reached boundary, new position has already been set

		if(Vy==0){
			ye = Py;
		}else if(fabs(v2-v1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 1" << endl; }
			double v1_v2v1 = v1/(v2-v1);
			ye = y1 - v1_v2v1 + (v1_v2v1 + (Py-y1)) * exp((v2-v1)*time);	// = y1 - (v1*dy/(v2-v1)) + ((v1*dy/(v2-v1))+(Py-y1)) * exp((v2-v1)*time/dy);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 2" << endl; }
			ye = Py + time*Vy;
		}

		P->Y = ye;
	}

	if(!ReachesBorderZ){		//If reached boundary, new position has already been set

		if(Vz==0){
			ze = Pz;
		}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = z1 - w1_w2w1 + (w1_w2w1 + (Pz-z1)) * exp((w2-w1)*time);	// = z1 - (w1*dz/(w2-w1)) + ((w1*dz/(w2-w1))+(Pz-z1)) * exp((w2-w1)*time/dz);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
			ze = Pz + time*Vz;
		}

		P->Z = ze;
	}

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Case 1 Solid:" << endl;
		AdvectDebug << "u1 = " << u1 << endl << "v1 = " << v1 << endl << "v2 = " << v2 << endl << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	if(u<0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 1" << endl; }
			Tx = (1 - (1/(1-Px)))/u1;					//Calculate tau_x = (dx*dx/u1)*((1/dx) - (1/(x2-Px)));
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 2" << endl; }
			Tx = DBL_MAX;
	}

	if(v==0 || (v2<=0 && v1>=0) || (v>0 && v2<=0) || (v<0 && v1>=0)){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 1" << endl; }
				Ty = DBL_MAX;
	}else if(fabs(v2-v1)>Min_dV){						//Calculate tau_y
		if(v>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2" << endl; }
				Ty = (1-Py)/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 3" << endl; }
				Ty = pow( v2/(v1 + (v2-v1)*Py) , 0.5*u1/(v2-v1) )/(u1*(1-Px)) - 1/(u1*(1-Px));	// = fabs((dx*dx/(u1*(x2-Px))) * pow( v2*dy/(v1*dy+(v2-v1)*(Py-y1)), u1*dy/(2*dx*(v2-v1)) ) - dx*dx/(u1*(x2-Px)));
				if(isInfNaN(Ty)){ Ty = DBL_MAX; }	//Simply suppress an inf in this case
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 4" << endl; }
				Ty = -Py/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 5" << endl; }
				Ty = pow( v1/(v1 + (v2-v1)*Py) , 0.5*u1/(v2-v1) )/(u1*(1-Px)) - 1/(u1*(1-Px));	// = fabs((dx*dx/(u1*(x2-Px))) * pow( v2*dy/(v1*dy+(v2-v1)*(y2-Py)), u1*dy/(2*dx*(v2-v1)) ) - dx*dx/(u1*(x2-Px)));
				if(isInfNaN(Ty)){ Ty = DBL_MAX; }	//Simply suppress an inf in this case
			}
		}
	}else{
		if(v>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 6" << endl; }
				Ty = (1-Py)/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 7" << endl; }
				Ty = (exp( (u1*(1-Py))/(2*v_face) ) - 1)/(u1*(1-Px));	//Modification: (v1+v2) -> v_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(y2-Py))/(dx*(v1+v2))) - 1);
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 8" << endl; }
				Ty = -Py/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 9" << endl; }
				Ty = (exp( -(u1*Py)/(2*v_face) ) - 1)/(u1*(1-Px));	//Modification: (v1+v2) -> v_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(y1-Py))/(dx*(v1+v2))) - 1);
			}
		}
	}

	if(w==0 || (w2<=0 && w1>=0) || (w>0 && w2<=0) || (w<0 && w1>=0)){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
				Tz = DBL_MAX;
	}else if(fabs(w2-w1)>Min_dV){								//Calculate tau_z
		if(w>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
				Tz = (1-Pz)/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 3" << endl << "w2/(w1 + (w2-w1)*Pz) = " << w2/(w1 + (w2-w1)*Pz) << endl << "0.5*u1/(w2-w1) = " << 0.5*u1/(w2-w1) << endl << "1/(u1*(1-Px)) = " << 1/(u1*(1-Px)) << endl; }
				Tz = pow( w2/(w1 + (w2-w1)*Pz) , 0.5*u1/(w2-w1) )/(u1*(1-Px)) - 1/(u1*(1-Px));		// = fabs((dx*dx/(u1*(x2-Px))) * pow( w2*dz/(w1*dz+(w2-w1)*(Pz-z1)), u1*dz/(2*dx*(w2-w1)) ) - dx*dx/(u1*(x2-Px)));
				if(isInfNaN(Tz)){ Tz = DBL_MAX; }	//Simply suppress an inf in this case
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 4" << endl; }
				Tz = -Pz/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 5" << endl << "w1/(w1 + (w2-w1)*Pz) = " << w1/(w1 + (w2-w1)*Pz) << endl << "0.5*u1/(w2-w1) = " << 0.5*u1/(w2-w1) << endl << "1/(u1*(1-Px)) = " << 1/(u1*(1-Px)) << endl; }
				Tz = pow( w1/(w1 + (w2-w1)*Pz) , 0.5*u1/(w2-w1) )/(u1*(1-Px)) - 1/(u1*(1-Px));	// = fabs((dx*dx/(u1*(x2-Px))) * pow( w2*dz/(w1*dz+(w2-w1)*(z2-Pz)), u1*dz/(2*dx*(w2-w1)) ) - dx*dx/(u1*(x2-Px)));
				if(isInfNaN(Tz)){ Tz = DBL_MAX; }	//Simply suppress an inf in this case
			}
		}
	}else{
		if(w>0){
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 6" << endl; }
				Tz = (1-Pz)/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 7" << endl; }
				Tz = (exp( (u1*(1-Pz))/(2*w_face) ) - 1)/(u1*(1-Px));	//Modification: (w1+w2) -> w_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(z2-Pz))/(dx*(w1+w2))) - 1);
			}
		}else{
			if(u1==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 8" << endl; }
				Tz = -Pz/w;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 9" << endl; }
				Tz = (exp( -(u1*Pz)/(2*w_face) ) - 1)/(u1*(1-Px));	//Modification: (w1+w2) -> w_face*2; = (dx*dx/(u1*(x2-Px))) * (exp((u1*(z1-Pz))/(dx*(w1+w2))) - 1);
			}
		}
	}

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(Tx)||Tx<0||isInfNaN(Ty)||Ty<0||isInfNaN(Tz)||Tz<0){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(Tx)||Tx<0){
		Tx = DBL_MAX;
	}
	if(isInfNaN(Ty)||Ty<0){
		Ty = DBL_MAX;
	}
	if(isInfNaN(Tz)||Tz<0){
		Tz = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 1" << endl; }
			xe = 0;
	}else if(u1==0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 2" << endl; }
			xe = Px;
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 3" << endl; }
			xe = 1 - 1/(t*u1 + 1/(1-Px));		// = x2 - 1/(t*u1/(dx*dx)+1/(x2-Px));
	}

	if(Case->ReachedNegativeBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 1" << endl; }
			ye = 0;
	}else if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 2" << endl; }
			ye = 1;
	}else if(fabs(v2-v1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 3" << endl; }
			double v1_v2v1 = v1/(v2-v1);
			ye = -v1_v2v1 + (v1_v2v1 + Py) * pow( 1+(u1*(1-Px)*t) , 2*(v2-v1)/u1 );		// = y1 - (v1*dy/(v2-v1)) + ((v1*dy+(v2-v1)*(Py-y1))/(v2-v1)) * pow( 1+(u1*t*(x2-Px)/(dx*dx)) , 2*dx*(v2-v1)/(u1*dy) );
	}else{
		if(u1==0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 4" << endl; }
			ye = Py + v*t;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 5" << endl; }
			ye = Py + log( 1 + (u1*(1-Px)*t) )*(v1+v2)/u1;		// = Py + (dx*(v1+v2)/u1) * log((u1*(x2-Px)*t)/(dx*dx) + 1);
		}
	}

	if(Case->ReachedNegativeCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
			ze = 0;
	}else if(Case->ReachedPositiveCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
			ze = 1;
	}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 3" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = -w1_w2w1 + (w1_w2w1 + Pz) * pow( 1+(u1*(1-Px)*t) , 2*(w2-w1)/u1 );	// = z1 - (w1*dz/(w2-w1)) + ((w1*dz+(w2-w1)*(Pz-z1))/(w2-w1)) * pow((1+(u1*t*(x2-Px)/(dx*dx))) , (2*dx*(w2-w1)/(u1*dz)));
	}else{
		if(u1==0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 4" << endl; }
			ze = Pz + w*t;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 5" << endl; }
			ze = Pz + log( 1 + (u1*(1-Px)*t) )*(w1+w2)/u1;		// = Pz + (dx*(w1+w2)/u1) * log((u1*(x2-Px)*t)/(dx*dx) + 1);
		}
	}

	/*

	if(xe<0||xe>1.00000001||ye<0||ye>1.00000001||ze<0||ze>1.00000001){
		cout << "1 Solid Fault - Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}
	*/

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

	Case->Pa_a1 = xe;
	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

}

//Advection function used when one neighbouring voxel is solid (6 cases)
void CaseOneSolid(Particle* P, Coords* C, Voxel* LatticeElement, double* dt){

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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Case 2 Opposite Solids:" << endl;
		AdvectDebug << "v1 = " << v1 << endl << "v2 = " << v2 << endl << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 1" << endl; }
			Tx = DBL_MAX;			//X will not reach border
	
	if(v==0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 1" << endl; }
			Ty = DBL_MAX;
	}else if(fabs(v2-v1)>Min_dV){	//Expression fails if denominator too small, treat as a special case if velocities close
		if(v>0 && v2<=0){			//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2" << endl; }
			Ty = DBL_MAX;
		}else if(v>0){				//Exit time to v2 face
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 3" << endl; }
			Ty = log( v2/(v1 + (v2-v1)*Py) )/(6*(v2-v1)*_1_PxPx);	// = (dx*dx*dy/(6*(v2-v1)*(x2-Px)*(Px-x1))) * log(v2*dy/(v1*dy+(v2-v1)*(Py-y1)));
		}else if(v<0 && v1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 4" << endl; }
			Ty = DBL_MAX;
		}else{						//Exit time to v1 face if v<0
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 5" << endl; }
			Ty = log( v1/(v1 + (v2-v1)*Py) )/(6*(v2-v1)*_1_PxPx);	// = (dx*dx*dy/(6*(v2-v1)*(x2-Px)*(Px-x1))) * -log(v2*dy/(v1*dy+(v2-v1)*(y2-Py)));
		}
	}else{
		if(v>0 && v2<=0){			//Take care of bizarre cases for rigour's sake - outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 6" << endl; }
			Ty = DBL_MAX;
		}else if(v>0){				//Exit time to v2 face
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 7" << endl; }
			Ty = (1-Py)/v;
		}else if(v<0 && v1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 8" << endl; }
			Ty = DBL_MAX;
		}else{						//Exit time to v1 face if v<0
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 9" << endl; }
			Ty = -Py/v;
		}
	}


	if(w==0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
			Tz = DBL_MAX;
	}else if(fabs(w2-w1)>Min_dV){
		if(w>0 && w2<=0){			//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
			Tz = DBL_MAX;
		}else if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 3" << endl; }
			Tz = log( w2/(w1 + (w2-w1)*Pz) )/(6*(w2-w1)*_1_PxPx);	// = (dx*dx*dz/(6*(w2-w1)*(x2-Px)*(Px-x1))) * log(w2*dz/(w1*dz+(w2-w1)*(Pz-z1)));
		}else if(w<0 && w1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 4" << endl; }
			Tz = DBL_MAX;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 5" << endl; }
			Tz = log( w1/(w1 + (w2-w1)*Pz) )/(6*(w2-w1)*_1_PxPx);	// = (dx*dx*dz/(6*(w2-w1)*(x2-Px)*(Px-x1))) * -log(w2*dz/(w1*dz+(w2-w1)*(z2-Pz)));
		}
	}else{
		if(w>0 && w2<=0){			//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 6" << endl; }
			Tz = DBL_MAX;
		}else if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 7" << endl; }
			Tz = (1-Pz)/w;
		}else if(w<0 && w1>=0){		//Outward direction but inward face velocity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 8" << endl; }
			Tz = DBL_MAX;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 9" << endl; }
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(Tx)||Tx<0||isInfNaN(Ty)||Ty<0||isInfNaN(Tz)||Tz<0){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(Tx)||Tx<0){
		Tx = DBL_MAX;
	}
	if(isInfNaN(Ty)||Ty<0){
		Ty = DBL_MAX;
	}
	if(isInfNaN(Tz)||Tz<0){
		Tz = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

	//cout << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 1" << endl; }

	double xe = Px;		//Positions after advection
	double ye;
	double ze;	

	if(Case->ReachedNegativeBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 1" << endl; }
		ye = 0;

	}else if(Case->ReachedPositiveBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 2" << endl; }
		ye = 1;

	}else if(fabs(v2-v1)>Min_dV){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 3" << endl; }
		double v1_v2v1 = v1/(v2-v1);
		ye = y1 - v1_v2v1 + (v1_v2v1 + Py) * exp( 6*(v2-v1)*_1_PxPx*t );	// = y1 - (v1*dy)/(v2-v1) + ( (v1*dy + (v2-v1)*(Py-y1))/((v2-v1)*dy) ) * exp( (6*(v2-v1)*(Px-x1)*(x2-Px)/(dx*dx*dy))*t );

	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 4" << endl; }	
		ye = Py + 3*(v1+v2)*_1_PxPx*t;		// = Py + ( 3*(v1+v2)*(x2-Px)*(Px-x1)/(dx*dx) )*t;

	}

	if(Case->ReachedNegativeCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
		ze = 0;

	}else if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
		ze = 1;

	}else if(fabs(w2-w1)>Min_dV){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 3" << endl; }
		double w1_w2w1 = w1/(w2-w1);
		ze = z1 - w1_w2w1 + (w1_w2w1 + Pz) * exp( 6*(w2-w1)*_1_PxPx*t );	// = z1 - (w1*dz)/(w2-w1) + ( (w1*dz + (w2-w1)*(Pz-z1))/((w2-w1)*dz) ) * exp( (6*(w2-w1)*(Px-x1)*(x2-Px)/(dx*dx*dz))*t );

	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 4" << endl; }
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Case 2 Adjacent Solids:" << endl;
		AdvectDebug << "u2 = " << u2 << endl << "v2 = " << v2 << endl << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	//cout << "Velocity: (" << u << ", " << v << ", " << w << ")" <<endl;

	double D = (u2 + v2)*_2PxPy;					// = 2*(u2*dy + v2*dx)*(Px-x1)*(Py-y1)/(dx*dx*dy*dy);
	double invD = 1/D;

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "D = " << D << endl << "invD = " << invD << endl;
	}

	double Tx;
	double Ty;
	double Tz;

	//Calculate time to reach the border, Tau
	
	//Tau x
	
	if(u > 0){
		if(u2 < Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 1" << endl; }
			Tx = (1-Px)/u;
		}else if(fabs(u2 + v2)>Min_dV){									// fabs(u2*dy + v2*dx)>Min_dV
			double Power = (u2 + v2)/u2;
		if(Power > 50){						//Power goes to 0
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 2a" << endl << "(u2 + v2)/u2 = " << (u2 + v2)/u2 << endl; }
			Tx = invD;
		}else if(Power < -50){				//Power goes to infinity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 2b" << endl << "(u2 + v2)/u2 = " << (u2 + v2)/u2 << endl; }
			Tx = DBL_MAX;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 2c" << endl << "(u2 + v2)/u2 = " << (u2 + v2)/u2 << endl; }
			Tx = invD * (1 - pow( Px , Power ));				// = (1/D) * (1 - pow( dx/(Px-x1) , -(u2*dy + v2*dx)/(u2*dy) ));
		}
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 3" << endl; }
			Tx = -log( Px )/(u2*_2PxPy);							// = ((dx*dx*dy*dy)/(2*u2*(Px-x1)*(Py-y1))) * log(dx/(Px-x1));
		}
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 4" << endl; }
			Tx = DBL_MAX;
	}

	//Tau y

	if(v > 0){
		if(v2 < Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 1" << endl; }
			Ty = (1-Py)/v;
		}else if(fabs(u2 + v2)>Min_dV){									// fabs(u2*dy + v2*dx)>Min_dV
			double Power = (u2 + v2)/v2;
		if(Power > 50){						//Power goes to 0
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2a" << endl << "(u2 + v2)/v2 = " << (u2 + v2)/v2 << endl; }
			Ty = invD;
		}else if(Power < -50){				//Power goes to infinity
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2b" << endl << "(u2 + v2)/v2 = " << (u2 + v2)/v2 << endl; }
			Ty = DBL_MAX;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2c" << endl << "(u2 + v2)/v2 = " << (u2 + v2)/v2 << endl; }
			Ty = invD * (1 - pow( Py , Power ));
		}
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 3" << endl; }
			Ty = -log( Py )/(v2*_2PxPy);							// = ((dx*dx*dy*dy)/(2*v2*(Px-x1)*(Py-y1))) * log(dy/(Py-y1));
		}
	}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 4" << endl; }
			Ty = DBL_MAX;
	}

	//Tau z
	if(w==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
				Tz = DBL_MAX;
	}else if(w>0){
		if(w2<=0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
				Tz = DBL_MAX;
		}else if(fabs(w2-w1)>Min_dV){
			if(fabs(u2 + v2)>Min_dV){														// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 3" << endl; }
				Tz = invD * (1 - pow( w2/(Pz*(w2-w1) + w1) , -(u2+v2)/(2*(w2-w1)) ));		// = (1/D)*(1 - pow( (dx*(w2-w1) + w1*dz)/((Pz-z1)*(w2-w1) + w1*dz) , -(dz*(u2*dy+v2*dx)/(2*(w2-w1)*dx*dy)) ));
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 4" << endl; }
				Tz = log( w2/(Pz*(w2-w1) + w1) )/(2*(w2-w1)*_2PxPy);						// = ((dx*dy*dz)/(4*(w2-w1)*(Px-x1)*(Py-y1))) * log((dx*(w2-w1)+(w1*dz))/((Pz-z1)*(w2-w1)+w1*dz));
			}
		}else{
			if(fabs(u2 + v2)>Min_dV){													// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 5" << endl; }
				Tz = invD * (1 - exp( -(1-Pz)*(u2+v2)/(2*w2) ));							// = (1/D) * (1 - exp( -(z2-Pz)*(u2*dy+v2*dx)/(2*w1*dx*dy) ));
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 6" << endl; }
				Tz = (1-Pz)/(2*w2*_2PxPy);													// = (x2-Pz)*dx*dy/(4*w1*(Px-x1)*(Py-y1));
			}
		}
	}else{
		if(w1>=0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 7" << endl; }
				Tz = DBL_MAX;
		}else if(fabs(w2-w1)>Min_dV){
			if(fabs(u2 + v2)>Min_dV){														// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 8" << endl; }
				Tz = invD * (1 - pow( w1/(Pz*(w2-w1) + w1) , -(u2+v2)/(2*(w2-w1)) ));		// = (1/D)*(1 - pow( (w1*dz)/((Pz-z1)*(w2-w1) + w1*dz) , -(dz*(u2*dy+v2*dx)/(2*(w2-w1)*dx*dy)) ));
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 9" << endl; }
				Tz = log( w1/(Pz*(w2-w1) + w1) )/(2*(w2-w1)*_2PxPy);						// = ((dx*dy*dz)/(4*(w2-w1)*(Px-x1)*(Py-y1))) * log((w1*dz)/((Pz-z1)*(w2-w1)+w1*dz));
			}
		}else{
			if(fabs(u2 + v2)>Min_dV){														// fabs(u2*dy + v2*dx)>Min_dV
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 10" << endl; }
				Tz = invD * (1 - exp( Pz*(u2+v2)/(2*w1) ));									// = (1 - exp( -(z1-Pz)*(u2*dy+v2*dx)/(2*w1*dx*dy) ))/D;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 11" << endl; }
				Tz = -Pz/(2*w1*_2PxPy);														// = (z1-Pz)*dx*dy/(4*w1*(Px-x1)*(Py-y1));
			}
		}
	}

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(Tx)||Tx<0||isInfNaN(Ty)||Ty<0||isInfNaN(Tz)||Tz<0){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(Tx)||Tx<0){
		Tx = DBL_MAX;
	}
	if(isInfNaN(Ty)||Ty<0){
		Ty = DBL_MAX;
	}
	if(isInfNaN(Tz)||Tz<0){
		Tz = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 1" << endl; }
			xe = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 2" << endl; }
			xe = Px * pow( 1-D*t , -u2/(u2 + v2) );								// = x1 + (Px-x1) * pow( 1-D*t , -(u2*dy)/(u2*dy + v2*dx) );
		}

		if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 1" << endl; }
			ye = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 2" << endl; }
			ye = Py * pow( 1-D*t , -v2/(u2 + v2) );								// = y1 + (Py-y1) * pow( 1-D*t , -(v2*dy)/(u2*dy + v2*dx) );
		}

		if(Case->ReachedNegativeCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
			ze = 0;
		}else if(Case->ReachedPositiveCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
			ze = 1;
		}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 3" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = -w1_w2w1 + (w1_w2w1 + Pz) * pow( 1-D*t , -2*(w2-w1)/(u2 + v2) );	// = z1 - ((w1*dz)/(w2-w1)) + ((w1*dz)/(w2-w1) + (Pz-z1)) * pow( 1-D*t , -(2*(w2-w1)*dx*dy)/(dz*(u2*dy + v2*dx)) );
			
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 4" << endl; }
			ze = Pz - ((2*w1)/(u2+v2)) * log( 1 - D*t );								// = Pz - ((2*w1*dx*dy)/(u2*dy+v2*dx)) * log( 1 - D*t );

		}

	}else{

		if(Case->ReachedPositiveABorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 3" << endl; }
			xe = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 4" << endl; }
			xe = Px * exp( u2*_2PxPy*t );											// = x1 + (Px-x1) * exp( (2*u2*(Px-x1)*(Py-y1)/(dx*dx*dy))*t );
		}

		if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 3" << endl; }
			ye = 1;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 4" << endl; }
			ye = Py * exp( v2*_2PxPy*t );											// = y1 + (Py-y1) * exp( (2*v2*(Px-x1)*(Py-y1)/(dx*dy*dy))*t );
		}

		if(Case->ReachedNegativeCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 5" << endl; }
			ze = 0;
		}else if(Case->ReachedPositiveCBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 6" << endl; }
			ze = 1;
		}else if(fabs(w2-w1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 7" << endl; }
			double w1_w2w1 = w1/(w2-w1);
			ze = -w1_w2w1 + (w1_w2w1 + Pz) * exp( 2*(w2-w1)*_2PxPy*t );				// = z1 - (w1*dz/(w2-w1)) + ( w1*dz/(w2-w1) + (Pz-z1) ) * exp( (4*(w2-w1)*(Px-x1)*(Py-y1)/(dx*dy*dz))*t );

		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 8" << endl; }
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

}

//Advection function used when two neighbouring voxels are solid (15 cases)
void CaseTwoSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt){

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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 1" << endl; }
		Tx = log( 1-Px )/(_4PxPyPz*u1);		// = (dx*dx*dy*dz*log((x2-Px)/dx)/(4*(x2-Px)*(y2-Py)*(Pz-z1)*u1));
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 2" << endl; }
		Tx = DBL_MAX;
	}

	//Tau y

	if(v<0){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 1" << endl; }
		Ty = log( 1-Py )/(_4PxPyPz*v1);		// = (dy*dy*dx*dz*log((y2-Py)/dy)/(4*(x2-Px)*(y2-Py)*(Pz-z1)*v1));
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2" << endl; }
		Ty = DBL_MAX;
	}

	//Tau z

	if(w>0){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
		Tz = -log( Pz )/(_4PxPyPz*w2);		// = (dz*dz*dx*dy*log(dz/(Pz-z1))/(4*(x2-Px)*(y2-Py)*(Pz-z1)*w2));
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
		Tz = DBL_MAX;
	}

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(Tx)||Tx<0||isInfNaN(Ty)||Ty<0||isInfNaN(Tz)||Tz<0){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(Tx)||Tx<0){
		Tx = DBL_MAX;
	}
	if(isInfNaN(Ty)||Ty<0){
		Ty = DBL_MAX;
	}
	if(isInfNaN(Tz)||Tz<0){
		Tz = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 1" << endl; }
		xe = 0;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 2" << endl; }
		xe = 1 - (1-Px) * exp(-u1*_4PxPyPz*t);		// = x2-(x2-Px)*exp(-4*u1*(x2-Px)*(y2-Py)*(Pz-z1)*t/(dx*dx*dy*dz));
	}

	if(Case->ReachedNegativeBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 1" << endl; }
		ye = 0;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 2" << endl; }
		ye = 1 - (1-Py) * exp(-v1*_4PxPyPz*t);		// = y2-(y2-Py)*exp(-4*v1*(x2-Px)*(y2-Py)*(Pz-z1)*t/(dx*dy*dy*dz));
	}

	if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
		ze = 1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
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

	
	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tx 1" << endl; }
		Tx = DBL_MAX;

	//Tau y

	if(v==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 0" << endl; }
				Ty = DBL_MAX;
	}else if(v>0){
		if(v2<=0 || (v2<0 && v1>0)){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 1" << endl; }
				Ty = DBL_MAX;
		}else if(fabs(v2-v1)>Min_dV){
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2" << endl; }
				Ty = (1-Py)/v;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 3" << endl; }
				Ty = (1 - pow( v2/v_face , -w2/(2*(v2-v1)) ))/(w2*_6PxPxPz);	// = (dx*dx*dz*dz/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1))) * (1 - pow( (v1*dy+(v2-v1)*dy)/(v1*dy+(v2-v1)*(Py-y1)) , -w2*dy/(2*(v2-v1)*dz) ));
			}
		}else{
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 4" << endl; }
				Ty = (1-Py)/v ;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 5" << endl; }
				Ty = (1 - exp( -w2*(1-Py)/(2*v_face) ))/(w2*_6PxPxPz);						// = ((dx*dx*dz*dz)/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1))) * (1 - exp( -w2*(y2-Py)/(2*v1*dz) ));
			}
		}
	}else{
		if(v1>=0 || (v2<0 && v1>0)){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 6" << endl; }
				Ty = DBL_MAX;
		}else if(fabs(v2-v1)>Min_dV){
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 7" << endl; }
				Ty = -Py/v ;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 8" << endl; }
				Ty = (1 - pow( v1/v_face , -w2/(2*(v2-v1)) ))/(w2*_6PxPxPz);	// = (dx*dx*dz*dz/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1)))*(1-pow((v1*dy)/(v1*dy+(v2-v1)*(Py-y1)), -w2*dy/(2*(v2-v1)*dz)));
			}
		}else{
			if(w2==0){
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 9" << endl; }
				Ty = -Py/v ;
			}else{
				if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 10" << endl; }
				Ty = (1 - exp( w2*Py/(2*v_face) ))/(w2*_6PxPxPz);							// = ((dx*dx*dz*dz)/(6*w2*(x2-Px)*(Px-x1)*(Pz-z1)))*(1-exp(-w2*(y1-Py)/(2*v1*dz)));
			}
		}
	}

	//Tau z

	if(w>0){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
		Tz = ( 1/(Pz-z1) - 1 )/( 6*(1-Px)*Px*w2 );										// = (dx*dx*dz*dz/(6*(x2-Px)*(Px-x1)*w2))*(1/(Pz-z1)-1/dz);
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(Tx)||Tx<0||isInfNaN(Ty)||Ty<0||isInfNaN(Tz)||Tz<0){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(Ty)||Ty<0){
		Ty = DBL_MAX;
	}
	if(isInfNaN(Tz)||Tz<0){
		Tz = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;
	}

	//cout << "Analytical Tau:      (" << Tx << ", " << Ty << ", " << Tz << ")" << endl;

	//Positions

	double xe, ye, ze;

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "xe 1" << endl; }
	xe = Px;

	if(Case->ReachedNegativeBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 1" << endl; }
			ye = 0;
	}else if(Case->ReachedPositiveBBorder){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 2" << endl; }
			ye = 1;
	}else if(fabs(v2-v1)>Min_dV){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 3" << endl; }
			double v1_v2v1 = v1/(v2-v1);
			ye = -v1_v2v1 + (v1_v2v1 + Py) * pow( 1-w2*_6PxPxPz*t , -2*(v2-v1)/w2 );		// = y1 - (v1*dy/(v2-v1)) + pow( 1-6*w2*(x2-Px)*(Px-x1)*(Pz-z1)*t/(dx*dx*dz*dz) , -2*(v2-v1)*dz/(w2*dy))*(v1*dy/(v2-v1)+(Py-y1));
	}else{
		if(w2==0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 4" << endl; }
			ye = Py + v*t;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 5" << endl; }
			ye = Py - (2*v1/w2) * log( 1 - w2*_6PxPxPz*t );											// = Py - (2*v1*dz/w2) * log( 1 - (6*w2*(x2-Px)*(Px-x1)*(Pz-z1)*t/(dx*dx*dz*dz)) );
		}
	}

	if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
		ze = 1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Raw Position: (" << xe << ", " << ye << ", " << ze << ")" << endl;
	}

}

//Advection function used when three neighbouring voxels are solid (20 cases)
void CaseThreeSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt){

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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
	
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
			AdvectDebug << "Positive Y Solid, Negative Y Solid, Negative X Solid" << endl; 
		}

		//a=y, b=z, c=x
		oppcoords.Velocity_b1 = w1;	//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = w2;	//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = u2;	//Outward velocity opposite solid c direction
		
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
			AdvectDebug << "Positive Z Solid, Negative Z Solid, Negative X Solid" << endl; 
		}

		//a=-z, b=y, c=x
		oppcoords.Velocity_b1 = v1;		//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = v2;		//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = u2;		//Outward velocity opposite solid c direction
		
		oppcoords.Pa_a1 = z2 - Pz;
		oppcoords.Pb_b1 = Py - y1;
		oppcoords.Pc_c1 = Px - x1;


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

			P->X = (double)(C->xp);
			P->VoxelX = C->xp;
			
			ReachesBorderX=true;
		}


		xe = oppcoords.Pc_c1 + x1;
		ye = oppcoords.Pb_b1 + y1;
		ze = z2 - oppcoords.Pa_a1;
		

	}else if(Solids&PositiveZSolid && Solids&NegativeZSolid && Solids&NegativeYSolid){	//Solid on z2,z1 &y1 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
			AdvectDebug << "Positive Y Solid, Negative Y Solid, Positive Z Solid" << endl; 
		}

		//a=y, b=x, c=-z
		oppcoords.Velocity_b1 = u1;		//Inward velocity in b direction, free of solids
		oppcoords.Velocity_b2 = u2;		//Outward velocity in b direction, free of solids
		oppcoords.Velocity_c2 = -w1;	//Outward velocity opposite solid c direction
		
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


	}else if(Solids&PositiveXSolid && Solids&NegativeXSolid && Solids&PositiveYSolid){		//Solid on x2,x1 &y2 sides
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Case 4 Opposite Solids:" << endl;
		AdvectDebug << "w1 = " << w1 << endl << "w2 = " << w2 << endl;
		AdvectDebug << "(Px,Py,Pz) = (" << Px << ", " << Py << ", " << Pz << ")" << endl;
		AdvectDebug << "(u,v,w) = (" << u << ", " << v << ", " << w << ")" << endl;
	}

	
	if(w==0 || (w>0 && w2<=0) || (w<0 && w1>=0)){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
		Tz = DBL_MAX;
	}else if(fabs(w2-w1) > Min_dV){	//Extension to the model on account of LB flux non-conservation
		if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
			Tz = log( w2/(w1 + (w2-w1)*Pz) ) / (36*(w2-w1)*(1-Px)*Px*(1-Py)*Py);
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 3" << endl; }
			Tz = log( w1/(w1 + (w2-w1)*Pz) ) / (36*(w2-w1)*(1-Px)*Px*(1-Py)*Py);
		}
	}else{
		if(w>0){
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 4" << endl; }
			Tz = (1-Pz)/w;
		}else{
			if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 5" << endl; }
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: ( [DBL_MAX], [DBL_MAX], " << Tz << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(Tz)||Tz<0){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(Tz)||Tz<0){
		Tz = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: ( [DBL_MAX], [DBL_MAX], " << Tz << ")" << endl;
	}

	double ze;

	if(Case->ReachedNegativeCBorder){
		ze = 0;
	}else if(Case->ReachedPositiveCBorder){
		ze = 1;
	}else if(fabs(w2-w1) > Min_dV){	//Extension to the model on account of LB flux non-conservation
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
		double _w1_w2w1 = w1/(w2-w1);
		ze = (Pz + _w1_w2w1)*exp( 36*(w2-w1)*Px*(1-Px)*Py*(1-Py)*t ) - _w1_w2w1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
		ze = Pz + w*t;
	}

	Case->Pc_c1 = ze;

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 1" << endl; }
		Ty = log( 1-Py )/(A*v1);						// = ( (dx*dx*dy*dy*dz)/(12*v1*(y2-Py)*(Pz-z1)*(x2-Px)*(Px-x1)) ) * log( (y2-Py)/dy );
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Ty 2" << endl; }
		Ty = DBL_MAX;
	}

	//Tau z
	if(w2>0){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 1" << endl; }
		Tz = -log( Pz )/(A*w2);							// = ( (dx*dx*dy*dz*dz)/(12*w2*(y2-Py)*(Pz-z1)*(x2-Px)*(Px-x1)) ) * log( dz/(Pz-z1) );
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "Tz 2" << endl; }
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

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: ( [DBL_MAX], " << Ty << ", " << Tz << ")" << endl;
	}

#if DebugAdvectionMode>=3
	if(isInfNaN(Ty)||Ty<0||isInfNaN(Tz)||Tz<0){
		TauError = true;
		return;
	}
#else
	if(isInfNaN(Ty)||Ty<0){
		Ty = DBL_MAX;
	}
	if(isInfNaN(Tz)||Tz<0){
		Tz = DBL_MAX;
	}
#endif

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Taus: ( [DBL_MAX], " << Ty << ", " << Tz << ")" << endl;
	}

	//Positions
	double ye, ze;

	if(Case->ReachedNegativeBBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 1" << endl; }
		ye = 0;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ye 2" << endl; }
		ye = 1 - (1 - Py) * exp( -A*v1*t );	// = y2 - (y2 - Py) * exp( -((12*v1*(Px-x1)*(x2-Px)*(y2-Py)*(Pz-z1))/(dx*dx*dy*dy*dz))*t );
	}

	if(Case->ReachedPositiveCBorder){
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 1" << endl; }
		ze = 1;
	}else{
		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){ AdvectDebug << "ze 2" << endl; }
		ze = Pz * exp( A*w2*t );				// = z1 + (Pz - z1) * exp( ((12*w2*(Px-x1)*(x2-Px)*(y2-Py)*(Pz-z1))/(dx*dx*dy*dz*dz))*t );
	}

	Case->Pb_b1 = ye;
	Case->Pc_c1 = ze;

	if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
		AdvectDebug << "Raw Position: (" << Px << ", " << ye << ", " << ze << ")" << endl;
	}

}

//Advection function used when four neighbouring voxels are solid (15 cases)
void CaseFourSolids(Particle* P, Coords* C, Voxel* LatticeElement, double* dt){

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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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

		if(DebugAdvectionMode>=2 && DebugParticle && OutputDebugData){
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