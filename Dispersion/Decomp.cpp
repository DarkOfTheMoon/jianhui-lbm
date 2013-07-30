//MPI Multithreading options
bool MPIDecomposeX = true;					//Whether or not to decompose the lattice in the X direction for MPI mode
bool MPIDecomposeY = true;					//Whether or not to decompose the lattice in the Y direction for MPI mode
bool MPIDecomposeZ = false;					//Whether or not to decompose the lattice in the Z (flow) direction for MPI mode

//Lattice
//short FlowDirection = 2;					//X = 0, Y = 1, Z = 2 for Propagator

short NLattice_x = 320;						//Lattice of NLattice_x x NLattice_y x NLattice_z
short NLattice_y = 320;
short NLattice_z = 640;

int BoundaryConditionX = 0;					//Boundary conditions in X,Y,Z. Loop boundary = 0, Solid boundary = 1
int BoundaryConditionY = 0;
int BoundaryConditionZ = 0;

//////////////

struct CPUDomain{
	short x0;				//Indeces of CPU domain voxel range
	short x1;				//Index of last voxel in x range
	short y0;
	short y1;
	short z0;
	short z1;

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

	int ThreadID;				//Thread ID associated with this CPU Domain
};

//Lattice Decomposition declarations

char* _DecompositionLattice;		//Temporary geometry array

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

/////////////////////////////////////////

void _DecompositionReadSolids(){	//Reads solid nodes from file

	FILE* InFile = fopen(InputFilePath(SolidsFileName), "rb");
	
	if(InFile!=0){
		dcout << "Reading in solids from '" << SolidsFileName << "'" << endl;
	}else{
		dcout << "Could not open solids file" << endl;
		return;
	}

	__int64 flength;

#ifdef WIN32
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

	while(true){

		if(readpos == readl){
			readl = buffersize;

			if(buffersize > (flength-nRead)){
				readl = flength-nRead;
			}
			if(readl == 0){
				dcout << "Error: reached end of geometry file before lattice entirely read" << endl;
				goto DecompositionReadSolidsEnd;
			}
			fread(buf,1,readl,InFile);
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
						goto DecompositionReadSolidsEnd;
					}
				}
			}

		}

		readpos++;
	}

DecompositionReadSolidsEnd:
	delete[] buf;
	fclose(InFile);
	return;
}



//Decomposes lattice and sends data to all threads
int LatticeDecomposition(int NProcessors, int ThreadID){

	ofstream DecompositionLog;
	DecompositionLog.open(OutputFilePath(DecompositionLogFile));

	int NArray = ((NLattice_x*NLattice_y*NLattice_z)>>3) + 1;
	_DecompositionLattice = new char[NArray];

	memset(_DecompositionLattice,0,NArray);
	_DecompositionReadSolids();

	int nZ;			//Number of partitions in Z

	if(MPIDecomposeX && MPIDecomposeY && MPIDecomposeZ){			//Decompose all three directions

		double ycoefficient = (double)NLattice_y / (double)NLattice_x;
		double zcoefficient = (double)NLattice_z / (double)NLattice_x;
		double N = pow((double)NProcessors / (ycoefficient*zcoefficient) , 1.0/3.0);

		if(FlowDirection!=2){
			nZ = (int)ceil( zcoefficient*N );
		}else{
			nZ = (int)floor( zcoefficient*N );		//Prefer division in non-flow directions
		}

	}else if(MPIDecomposeX && !MPIDecomposeY && MPIDecomposeZ){		//Decompose X and Z

		double zcoefficient = (double)NLattice_z / (double)NLattice_x;
		double N = pow((double)NProcessors / zcoefficient , 1.0/2.0);

		if(FlowDirection!=2){
			nZ = (int)ceil( zcoefficient*N );
		}else{
			nZ = (int)floor( zcoefficient*N );		//Prefer division in non-flow directions
		}

	}else if(!MPIDecomposeX && MPIDecomposeY && MPIDecomposeZ){		//Decompose Y and Z

		double zcoefficient = (double)NLattice_z / (double)NLattice_y;
		double N = pow((double)NProcessors / zcoefficient , 1.0/2.0);

		if(FlowDirection!=2){
			nZ = (int)ceil( zcoefficient*N );
		}else{
			nZ = (int)floor( zcoefficient*N );		//Prefer division in non-flow directions
		}

	}else if(!MPIDecomposeX && !MPIDecomposeY && MPIDecomposeZ){	//Decompose Z only

		nZ = NProcessors;

	}else if(!MPIDecomposeZ){										//No Z decomposition

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

	if(MPIDecomposeZ && nZ>1){

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

					if(Plane[i].UpperRound){
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

				if(Plane[i].UpperRound){
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

		int nY;		//Number of blocks in y

		if(MPIDecomposeX && MPIDecomposeY){			//Decompose X and Y (Z doesn't matter here)

			double ycoefficient = (double)NLattice_y / (double)NLattice_x;
			double N = sqrt( (double)Plane[zBlock].NSquares/ycoefficient );

			nY = (int)ceil( ycoefficient*N );

		}else if(!MPIDecomposeX){					//No X decomposition

			nY = Plane[zBlock].NSquares;

		}else if(!MPIDecomposeY){					//No Y decomposition

			nY = 1;

		}

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
		
		if(MPIDecomposeY && nY>1){

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

						if(Plane[zBlock].YBlocks[0].UpperRound){
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

						if(Plane[zBlock].YBlocks[i].UpperRound){
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

					if(Plane[zBlock].YBlocks[i].UpperRound){
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

		}else{		//Not decomposing Y

			Plane[zBlock].YBlocks[0].y0Rounded = 0;
			Plane[zBlock].YBlocks[0].y1Rounded = NLattice_y;

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

			if(MPIDecomposeX && nX>1){

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

							if(Block->XBlocks[i].UpperRound){
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

						if(Block->XBlocks[i].UpperRound){
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

			}else{

				Block->XBlocks[0].x0Rounded = 0;
				Block->XBlocks[0].x1Rounded = NLattice_x;

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

		Domains[i].ThreadID = i;

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

			if(Domains[i2].x1 == Cn.xn){
				if(((Domains[i2].y0>=Domains[i].y0&&Domains[i2].y0<=Domains[i].y1)||(Domains[i2].y1>=Domains[i].y0&&Domains[i2].y1<=Domains[i].y1)||(Domains[i2].y0<Domains[i].y0&&Domains[i2].y1>Domains[i].y1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsXn[Domains[i].nNeighbourDomainsXn] = i2;
					Domains[i].nNeighbourDomainsXn++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in xn" << endl;

				}
			}

			if(Domains[i2].x0 == Cp.xp){
				if(((Domains[i2].y0>=Domains[i].y0&&Domains[i2].y0<=Domains[i].y1)||(Domains[i2].y1>=Domains[i].y0&&Domains[i2].y1<=Domains[i].y1)||(Domains[i2].y0<Domains[i].y0&&Domains[i2].y1>Domains[i].y1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsXp[Domains[i].nNeighbourDomainsXp] = i2;
					Domains[i].nNeighbourDomainsXp++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in xp" << endl;
				}
			}

			if(Domains[i2].y1 == Cn.yn){
				if(((Domains[i2].x0>=Domains[i].x0&&Domains[i2].x0<=Domains[i].x1)||(Domains[i2].x1>=Domains[i].x0&&Domains[i2].x1<=Domains[i].x1)||(Domains[i2].x0<Domains[i].x0&&Domains[i2].x1>Domains[i].x1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsYn[Domains[i].nNeighbourDomainsYn] = i2;
					Domains[i].nNeighbourDomainsYn++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in yn" << endl;
				}
			}

			if(Domains[i2].y0 == Cp.yp){
				if(((Domains[i2].x0>=Domains[i].x0&&Domains[i2].x0<=Domains[i].x1)||(Domains[i2].x1>=Domains[i].x0&&Domains[i2].x1<=Domains[i].x1)||(Domains[i2].x0<Domains[i].x0&&Domains[i2].x1>Domains[i].x1))
					&& ((Domains[i2].z0>=Domains[i].z0&&Domains[i2].z0<=Domains[i].z1)||(Domains[i2].z1>=Domains[i].z0&&Domains[i2].z1<=Domains[i].z1)||(Domains[i2].z0<Domains[i].z0&&Domains[i2].z1>Domains[i].z1))){
					
					Domains[i].NeighbourDomainsYp[Domains[i].nNeighbourDomainsYp] = i2;
					Domains[i].nNeighbourDomainsYp++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in yp" << endl;
				}
			}

			if(Domains[i2].z1 == Cn.zn){
				if(((Domains[i2].x0>=Domains[i].x0&&Domains[i2].x0<=Domains[i].x1)||(Domains[i2].x1>=Domains[i].x0&&Domains[i2].x1<=Domains[i].x1)||(Domains[i2].x0<Domains[i].x0&&Domains[i2].x1>Domains[i].x1))
					&& ((Domains[i2].y0>=Domains[i].y0&&Domains[i2].y0<=Domains[i].y1)||(Domains[i2].y1>=Domains[i].y0&&Domains[i2].y1<=Domains[i].y1)||(Domains[i2].y0<Domains[i].y0&&Domains[i2].y1>Domains[i].y1))){
					
					Domains[i].NeighbourDomainsZn[Domains[i].nNeighbourDomainsZn] = i2;
					Domains[i].nNeighbourDomainsZn++;
					DecompositionLog << "\tAdjoins [" << i2 << "] in zn" << endl;
				}
			}

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
			dcout << "Memory allocation failed at (0): stopping execution" << endl;
			DecompositionLog << "Memory allocation failed at (0)" << endl;
			delete[] FilePointers;
			fclose(fin);
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
			dcout << "Memory allocation failed at thread [" << i << "]: stopping execution" << endl;
			delete[] Geometry;
			delete[] FilePointers;
			fclose(fin);
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
		dcout << "Sent geometry to thread " << i << endl;

		delete[] Geometry;

		i++;
	}

	//Thread 0 domain

	if(ThreadID==0){

		memcpy(&ProcessorDomain,&(Domains[0]),sizeof(CPUDomain));

		int nX = Domains[0].x1 - Domains[0].x0 + 1;
		int nY = Domains[0].y1 - Domains[0].y0 + 1;
		int nZ = Domains[0].z1 - Domains[0].z0 + 1;

		int nVoxels = nX*nY*nZ;

		LatticeBase = new unsigned int[nVoxels];

		_CPUDomainWidthX = nX;
		_CPUDomainWidthY = nY;
		
		//Geometry

		int nNonSolids = 0;

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
			dcout << "Memory allocation failed at thread [0]: stopping execution" << endl;
			delete[] FilePointers;
			fclose(fin);
			return 1;
		}

		Lattice[0].Solid = VoxelSolid;
		Lattice[0].Vx = 0;
		Lattice[0].Vy = 0;
		Lattice[0].Vz = 0;

		int LatticeIndex = 1;

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
		dcout << "Sent geometry to thread 0" << endl;

	}

	//

	//Free memory

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

	DecompositionLog.close();

	DecompositionInitialiseLattice(NProcessors, ThreadID);	//Initialise lattice velocities

	return 0;
}
