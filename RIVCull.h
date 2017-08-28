#include <stdio.h>
#include <stdlib.h>
#include <strsafe.h>
#include <windows.h>

#define GPUTHREADS 512
#define SEEDMASK 25214903917
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))
struct RIVData{
	int RIVsize;
	int nonZeros;
	long long int *masks;
	int *d_OpenSlot;
	int *d_SlotEnd;
	float *d_magnitudes;
	int thing =0;
}RIVKeyData;
typedef struct{
	char name[100];
	int *values;
	int *locations;
	int count;
	float magnitude;
	int boolean;
}sparseRIV;
static void HandleError(cudaError_t err, const char *file, int line){
	if(err !=cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		exit(EXIT_FAILURE);
	}
}

	
__global__ void squirt(float *d_magnitudes, int N){
	int id =(blockIdx.x*blockDim.x + threadIdx.x);
	if(id>=N) return;
	
	d_magnitudes[id] = sqrt(d_magnitudes[id]);
	
}
	
__global__ void generateLocations(int *d_seeds, long long int mask, int *d_locations, int RIVsize, int team, int seedCount, int nonZeros){

	int id =nonZeros*(blockIdx.x*blockDim.x + threadIdx.x)+team;
	if(id>=seedCount) return;
	d_locations[id] = ((d_seeds[id]^mask) & 2147483647) %(RIVsize);

}
__global__ void D2S( int* d_DenseRIV, int* d_SparseValues, int* d_SparseLocations, int *d_NZCount, int d_DenseSize){
	
	int id =(blockIdx.x*blockDim.x + threadIdx.x);
	if(id>=d_DenseSize) return;
	int value = *(d_DenseRIV+id);
	if(!value) return;
	
	int sparseSlot = atomicAdd(d_NZCount, 1);
	*(d_SparseValues+sparseSlot) = value;
	*(d_SparseLocations+sparseSlot) = id;
}
	
__global__ void S2D(int *d_locations, int *d_values, int *d_OpenSlot, int numberOfValues){
	int id = blockIdx.x*blockDim.x + threadIdx.x;
		if(id>=numberOfValues) return ;
		atomicAdd( d_OpenSlot + *(d_locations+id) , *(d_values+id));
		
}
__global__ void I2D(int *d_locations, int *d_OpenSlot, int numberOfValues){
	int id = blockIdx.x*blockDim.x + threadIdx.x;
		int value = (id%2) ? -1: 1;
		if(id>=numberOfValues) return ;
		atomicAdd( d_OpenSlot + *(d_locations+id) , value);
		
}
sparseRIV FileToL2(FILE *data);
sparseRIV compileD2S(int *denseInput);
void setKeyData(int RIVsize, int nonZeros);
int* mapS2D(sparseRIV input);
int* makeSparseLocations(int *seeds, int seedCount);
void makeSeeds(unsigned char* word, int **seeds, int *seedCount);
char *bin(void *shit, int size);
float* getCosSimsCPU(sparseRIV *inputs, int baseNumber, int RIVCount);
float* cosineCompare(sparseRIV baseRIV, sparseRIV *multipliers, int multiplierCount);
float *getMagnitudesCPU(sparseRIV *inputs, int RIVCount);
int *mapI2D(int *locations, int seedCount);


sparseRIV FileToL2(FILE *data){
	
	
	unsigned char *word = (unsigned char*)calloc(2000, 1);
	int *seeds = ( int*)malloc(RIVKeyData.nonZeros*sizeof( int));
	if(seeds == NULL) printf("malloc fail 146");
		
	
	int seedCount = 0;
	while(fscanf(data, "%s", word)){
		if(feof(data)){
			break;
		}
		if(!(*word)){
			break;
		}
		
		makeSeeds(word, &seeds, &seedCount);
		memset(word, 0, 2000);
		
	}
	
	int *locations = makeSparseLocations(seeds, seedCount);
	
	int *L2dense;
	L2dense = mapI2D(locations, seedCount);
	free(locations);
	
	sparseRIV output  = compileD2S(L2dense);	

	free(seeds);
	output.boolean = 1;
	return output;
}

float* getCosSimsCPU(sparseRIV *inputs, int baseNumber, int RIVCount){

	sparseRIV baseSparse = *(inputs+baseNumber);
	int *baseDenseRIV = mapS2D(baseSparse);
	float *outputs = (float*)malloc((RIVCount-(baseNumber+1))* sizeof(float));
	float *output_slider = outputs;
	sparseRIV *multipliers = inputs+baseNumber+1;
	sparseRIV *multipliersStop = inputs+RIVCount;
	float minsize = baseSparse.magnitude * .75;
	float maxsize = baseSparse.magnitude * 1.25;
	
	while(multipliers<multipliersStop){
		if(((*multipliers).boolean) && (((*multipliers).magnitude < maxsize) && ((*multipliers).magnitude > minsize))){
			
			int dot = 0;
			int *values = (*multipliers).values;
			int *locations = (*multipliers).locations;
			int *locationsStop = locations+(*multipliers).count;
			
			while(locations<locationsStop){
				dot += (*values)*(*(baseDenseRIV+(*locations)));
				locations++;
				values++;
			}
			//printf("%d\n", dot);
			*output_slider= dot/((baseSparse.magnitude)*((*multipliers).magnitude));
			if(*output_slider>=.98){
				printf("%s\t%s\t%f\n", (*multipliers).name, baseSparse.name, *output_slider);
				(*multipliers).boolean = 0;
				RIVKeyData.thing ++;
			}
		}
		multipliers++;
		output_slider++;
	}
	
	return outputs;
	
	
}
float* cosineCompare(sparseRIV baseRIV, sparseRIV *multipliers, int multiplierCount){

	int *baseDenseRIV = mapS2D(baseRIV);
	float *outputs = (float*)malloc((multiplierCount)* sizeof(float));
	float *output_slider = outputs;
	sparseRIV *multipliersStop = multipliers+multiplierCount;
	float minsize = baseRIV.magnitude * .75;
	float maxsize = baseRIV.magnitude * 1.25;
	
	while(multipliers<multipliersStop){
		if(((*multipliers).boolean) && (((*multipliers).magnitude < maxsize) && ((*multipliers).magnitude > minsize))){
			
			int dot = 0;
			int *values = (*multipliers).values;
			int *locations = (*multipliers).locations;
			int *locationsStop = locations+(*multipliers).count;
			
			while(locations<locationsStop){
				dot += (*values)*(*(baseDenseRIV+(*locations)));
				locations++;
				values++;
			}
			//printf("%d\n", dot);
			*output_slider= dot/((baseRIV.magnitude)*((*multipliers).magnitude));
			if(*output_slider>=.98){
				printf("%s\t%s\t%f\n", (*multipliers).name, baseSparse.name, *output_slider);
				(*multipliers).boolean = 0;
				RIVKeyData.thing ++;
			}
		}
		multipliers++;
		output_slider++;
	}
	
	return outputs;
	
	
}

float *getMagnitudesCPU(sparseRIV *inputs, int RIVCount){
	float *magnitudes = (float*) malloc(RIVCount*sizeof(float));
	float *magnitudes_slider = magnitudes;
	for(int i=0; i<RIVCount; i++){
		int temp = 0;
		int *values = inputs[i].values;
		int *values_stop = values+inputs[i].count;
		while(values<values_stop){
			temp += (*values)*(*values);
			values++;
		}
		*magnitudes_slider = temp;
		magnitudes_slider++;

	}
		HANDLE_ERROR (cudaMalloc((void**)&RIVKeyData.d_magnitudes, RIVCount*sizeof(float)));
		HANDLE_ERROR (cudaMemcpy (RIVKeyData.d_magnitudes, magnitudes, RIVCount*sizeof(float), cudaMemcpyHostToDevice));
		
		int blockSize;  
		int minGridSize = 0;
		int gridSize; 
		cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, squirt); 
		gridSize = ((RIVCount + blockSize -1) / blockSize)+1; 
		
		squirt<<<gridSize,blockSize >>> (RIVKeyData.d_magnitudes, RIVCount);

		HANDLE_ERROR (cudaMemcpy (magnitudes, RIVKeyData.d_magnitudes, RIVCount*sizeof(float), cudaMemcpyDeviceToHost));
		magnitudes_slider = magnitudes;
		for(int i=0; i<RIVCount; i++){
			inputs[i].magnitude =  *magnitudes_slider;
			magnitudes_slider++;
			
		}
	return magnitudes;	
}

int* mapS2D(sparseRIV input){
	int* valuesOut = (int*)malloc(RIVKeyData.RIVsize*sizeof(int));

	int *d_locations = RIVKeyData.d_OpenSlot+RIVKeyData.RIVsize;
	
	int *d_values = d_locations+input.count;
						

	HANDLE_ERROR (cudaMemset (RIVKeyData.d_OpenSlot, 0, RIVKeyData.RIVsize*sizeof(int)));
	HANDLE_ERROR (cudaMemcpy (d_locations, input.locations, input.count*sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR (cudaMemcpy (d_values, input.values, input.count*sizeof(int), cudaMemcpyHostToDevice));
	
	int blockSize;  
	int minGridSize = 0;
	int gridSize; 
	cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, S2D); 
	gridSize = ((input.count + blockSize -1) / blockSize)+1; 
		
	S2D <<<gridSize,blockSize>>> (d_locations, d_values, RIVKeyData.d_OpenSlot, input.count);
	
	HANDLE_ERROR (cudaMemcpy (valuesOut, RIVKeyData.d_OpenSlot, RIVKeyData.RIVsize*sizeof(int), cudaMemcpyDeviceToHost));
	return valuesOut;
	
}
int* mapI2D(int *locations, int valueCount){

	int* valuesOut = (int*)malloc(RIVKeyData.RIVsize*sizeof(int));

	int *d_locations = RIVKeyData.d_OpenSlot+RIVKeyData.RIVsize;
	
	HANDLE_ERROR (cudaMemset (RIVKeyData.d_OpenSlot, 0, RIVKeyData.RIVsize*sizeof(int)));
	HANDLE_ERROR (cudaMemcpy (d_locations, locations, valueCount*sizeof(int), cudaMemcpyHostToDevice));
	int blockSize;  
	int minGridSize = 0;
	int gridSize; 
	cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, I2D); 
	gridSize = ((valueCount + blockSize -1) / blockSize)+1; 
		
	I2D <<<gridSize,blockSize>>> (d_locations, RIVKeyData.d_OpenSlot, valueCount);
	
	HANDLE_ERROR (cudaMemcpy (valuesOut, RIVKeyData.d_OpenSlot, RIVKeyData.RIVsize*sizeof(int), cudaMemcpyDeviceToHost));
	return valuesOut;
	
}
sparseRIV compileD2S(int *denseInput){
	
	int *d_valueCount;
	HANDLE_ERROR (cudaMalloc((void**)&d_valueCount, sizeof(int)));
	HANDLE_ERROR(cudaMemset(d_valueCount, 0, sizeof(int)));
	
	HANDLE_ERROR (cudaMemcpy (RIVKeyData.d_OpenSlot, denseInput, RIVKeyData.RIVsize*sizeof(int), cudaMemcpyHostToDevice));
	int *d_outValues = RIVKeyData.d_OpenSlot+RIVKeyData.RIVsize;
	int *d_outLocations = d_outValues+RIVKeyData.RIVsize;
	int blockSize;  
	int minGridSize = 0;
	int gridSize; 
	cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, D2S); 
	
	gridSize = ((RIVKeyData.RIVsize + blockSize -1) / blockSize)+1; 
		D2S <<<gridSize,blockSize>>> (RIVKeyData.d_OpenSlot, d_outValues, d_outLocations, d_valueCount, RIVKeyData.RIVsize);
	cudaDeviceSynchronize();
	sparseRIV output;
	HANDLE_ERROR (cudaMemcpy (&output.count, d_valueCount, sizeof(int), cudaMemcpyDeviceToHost));

	output.values = (int*)malloc(output.count*sizeof(int));
	if(output.values ==NULL) printf("malloc fail 215");
	output.locations = (int*)malloc(output.count*sizeof(int));
	if(output.locations ==NULL) printf("malloc fail 217");
	HANDLE_ERROR (cudaMemcpy (output.values, d_outValues, (output.count)*sizeof(int), cudaMemcpyDeviceToHost));
	HANDLE_ERROR (cudaMemcpy (output.locations, d_outLocations, (output.count)*sizeof(int), cudaMemcpyDeviceToHost));
	cudaFree(d_valueCount);
	return output;
	
}


void setKeyData(int RIVsize, int nonZeros){
	RIVKeyData.RIVsize = RIVsize;
	RIVKeyData.nonZeros = nonZeros;
	RIVKeyData.masks = (long long int*)malloc(nonZeros*sizeof(long long int));
	
	for(int i = 0; i<nonZeros; i++){
		RIVKeyData.masks[i] = SEEDMASK>>(5*i);
	}
	//HANDLE_ERROR (cudaMalloc((void**)&RIVKeyData.d_masks, nonZeros*sizeof(long long int)));
	//HANDLE_ERROR (cudaMemcpy(RIVKeyData.d_masks, masks, nonZeros*sizeof(int), cudaMemcpyHostToDevice));
	
	HANDLE_ERROR (cudaMalloc((void**)&RIVKeyData.d_OpenSlot, 100000000*sizeof(int)));
	RIVKeyData.d_SlotEnd = RIVKeyData.d_OpenSlot+100000000;
	RIVKeyData.thing = 0;
	
}
void makeSeeds(unsigned char* word,  int **seeds, int *seedCount){
	
	int i=0;
	int seedbase = 0;
	(*seeds) = ( int*)realloc((*seeds), (*seedCount+RIVKeyData.nonZeros)*sizeof(int));
	while(*word){
		seedbase += (*(word))<<(i*5);
		word++;
		i++;
		
	}
	int *seedTrack = (*seeds)+*seedCount;
	for(i =0 ; i<RIVKeyData.nonZeros; i++){
		
		*seedTrack = (seedbase>>i)+(3*i);
		seedTrack++;
	
	}
	*seedCount+=RIVKeyData.nonZeros;
	return;
}

int* makeSparseLocations(int* seeds, int seedCount){

	int *d_locations = RIVKeyData.d_OpenSlot;
	int *d_seeds = d_locations+seedCount;
	HANDLE_ERROR (cudaMemcpy(d_seeds, seeds, seedCount*sizeof(int), cudaMemcpyHostToDevice));
	
	int blockSize;  
	int minGridSize = 0;
	int gridSize; 

	cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, generateLocations); 
	gridSize = ((seedCount + blockSize -1) / (RIVKeyData.nonZeros*blockSize))+1; 
	long long int *mask = RIVKeyData.masks;
	for(int team=0; team<RIVKeyData.nonZeros; team++){
		generateLocations <<<gridSize,blockSize>>> (d_seeds, *mask, d_locations, RIVKeyData.RIVsize, team, seedCount, RIVKeyData.nonZeros);
		mask++;
	}
	
	cudaDeviceSynchronize();
	int *locations = (int*) malloc(seedCount*sizeof(int));
	HANDLE_ERROR (cudaMemcpy(locations, d_locations, seedCount*sizeof(int), cudaMemcpyDeviceToHost));
	return locations;
}
int* makeSparseValues(int seedCount){
	int *values = (int*)malloc(seedCount * sizeof(int));
	if(values ==NULL) printf("malloc fail 330");
	for(int i=0; i<seedCount; i++){
		(i%2)? (*(values+i) = -1) : (*(values+i) = 1);
	}
	return values;
	
	
}

