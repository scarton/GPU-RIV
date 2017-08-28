
__global__ void dotProduct(int *d_baseDenseRIV,
 							int *d_multiplierBlock, int *d_multiplierValueCount, 
							int *d_output, int multiplierCount, int displacement){
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if (id>=multiplierCount) return;
	int *d_multiplierLocations = d_multiplierBlock +(id*2*displacement);
	int *d_multiplierValues = d_multiplierLocations+displacement;
	//int *d_baseStop = d_baseLocations+baseValueCount;
	int *d_multiplierStop = d_multiplierLocations+(d_multiplierValueCount[id]);
	//printf("spacing: %d, %d\n", (d_multiplierStop-d_multiplierLocations), d_multiplierValueCount[id]);
	d_output+=id;
	*d_output= 0;	
	while(d_multiplierLocations< d_multiplierStop){
		*d_output += (*d_multiplierValues)*(d_baseDenseRIV[*d_multiplierLocations]);
		d_multiplierValues++;
		d_multiplierLocations++;
		
	}
}
__global__ void getMagnitude(float *d_magnitudes, int *d_values, int *valueCount, int RIVCount, int memSectionSize){
	//consider changing to single operation per thread
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if (id>=RIVCount) return;
	d_values+=(id*(memSectionSize));
	int *stop = d_values+ valueCount[id];
	float *magnitude = d_magnitudes+id;
	*magnitude = 0;
	for( ;d_values<stop; d_values++){
		*magnitude +=(*d_values)*(*d_values);
	}
	
	*magnitude = sqrt(*magnitude);
	
}
int* getDotProducts(sparseRIV *inputs, int baseNumber, int RIVCount, int maxSize){
	
	int remainingSet = RIVCount-(baseNumber+1);
	int *output = (int*)malloc(remainingSet*sizeof(int));
	int *output_slider=output;
	//badly written, fix it
	int *d_baseDenseRIV = RIVKeyData.d_OpenSlot;
	int *baseDenseRIV = mapS2D(inputs[baseNumber]);//as a byproduct, also places a denseRIV form of the input into RIVKeyData.d_OpenSlot at the beginning;
	
	HANDLE_ERROR (cudaMemcpy (d_baseDenseRIV, baseDenseRIV, RIVKeyData.RIVsize*sizeof(int), cudaMemcpyHostToDevice));
	int *d_slider = RIVKeyData.d_OpenSlot+RIVKeyData.RIVsize;

	
	int *valueCounts = (int*)malloc(remainingSet*sizeof(int));
	
	int *d_valueCounts = d_slider;
	d_slider+=remainingSet;
	int *d_output = d_slider;
	d_slider+=remainingSet;
	int *d_multiplierBlock = d_slider;
	
	int i=baseNumber+1;
	while(i<RIVCount){
		int doneSoFar = i;
		while((i<RIVCount) && (d_slider < RIVKeyData.d_SlotEnd)){
			if(inputs[i].boolean){
				//each set of locations and then values is layed out linear in GPU ram with buffer the size of the largest RIV
				HANDLE_ERROR (cudaMemcpy (d_slider, inputs[i].locations, inputs[i].count*sizeof(int), cudaMemcpyHostToDevice));
				d_slider +=maxSize;
				HANDLE_ERROR (cudaMemcpy (d_slider, inputs[i].values, inputs[i].count*sizeof(int), cudaMemcpyHostToDevice));
				d_slider +=maxSize;
				//printf("%p", d_slider);
				valueCounts[i-doneSoFar] = inputs[i].count;
			}else{
				valueCounts[i-doneSoFar] = 0;
			}
			i++;
		}
		int thisBlock = i-doneSoFar;
		HANDLE_ERROR (cudaMemcpy (d_valueCounts, valueCounts, thisBlock*sizeof(int), cudaMemcpyHostToDevice));
		//d_slider+= remainingSet;

		int blockSize;  
		int minGridSize = 0;
		int gridSize; 
		cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, dotProduct); 
		gridSize = ((thisBlock + blockSize -1) / blockSize)+1; 

			dotProduct<<<gridSize,blockSize>>>(d_baseDenseRIV, d_multiplierBlock, d_valueCounts, 
																d_output, thisBlock, maxSize);
		
		
		HANDLE_ERROR (cudaMemcpy (output_slider, d_output, thisBlock*sizeof(int), cudaMemcpyDeviceToHost));
		//printf("did a thing");
		output_slider+=thisBlock;
		d_slider = d_multiplierBlock;
		
	
	}
	return output;
	
}
float* getMagnitudes(sparseRIV *dataSet, int RIVCount, int maxSize){
	//int **values = (int**)malloc(RIVCount*sizeof(int*));
	/*int **d_values;
	HANDLE_ERROR (cudaMalloc((void***)&d_values, RIVCount*sizeof(int*)));
	*/
	int *valueCounts = (int*)malloc(RIVCount*sizeof(int));
	float *magnitudes = (float*)malloc(RIVCount*sizeof(float));
	HANDLE_ERROR (cudaMalloc((void**)&RIVKeyData.d_magnitudes, RIVCount*sizeof(float)));
	//HANDLE_ERROR(cudaMemset(RIVKeyData.d_magnitudes, 0, RIVCount*(sizeof(float))));
	float *magnitudes_slider = magnitudes;
	
	int *d_slider = RIVKeyData.d_OpenSlot;//+RIVCount;
	
	int *d_valueCounts = d_slider;
	d_slider+=RIVCount;
	int *d_valuesBlock = d_slider;
	//printf("magnitudesSlot: %d, %d\n", d_slider-RIVKeyData.d_OpenSlot, RIVCount*sizeof(float));
	//prepare for overflow?
	//printf("%d\n", RIVCount);
	int i=0;
	while(i<RIVCount){
		int doneSoFar = i;
		
		while((d_slider<RIVKeyData.d_SlotEnd) && (i<RIVCount)){
			HANDLE_ERROR (cudaMemcpy (d_slider, dataSet[i].values, dataSet[i].count*sizeof(int), cudaMemcpyHostToDevice));
			d_slider +=maxSize;
			valueCounts[i-doneSoFar] = dataSet[i].count;
			i++;
			//printf("%d, %d, %d, %d\n", d_slider-RIVKeyData.d_OpenSlot, doneSoFar, i, RIVCount);
		}
		int thisBlock = i-doneSoFar;
		
		//HANDLE_ERROR (cudaMemcpy (d_values, values, RIVCount*sizeof(int*), cudaMemcpyHostToDevice));
		HANDLE_ERROR (cudaMemcpy (d_valueCounts, valueCounts, thisBlock*sizeof(int), cudaMemcpyHostToDevice));
		int blockSize;  
		int minGridSize = 0;
		int gridSize; 
		cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, getMagnitude); 

		
		//HANDLE_ERROR (cudaMalloc((void**)&d_magnitudes, RIVCount*sizeof(float)));
		
		
		gridSize = ((thisBlock + blockSize -1) / blockSize)+1; 
			
			
			getMagnitude<<<gridSize,blockSize>>>(RIVKeyData.d_magnitudes, d_valuesBlock, d_valueCounts, thisBlock, maxSize);
		//printf("got here");
		HANDLE_ERROR (cudaMemcpy (magnitudes_slider, RIVKeyData.d_magnitudes, thisBlock*sizeof(float), cudaMemcpyDeviceToHost));
		magnitudes_slider+=thisBlock;
		d_slider =d_valuesBlock;
	}
	for(int i=0; i<RIVCount; i++){
		dataSet[i].magnitude = magnitudes[i];
		//printf("%f\n", dataSet[i].magnitude);
	}
	return magnitudes;
}
sparseRIV compileD2SOrdered(denseRIV input){
	
	//int *valueCount;
	//*RIVsize = 0;
	int *d_valueCount;
	HANDLE_ERROR(cudaMalloc((void**)&d_valueCount, sizeof(int)));
	HANDLE_ERROR(cudaMemset(d_valueCount, 0, sizeof(int)));
	int *d_locations = RIVKeyData.d_OpenSlot+RIVKeyData.RIVsize;
	//HANDLE_ERROR (cudaMemcpy (d_valueCount, valueCount, sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR (cudaMemcpy (RIVKeyData.d_OpenSlot, input.values, RIVKeyData.RIVsize*sizeof(int), cudaMemcpyHostToDevice));
	int blockSize;  
	int minGridSize = 0;
	int gridSize; 
	cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, D2SLocations); 
	
	gridSize = ((RIVKeyData.RIVsize + blockSize -1) / blockSize)+1; 
		D2SLocations <<<gridSize,blockSize>>> (RIVKeyData.d_OpenSlot, d_locations, d_valueCount, RIVKeyData.RIVsize);
		cudaDeviceSynchronize();
	sparseRIV output;
	HANDLE_ERROR (cudaMemcpy (&output.count, d_valueCount, sizeof(int), cudaMemcpyDeviceToHost));

	output.values = (int*)malloc(output.count*sizeof(int));
	if(output.values ==NULL) printf("malloc fail 246");
	output.locations = (int*)malloc(output.count*sizeof(int));
	if(output.locations ==NULL) printf("malloc fail 248");
	HANDLE_ERROR (cudaMemcpy (output.locations, d_locations, (output.count)*sizeof(int), cudaMemcpyDeviceToHost));
	qsort(output.locations, output.count, sizeof(int), compareLocations);
	
	
	for(int i=0; i<output.count; i++){
		output.values[i] = input.values[output.locations[i]];
	}
	free(input.values);
	cudaFree(d_valueCount);
	return output;
	
}
int compareLocations(const void *first, const void *second){
	int *f = (int*)first;
	int *s = (int*)second;
	return(*f - *s);
}
__global__ void D2SLocations(int *d_DenseRIV, int* d_SparseLocations, int* d_NZCount, int d_DenseSize){
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id>=d_DenseSize) return;
	if(!d_DenseRIV[id]) return;
	int sparseSlot = atomicAdd(d_NZCount, 1);
	d_SparseLocations[sparseSlot] = id;
}	
