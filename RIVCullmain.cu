#include <stdio.h>
#include <stdlib.h>
#include "RIVCull.h"
int main(int argc, char argv[]){
	clock_t beginprep = clock();
	HANDLE hfind = NULL;
	WIN32_FIND_DATA findFile;
	int fileCount = 0;
	setKeyData(25000, 4);
	sparseRIV *fileRIVs = (sparseRIV*) malloc(7000*sizeof(sparseRIV));	
	//switch to dynamic
	char pathString[2000];
	sprintf(pathString, "%s\\*.*","0");
	if((hfind = FindFirstFile(pathString, &findFile)) == INVALID_HANDLE_VALUE){
		printf("invalid path: %s\n", pathString);
		return 1;
	}
	int maxSize = 0;
	do{
		if(strcmp(findFile.cFileName, ".") != 0 && 
		strcmp(findFile.cFileName, "..") != 0){
			sprintf(pathString, "%s\\%s", "0", findFile.cFileName);
			FILE *input = fopen(pathString, "r");
			//printf("%s\n", pathString);
			fileRIVs[fileCount] = FileToL2(input);
			strcpy(fileRIVs[fileCount].name, pathString);
			if(fileRIVs[fileCount].count>maxSize){
				maxSize = fileRIVs[fileCount].count;
			}
			/*for(int i=0; i<fileRIVs[fileCount].count; i++){
				printf("%d, %d\n", fileRIVs[fileCount].locations[i], fileRIVs[fileCount].values[i]);
			}*/
			fclose(input);
			fileCount++;
		}
	}while(FindNextFile(hfind, &findFile));
	

	FindClose(hfind);
	float *magnitudes = getMagnitudesCPU(fileRIVs, fileCount); 

	clock_t beginnsquared = clock();

	float **cosSims = (float**) malloc(fileCount*sizeof(float*));
	for(int i=0; i<fileCount; i++){
		if(fileRIVs[i].boolean){
			*cosSims = cosineCompare(fileRIVs[i], fileRIVs[i+1], fileCount-(i+1));
		}
		cosSims++;
		
	}
	clock_t endnsquared = clock();
	double time = (double)(endnsquared - beginnsquared) / CLOCKS_PER_SEC;
	printf("nsquared time:%lf\n\n", time);
	printf("%d <", RIVKeyData.thing);
		clock_t endprep = clock();
	double time_spent = (double)(endprep - beginprep) / CLOCKS_PER_SEC;
	printf("total time:%lf\n\n", time_spent);
return 0;
}
