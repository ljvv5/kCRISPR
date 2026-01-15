#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main()
{
	FILE *read_state, *ratio_state;
	read_state=fopen("time_structure.dat","r+");
	ratio_state=fopen("ratio_state.dat","w+");
	int state;
	float t;
	int time_fold[1000][20];
	int time_total[1000];
	memset (time_fold, 0, sizeof(int)*1000*20);
	memset (time_total, 0, sizeof(int)*1000);
	int i=0;
	while(!feof(read_state))
	{
		fscanf(read_state, "%f %d\n", &t, &state);
		i=0;
		while(t*10.>(float)i)
		{
			i++;
			if(i==999)
				break;
		}
		time_fold[i][state-1]+=1;
		
	}//return 0;
	for(int m=0;m<1000;m++)
	{
		for(int n=0;n<20;n++)
		{
			time_total[m]+=time_fold[m][n];
		
		}
		
	}
	
	int total_s=0;// total_n=0;
	for(int m=0;m<1000;m++)
	{
		fprintf(ratio_state,"%f ",(float)m/10.);
		for(int n=0;n<20;n++)
		{
			total_s+=time_fold[m][n];
			
			//if(fmod(n+1,4)==0&&n!=0)
			{
				fprintf(ratio_state,"%f ",(float)total_s/(float)time_total[m]);
				total_s=0;
			}
			//fprintf(ratio_state,"%f ",(float)time_fold[m][n]/(float)time_total[m]);
		
		}
		fprintf(ratio_state,"\n");
	}
	fclose(read_state);
	fclose(ratio_state);
	return 0;
}
