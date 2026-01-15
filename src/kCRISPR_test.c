#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define k0 1

int main()
{
	FILE *out_time_structure;
	
	out_time_structure=fopen("time_structure.dat", "w+");
	float k_u[20];
	for(int i=0;i<20;i++)
	{
		k_u[i]=0.5*k0;
	}
	int state;
	k_u[0]=0.5;k_u[1]=0.9;k_u[2]=0.8;k_u[3]=0.5;k_u[4]=0.6;k_u[5]=0.4;k_u[6]=0.8;k_u[7]=0.7;k_u[8]=0.5;k_u[9]=0.6;k_u[10]=0.5;k_u[11]=2.1;k_u[12]=1.2;k_u[13]=1.1;k_u[14]=2.1;k_u[15]=0.9;k_u[16]=1.1;k_u[17]=0.9;k_u[18]=0.7;k_u[19]=0.9;k_u[20]=2.2;
	state=1;
	float r1,r2;
	float k_total;
	float r1_k;
	
	float t=0.0;
	
	srand((unsigned)time(NULL));
	for(int m=0;m<1000;m++)
	{
		t=0.0; state=1;
	
	for(int n=0;n<2000;n++)
	{
		if(n==0)
		{
			fprintf(out_time_structure, "%f %d\n", t, state);
		}
		
		r1=(double)(rand()/(double)RAND_MAX);
		r2=(double)(rand()/(double)RAND_MAX);
		
		if(state==1)
		{
			k_total=k0;//+k0;
		}
		else if(state==20)
		{
			k_total=k_u[state-1];//+k0;
		}
		else
		{
			k_total=k0+k_u[state-1];//+k0;
		}
		r1_k=r1*k_total;
		if(state!=20)
		{
			if(r1_k<k0)
			{	
				state=state+1;
			}
			/*else if(r1_k<k0*2)
			{
				state=state;
			}*/
			else
			{
				state=state-1;
			}
		}
		else
		{
			if(r1<k_u[state-1])
			{
				state=state-1;
			}
			else
			{
				state=state;
			}
		}
		
		t=t-(log(r2)/k_total);
		
		fprintf(out_time_structure, "%f %d\n", t, state);
	}
	}
	fclose(out_time_structure);

	return 0;
}

