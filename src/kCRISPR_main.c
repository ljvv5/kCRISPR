#include<stdio.h>  
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#define KMC_MAXT 1000
#define RUN_MAX 1000
#define K0 1.0
#define kbT 2.0 

typedef struct
{
	int position;
	char stacking_type[3];
	float energy;
}	Paraset;


typedef struct S_node
{
	int step_mc;
	char sequence_ntDNA[50];
	char sequence_sgRNA[50];
	int base_pairs[50];
	float energy_dsDNA;
	float energy_RDNA;
	float energy_ssRNA;
	float energy_ssDNA;
	S_node *neighbor[30];
}	structure_node;

typedef struct P_node
{
	structure_node *structure_node_C;
	P_node *node_before;
	P_node *node_next;
	P_node *last;
}	path_node;

extern int sttype();

Paraset Paraset_on_stenergy[30][16];
Paraset Paraset_off_stenergy[30][16];
extern void read_para();

char seq_tgDNA[30], seq_ntDNA[30], seq_sgRNA[30];
extern void set_seq();

path_node *path_all[2000];
extern void KMC_initialation();
extern void NODE_search();
extern void KMC_simulation();
extern void Data_output();

extern void KMC_run();

extern float calculate_rate(structure_node *structure_1);
extern float calculate_energy(structure_node *structure_1);


int main()
{
	FILE *Pop_step, *Pop_fin, *Path_step, *Path_avg;
	
	//Initialization//	
	read_para();
	set_seq();
		
	//KMC_simulation();
	
	//data_output();
	
	return 0;
}

int sttype(char type[3])
{
	if(!strcmp(type,"AA"))
		return 1;
	else if(!strcmp(type,"AC"))
		return 2;
	else if(!strcmp(type,"AG"))
		return 3;
	else if(!strcmp(type,"AT"))
		return 4;
	else if(!strcmp(type,"CA"))
		return 5;
	else if(!strcmp(type,"CC"))
		return 6;
	else if(!strcmp(type,"CG"))
		return 7;
	else if(!strcmp(type,"CT"))
		return 8;
	else if(!strcmp(type,"GA"))
		return 9;
	else if(!strcmp(type,"GC"))
		return 10;
	else if(!strcmp(type,"GG"))
		return 11;
	else if(!strcmp(type,"GT"))
		return 12;
	else if(!strcmp(type,"TA"))
		return 13;
	else if(!strcmp(type,"TC"))
		return 14;
	else if(!strcmp(type,"TG"))
		return 15;
	else if(!strcmp(type,"TT"))
		return 16;
	else return 0;
}

void read_para()
{
	FILE *PARA_st_energy;
	int pos;
	char type0[6], type[3];
	float En;
	
	printf("Reading On target Energy Parameters from parameters-ontarget.dat\n");
	PARA_st_energy=fopen("parameters-ontarget.dat","r+");
	while(!feof(PARA_st_energy))
	{
		memset (type0, 0, sizeof(type0));
		memset (type, 0, sizeof(type));
		En=0; pos=0;
		fscanf(PARA_st_energy, "%f %s\n", &En, type0);
		type[0]=type0[0];	
		type[1]=type0[1];
		if(type0[3]=='-')
			pos=0-((int)type0[4]-48);
		else
			pos=((int)type0[3]-48)*pow(10,(strlen(type0+3)-1))+(int)type0[4]-48*((strlen(type0+3)-1));
		if(sttype(type)==0)
		{
			printf("Stacking type error.");
			exit (-1);
		}
		Paraset_on_stenergy[pos+3][sttype(type)-1].position=pos;
		strcpy(Paraset_on_stenergy[pos+3][sttype(type)-1].stacking_type,type);
		Paraset_on_stenergy[pos+3][sttype(type)-1].energy=En;
	}
	fclose(PARA_st_energy);

	printf("Reading Off target Energy Parameters from parameters-offtarget.dat\n");
	PARA_st_energy=fopen("parameters-offtarget.dat","r+");
	while(!feof(PARA_st_energy))
	{
		memset (type0, 0, sizeof(type0));
		memset (type, 0, sizeof(type));
		En=0; pos=0;
		fscanf(PARA_st_energy, "%f %s\n", &En, type0);
		type[0]=type0[0];	
		type[1]=type0[1];
		if(type0[3]=='-')
			pos=0-((int)type0[4]-48);
		else
			pos=((int)type0[3]-48)*pow(10,(strlen(type0+3)-1))+(int)type0[4]-48*((strlen(type0+3)-1));
		if(sttype(type)==0)
		{
			printf("Stacking type error.");
			exit (-1);
		}
		Paraset_off_stenergy[pos-1][sttype(type)-1].position=pos;
		strcpy(Paraset_off_stenergy[pos-1][sttype(type)-1].stacking_type,type);
		Paraset_off_stenergy[pos-1][sttype(type)-1].energy=En;
	}
	fclose(PARA_st_energy);

	for(int i=0; i<30; i++)
		for(int j=0; j<16; j++)
		{
			if(Paraset_off_stenergy[i][j].stacking_type[0]==NULL) continue;
			printf("%d %s %f\n", Paraset_off_stenergy[i][j].position,Paraset_off_stenergy[i][j].stacking_type, Paraset_off_stenergy[i][j].energy);
		}
		
		
		
}

void set_seq()
{
	memset(seq_tgDNA, 0, sizeof(seq_tgDNA));
	memset(seq_ntDNA, 0, sizeof(seq_ntDNA));
	memset(seq_sgRNA, 0, sizeof(seq_sgRNA));
	FILE *SEQ;
	printf("Reading tgDNA Sequence from tgDNA.fasta\n");
	SEQ=fopen("tgDNA.fasta","r+");
	fscanf(SEQ,"%s\n", seq_tgDNA);
	fclose(SEQ);
	for (int i=0; i<(int)strlen(seq_tgDNA); i++)
	{
		if(seq_tgDNA[i]=='A')
			seq_ntDNA[i]='T';
		else if(seq_tgDNA[i]=='C')
			seq_ntDNA[i]='G';
		else if(seq_tgDNA[i]=='G')
			seq_ntDNA[i]='C';
		else if(seq_tgDNA[i]=='T')
			seq_ntDNA[i]='A';
		else
		{
			printf("Sequence error\n");
			exit(-1);
		}
	}
	if(0)
	{
		printf("Reading sgRNA Sequence from sgRNA.fasta");
		SEQ=fopen("sgRNA.fasta","r+");
		fscanf(SEQ,"%s\n", seq_sgRNA);
		fclose(SEQ);
	}
	else if(1)
	{
		for (int i=0; i<(int)strlen(seq_tgDNA); i++)
		{
			if(seq_tgDNA[i]=='A')
				seq_sgRNA[i]='U';
			else if(seq_tgDNA[i]=='C')
				seq_sgRNA[i]='G';
			else if(seq_tgDNA[i]=='G')
				seq_sgRNA[i]='C';
			else if(seq_tgDNA[i]=='T')
				seq_sgRNA[i]='A';
			else
			{
				printf("Sequence error\n");
				exit(-1);
			}
		}
	}
	//printf("%s\n%s\n%s\n",seq_tgDNA,seq_ntDNA,seq_sgRNA);
}

float calculate_rate(float dGij)
{
	float kij;
	
	kij=0.;
	
	kij=K0*exp((0.-dGij)/kbT);
	
	if(K0<kij)
		return K0;
	else
		return kij;
}

void KMC_initialation(int KMC_trail)
{
	path_all[KMC_trail]=(path_node *) malloc(sizeof(path_node));
	path_all[KMC_trail]->structure_node_C=(structure_node *) malloc(sizeof(structure_node));
	path_all[KMC_trail]->node_before=NULL;
	path_all[KMC_trail]->node_next=(path_node *) malloc(sizeof(path_node));
	path_all[KMC_trail]->last=path_all[KMC_trail];
	path_all[KMC_trail]->structure_node_C->step_mc=0;
	memset(path_all[KMC_trail]->structure_node_C->sequence_ntDNA, 0, sizeof(path_all[KMC_trail]->structure_node_C->sequence_ntDNA));
	memset(path_all[KMC_trail]->structure_node_C->sequence_sgRNA, 0, sizeof(path_all[KMC_trail]->structure_node_C->sequence_sgRNA));
	sprintf(path_all[KMC_trail]->structure_node_C->sequence_ntDNA, "%s", seq_ntDNA);
	sprintf(path_all[KMC_trail]->structure_node_C->sequence_sgRNA, "%s", seq_sgRNA);
	memset(path_all[KMC_trail]->structure_node_C->base_pairs, 0, sizeof(path_all[KMC_trail]->structure_node_C->base_pairs));
	path_all[KMC_trail]->structure_node_C->energy_dsDNA=0.;
	path_all[KMC_trail]->structure_node_C->energy_RDNA=0.;
	path_all[KMC_trail]->structure_node_C->energy_ssRNA=0.;
	path_all[KMC_trail]->structure_node_C->energy_ssDNA=0.;
	memset(path_all[KMC_trail]->structure_node_C->neighbor, 0, sizeof(path_all[KMC_trail]->structure_node_C->neighbor));

}

void KMC_simulation()
{
	for(int KMC_trail=0; KMC_trail<KMC_MAXT; KMC_trail++)
	{
		KMC_initialation(KMC_trail);
		for(int run_step=0; run_step<RUN_MAX; run_step++)
		{
			//NODE_search();
			//KMC_run();
		}
		
	}
}


void NODE_search()
{
	for(int i=1, j=0;i<strlen(path_all[KMC_trail]->structure_node_C->sequence_sgRNA)+1;i++)
	{
		if(path_all[KMC_trail]->structure_node_C->base_pairs[i]>0) continue; //paired
		if(path_all[KMC_trail]->structure_node_C->base_pairs[i]<0) 
		{
			j=j+path_all[KMC_trail]->structure_node_C->base_pairs[i];
			continue;
		} //Bulge unpaired
		if(path_all[KMC_trail]->structure_node_C->base_pairs[i]==0) //End unpaired
		{
		
		
			path_all[KMC_trail]->structure_node_C->neighbor[1]=(structure_node *) malloc(sizeof(structure_node));
			path_all[KMC_trail]->structure_node_C->neighbor[1]->step_mc++;
			strcpy(path_all[KMC_trail]->structure_node_C->neighbor[1]->sequence_ntDNA,path_all[KMC_trail]->structure_node_C->sequence_ntDNA);
			strcpy(path_all[KMC_trail]->structure_node_C->neighbor[1]->sequence_sgRNA,path_all[KMC_trail]->structure_node_C->sequence_sgRNA);
			memset(path_all[KMC_trail]->structure_node_C->neighbor[1]->neighbor, 0, sizeof(path_all[KMC_trail]->structure_node_C->neighbor));
			for(int k=0;k<50;k++)	path_all[KMC_trail]->structure_node_C->neighbor[1]->base_pairs[k]=path_all[KMC_trail]->structure_node_C->base_pairs[k];
			path_all[KMC_trail]->structure_node_C->neighbor[1]->base_pairs[i]=i+j;
			calculate_energy(structure_node path_all[KMC_trail]->structure_node_C->neighbor[1]);
			
			if(i==0) break;
			path_all[KMC_trail]->structure_node_C->neighbor[2]=(structure_node *) malloc(sizeof(structure_node));
			path_all[KMC_trail]->structure_node_C->neighbor[2]->step_mc++;
			strcpy(path_all[KMC_trail]->structure_node_C->neighbor[2]->sequence_ntDNA,path_all[KMC_trail]->structure_node_C->sequence_ntDNA);
			strcpy(path_all[KMC_trail]->structure_node_C->neighbor[2]->sequence_sgRNA,path_all[KMC_trail]->structure_node_C->sequence_sgRNA);
			memset(path_all[KMC_trail]->structure_node_C->neighbor[2]->neighbor, 0, sizeof(path_all[KMC_trail]->structure_node_C->neighbor));
			for(int k=0;k<50;k++)	path_all[KMC_trail]->structure_node_C->neighbor[1]->base_pairs[k]=path_all[KMC_trail]->structure_node_C->base_pairs[k];
			path_all[KMC_trail]->structure_node_C->neighbor[2]->base_pairs[i-1]=0;
			calculate_energy(structure_node path_all[KMC_trail]->structure_node_C->neighbor[1]);
			
			break;
		}
	}
}


void KMC_run()
{
	float k_p[30];
	memset(k_p,0.,sizeof(k_p));
	float dG=0.;
	int Nneighbor;
	for(int i=0; i<29; i++)
	{
		if (path_all[KMC_trail]->structure_node_C->neighbor[i]==NULL)
		{
			Nneighbor=i;
			break;
		}
		dG=path_all[KMC_trail]->structure_node_C->neighbor[i+1]->energy_RDNA-path_all[KMC_trail]->structure_node_C->energy_RDNA;
		k_p[i+1]=calculate_rate(dG);
	}
	for(int i=0; i<Nneighbor; i++)
	{
		
	}

}

void calculate_energy(structure_node *structure_1)
{

	float En0=0.;
	//PAM region
	
	
	
	for(int i=0; i<strlen(structure_1->sequence_sgRNA); i++)
	{
		if (i==0)
		{
			if(structure_1->base_pairs[i+1]==0)
			{
				En0+=0.;
				break;			
			}
			else if(structure_1->base_pairs[i+1]==1)
			{
				En0+=0.;
				continue;			
			}
		}
		else if (structure_1->base_pairs[i+1]==1)
		{
			//En0+=Paraset_on_stenergy[][]
			if (structure_1->sequence_ntDNA[i]=='A')
			{
				if (structure_1->sequence_ntDNA[i+1]=='A')
					En0+=Paraset_on_stenergy[i][0];
				else if(structure_1->sequence_ntDNA[i+1]=='C')
					En0+=Paraset_on_stenergy[i][1];
				else if(structure_1->sequence_ntDNA[i+1]=='G')
					En0+=Paraset_on_stenergy[i][2];
				else if(structure_1->sequence_ntDNA[i+1]=='T')
					En0+=Paraset_on_stenergy[i][3];
			}
			else if (structure_1->sequence_ntDNA[i]=='C')
			{
				if (structure_1->sequence_ntDNA[i+1]=='A')
					En0+=Paraset_on_stenergy[i][4];
				else if(structure_1->sequence_ntDNA[i+1]=='C')
					En0+=Paraset_on_stenergy[i][5];
				else if(structure_1->sequence_ntDNA[i+1]=='G')
					En0+=Paraset_on_stenergy[i][6];
				else if(structure_1->sequence_ntDNA[i+1]=='T')
					En0+=Paraset_on_stenergy[i][7];
			}
			if (structure_1->sequence_ntDNA[i]=='G')
			{
				if (structure_1->sequence_ntDNA[i+1]=='A')
					En0+=Paraset_on_stenergy[i][8];
				else if(structure_1->sequence_ntDNA[i+1]=='C')
					En0+=Paraset_on_stenergy[i][9];
				else if(structure_1->sequence_ntDNA[i+1]=='G')
					En0+=Paraset_on_stenergy[i][10];
				else if(structure_1->sequence_ntDNA[i+1]=='T')
					En0+=Paraset_on_stenergy[i][11];
			}
			if (structure_1->sequence_ntDNA[i]=='T')
			{
				if (structure_1->sequence_ntDNA[i+1]=='A')
					En0+=Paraset_on_stenergy[i][12];
				else if(structure_1->sequence_ntDNA[i+1]=='C')
					En0+=Paraset_on_stenergy[i][13];
				else if(structure_1->sequence_ntDNA[i+1]=='G')
					En0+=Paraset_on_stenergy[i][14];
				else if(structure_1->sequence_ntDNA[i+1]=='T')
					En0+=Paraset_on_stenergy[i][15];
			}
		
		}
		else if (structure_1->base_pairs[i+1]==0)
		{
			break;
		}
	
	}
	structure_1->energy_RDNA=En0;

}









