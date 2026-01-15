#include<stdio.h>  
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

typedef struct
{
	int position;
	char stacking_type[3];
	float energy;
}	Paraset;


typedef struct
{
	int step_mc;
	char sequence_ntDNA[50];
	char sequence_sgRNA[50];
	int base_pairs[50];
	float energy_dsDNA;
	float energy_RDNA;
	float energy_ssRNA;
	float energy_ssDNA;
	structure_node *neighbor[30]
}	structure_node;

typedef struct
{
	structure_node *structure_node_C;
	path_node *node_before;
	path_node *node_next;
	path_node *last;
}	path_node;

extern int sttype();

Paraset Paraset_on_stenergy[30][16];
Paraset Paraset_off_stenergy[30][16];
extern void read_para();

char seq_tgDNA[30], seq_ntDNA[30], seq_sgRNA[30];
extern void set_seq();

path_node **path_all;
extern void KMC_initialation();
extern void KMC_simulation();
extern void Data_output();

extern float calculate_rate(structure_node structure_1);


int main()
{
	FILE *Pop_step, *Pop_fin, *Path_step, *Path_avg;
	
	//Initialization//	
	read_para();
	set_seq();
	
	KMC_simulation();
	
	data_output();
	
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

float calculate_rate()
{
	return 0;
}










