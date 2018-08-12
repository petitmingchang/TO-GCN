#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int num_of_genes;
int num_of_TFs;

int num_of_point_LD;
int num_of_point_TD;

typedef struct Time_Course{
        char gene_ID[20];
        double *LD_exp;
        double *TD_exp;
        int network[10];
}Time_Course;

Time_Course *gene_exp_table;
Time_Course *TF_exp_table;

float pos_cutoff_LD;
float pos_cutoff_TD;
float neg_cutoff_LD;
float neg_cutoff_TD;
float pos_no_cutoff = 0.5;
float neg_no_cutoff = -0.5;

void Read_Time_Course_Data_TFs (char *input) {
	FILE *fptr = fopen(input, "r");
	char GID[20];
	double LDE[num_of_point_LD];
	double TDE[num_of_point_TD];

	num_of_TFs = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
		}
		for(int j=0; j<num_of_point_TD; j++) {
			fscanf(fptr,"\t%lf", &TDE[j]);
		}
		num_of_TFs++;
	}		

	rewind(fptr);
	//table initialization
	
	TF_exp_table = new Time_Course[num_of_TFs];
    
    for(int i=0; i<num_of_TFs; i++) {
        TF_exp_table[i].LD_exp = new double[num_of_point_LD];
        TF_exp_table[i].TD_exp = new double[num_of_point_TD];
    }
    
	for(int i=0; i<num_of_TFs; i++) {

    	for(int k=0; k<10; k++) {
    		TF_exp_table[i].network[k] = 0;
    	}
        for(int j=0; j<num_of_point_LD; j++) {
            TF_exp_table[i].LD_exp[j] = 0;
        }
        for(int j=0; j<num_of_point_TD; j++) {
            TF_exp_table[i].TD_exp[j] = 0;
        }
    }
    

	int index = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		strcpy(TF_exp_table[index].gene_ID, GID);
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
			TF_exp_table[index].LD_exp[i] = LDE[i];
		}
		for(int j=0; j<num_of_point_TD; j++) {
			fscanf(fptr,"\t%lf", &TDE[j]);
			TF_exp_table[index].TD_exp[j] = TDE[j];
		}
		index++;
	}		
	
    fclose(fptr);    
}

void Read_Time_Course_Data_genes (char *input) {
	FILE *fptr = fopen(input, "r");
	char GID[20];
	double LDE[num_of_point_LD];
	double TDE[num_of_point_TD];
	
	num_of_genes = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
		}
		for(int j=0; j<num_of_point_TD; j++) {
			fscanf(fptr,"\t%lf", &TDE[j]);
		}
		num_of_genes++;
	}	

	rewind(fptr);
	//table initialization
	gene_exp_table = new Time_Course[num_of_genes];
    
    for(int i=0; i<num_of_genes; i++) {
        gene_exp_table[i].LD_exp = new double[num_of_point_LD];
        gene_exp_table[i].TD_exp = new double[num_of_point_TD];
    }
    
	for(int i=0; i<num_of_genes; i++) {

    	for(int k=0; k<10; k++) {
    		gene_exp_table[i].network[k] = 0;
    	}
        for(int j=0; j<num_of_point_LD; j++) {
            gene_exp_table[i].LD_exp[j] = 0;
        }
        for(int j=0; j<num_of_point_TD; j++) {
            gene_exp_table[i].TD_exp[j] = 0;
        }
    }	

	int index = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		strcpy(gene_exp_table[index].gene_ID, GID);
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
			gene_exp_table[index].LD_exp[i] = LDE[i];
		}
		for(int j=0; j<num_of_point_TD; j++) {
			fscanf(fptr,"\t%lf", &TDE[j]);
			gene_exp_table[index].TD_exp[j] = TDE[j];
		}
		index++;
	}		
	
    fclose(fptr);
}

double r_calculator (int x, int y, int z) {
	double N_LD = num_of_point_LD;
	double N_TD = num_of_point_TD;
	double R, SUM_XY, SUM_X, SUM_Y, SUM_X2, SUM_Y2;
	
	double temp_XY, temp_X, temp_Y, temp_X2, temp_Y2;
	temp_XY = temp_X = temp_Y = temp_X2 = temp_Y2 = 0;
	if(z == 0) { //calculate r of LD expression data
		for(int i=0; i<N_LD; i++) {
			temp_XY = temp_XY + (TF_exp_table[x].LD_exp[i] * gene_exp_table[y].LD_exp[i]);
			temp_X = temp_X + TF_exp_table[x].LD_exp[i];
			temp_Y = temp_Y + gene_exp_table[y].LD_exp[i];
			temp_X2 = temp_X2 + (TF_exp_table[x].LD_exp[i] * TF_exp_table[x].LD_exp[i]);
			temp_Y2 = temp_Y2 + (gene_exp_table[y].LD_exp[i] * gene_exp_table[y].LD_exp[i]);
		}
	} else { //calculate r of TD expression data
		for(int j=0; j<N_TD; j++) {
			temp_XY = temp_XY + (TF_exp_table[x].TD_exp[j] * gene_exp_table[y].TD_exp[j]);
			temp_X = temp_X + TF_exp_table[x].TD_exp[j];
			temp_Y = temp_Y + gene_exp_table[y].TD_exp[j];
			temp_X2 = temp_X2 + (TF_exp_table[x].TD_exp[j] * TF_exp_table[x].TD_exp[j]);
			temp_Y2 = temp_Y2 + (gene_exp_table[y].TD_exp[j] * gene_exp_table[y].TD_exp[j]);
		}	
	}
	
	SUM_XY = temp_XY;
	SUM_X = temp_X;
	SUM_Y = temp_Y;
	SUM_X2 = temp_X2;
	SUM_Y2 = temp_Y2;
    
    if (SUM_X == 0 || SUM_Y == 0) {
        R = 0;
    } else {
        if(z == 0) { //calculate r of LD expression data
            R = (SUM_XY-(SUM_X)*(SUM_Y)/N_LD)/(sqrt((SUM_X2-(SUM_X)*(SUM_X)/N_LD)*(SUM_Y2-(SUM_Y)*(SUM_Y)/N_LD)));
        } else { //calculate r of TD expression data
            R = (SUM_XY-(SUM_X)*(SUM_Y)/N_TD)/(sqrt((SUM_X2-(SUM_X)*(SUM_X)/N_TD)*(SUM_Y2-(SUM_Y)*(SUM_Y)/N_TD)));
        }
    }
	return R;
}

void node_pair_generator_LD_or_TD() {
    FILE *fout1, *fout2, *fout3, *fout4, *fout5, *fout6, *fout7, *fout8;
    
    double R_LD;
    double R_TD;

    for(int i=0; i<num_of_TFs; i++) {
        for(int j=0; j<num_of_genes; j++) {
			for(int k=0; k<10; k++) {
				TF_exp_table[i].network[k] = 0;
				gene_exp_table[j].network[k] = 0;
			}
		}
	}
    
    fout1 = fopen("C1+C2+.csv","w");
    fout2 = fopen("C1+C20.csv","w");
    fout3 = fopen("C10C2+.csv","w");
    fout4 = fopen("C1-C2-.csv","w");
    fout5 = fopen("C2-C20.csv","w");
    fout6 = fopen("C10C2-.csv","w");
    fout7 = fopen("C1+C2-.csv","w");
    fout8 = fopen("C1-C2+.csv","w");

    fprintf(fout1, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    fprintf(fout2, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    fprintf(fout3, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    fprintf(fout4, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    fprintf(fout5, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    fprintf(fout6, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    fprintf(fout7, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    fprintf(fout8, "TF gene ID, coexpression type, gene ID, PCC under C1, PCC under C2\n");
    
        for(int i=0; i<num_of_TFs; i++) {
            for(int j=0; j<num_of_genes; j++) {
                if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                    R_LD = r_calculator(i,j,0);
                    R_TD = r_calculator(i,j,1);
                    
                    if(R_LD >= pos_cutoff_LD && R_TD >= pos_cutoff_TD) { // (3) LD+TD+
                        fprintf(fout1, "%s,C1+C2+,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[3]++;
                        gene_exp_table[j].network[3]++;
                    } else if(R_LD >= pos_cutoff_LD && R_TD < pos_no_cutoff && R_TD >= neg_no_cutoff) { // (1)LD+TD0
                        fprintf(fout2, "%s,C1+C20,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[1]++;
                        gene_exp_table[j].network[1]++;
                    } else if(R_TD >= pos_cutoff_TD && R_LD < pos_no_cutoff && R_LD >= neg_no_cutoff) { // (5)LD0TD+
                        fprintf(fout3, "%s,C10C2+,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[5]++;
                        gene_exp_table[j].network[5]++;
                    } else if(R_LD < neg_cutoff_LD && R_TD < neg_cutoff_TD) { // (4)LD-TD-
                        fprintf(fout4, "%s,C1-C2-,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[4]++;
                        gene_exp_table[j].network[4]++;
                    } else if(R_LD < neg_cutoff_LD && R_TD < pos_no_cutoff && R_TD >= neg_no_cutoff) { // (2)LD-TD0
                        fprintf(fout5, "%s,C1-C20,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[2]++;
                        gene_exp_table[j].network[2]++;
                    } else if(R_TD < neg_cutoff_TD && R_LD < pos_no_cutoff && R_LD >= neg_no_cutoff) { // (6)LD0TD-
                        fprintf(fout6, "%s,C10C2-,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[6]++;
                        gene_exp_table[j].network[6]++;
                    } else if(R_LD >= pos_cutoff_LD && R_TD < neg_cutoff_TD) { // (7)LD+TD-
                        fprintf(fout7, "%s,C1+C2-,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[7]++;
                        gene_exp_table[j].network[7]++;
                    } else if(R_TD >= pos_cutoff_TD && R_LD < neg_cutoff_LD) { // (8)LD-TD+
                        fprintf(fout8, "%s,C1-C2+,%s,%lf,%lf\n", TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID, R_LD, R_TD);
                        TF_exp_table[i].network[8]++;
                        gene_exp_table[j].network[8]++;
                    }
                }
            }
        }
    
    fclose(fout1);
    fclose(fout2);
    fclose(fout3);    
    fclose(fout4);
    fclose(fout5);
    fclose(fout6);
    fclose(fout7);
    fclose(fout8);
}


void function_two () {
    FILE *fout1, *fout2;
    fout1 = fopen("TF_edges_in_eight_coexpression_types.csv","w");
    fout2 = fopen("Gene_edges_in_eight_coexpression_types.csv","w");

    fprintf(fout1, "Gene ID,C1+C20,C1-C20,C1+C2+,C1-C2-,C10C2+,C10C2-,C1+C2-,C1-C2+\n");
    fprintf(fout2, "Gene ID,C1+C20,C1-C20,C1+C2+,C1-C2-,C10C2+,C10C2-,C1+C2-,C1-C2+\n");
    
    for(int i=0; i<num_of_TFs; i++) {
        fprintf(fout1, "%s", TF_exp_table[i].gene_ID);
        for(int k=1; k<9; k++) {
            fprintf(fout1, ",%d", TF_exp_table[i].network[k]);
        }
        fprintf(fout1, "\n");
    }
    
    
    for(int i=0; i<num_of_genes; i++) {
        fprintf(fout2, "%s", gene_exp_table[i].gene_ID);
        for(int k=1; k<9; k++) {
            fprintf(fout2, ",%d", gene_exp_table[i].network[k]);
        }
        fprintf(fout2, "\n");
    }
    
    fclose(fout1);
    fclose(fout2);
}

int main(int argc, char* argv[]) {
    char input_file1[100];   //TF gene list
    char input_file2[100];   //All gene list
    
    if (argc != 9) {
        printf("\nUsage: GCN #Cond1_samples #Cond2_samples file_of_TF_genes file_of_all_genes Cutoff_pos_C1 Cutoff_pos_C2 Cutoff_neg_C1 Cutoff_neg_C2\n\n");
        
    } else {
        
        num_of_point_LD = atoi(argv[1]);
        num_of_point_TD = atoi(argv[2]);
        
        strcpy(input_file1, argv[3]);
        strcpy(input_file2, argv[4]);
        
        pos_cutoff_LD = atof(argv[5]);
        pos_cutoff_TD = atof(argv[6]);
        neg_cutoff_LD = atof(argv[7]);
        neg_cutoff_TD = atof(argv[8]);

        FILE *fptr1 = fopen(input_file1, "r");
        FILE *fptr2 = fopen(input_file2, "r");
        
        if(fptr1 == NULL || fptr2 == NULL) {
            
            printf("\nCan't find the input file. Please check the inupt file again!\n\n");
            
        } else {
        
            Read_Time_Course_Data_TFs(input_file1);
            Read_Time_Course_Data_genes(input_file2);

            printf("NO. of TFs: %d\n", num_of_TFs);
            printf("NO. of Genes: %d\n", num_of_genes);
            printf("No. of samples under Cond. 1: %d\n", num_of_point_LD);
            printf("No. of samples under Cond. 2: %d\n", num_of_point_TD);
            printf("Cutoff values for (Pos_C1, Pos_C2, Neg_C1, Neg_C2): (%1.2lf, %1.2lf, %1.2lf, %1.2lf)\n\n", pos_cutoff_LD, pos_cutoff_TD, neg_cutoff_LD, neg_cutoff_TD);
            printf("Generating tables of eight types of coexpression......\n");
        
            node_pair_generator_LD_or_TD();
        
            printf("Done!\n");
        }
    }
    
    return 0;
}
