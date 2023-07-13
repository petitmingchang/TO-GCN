#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

int num_of_genes;
int num_of_TFs;
int num_of_seeds;

int num_of_point_LD;
int num_of_point_TD;

typedef struct Time_Course{
        char gene_ID[20];
        double *LD_exp;
        double *TD_exp;
        int level;
        int network[10];
}Time_Course;

Time_Course *gene_exp_table;
Time_Course *TF_exp_table;
Time_Course *gene_seed_table;

typedef struct Relation_table{
    char query_gene_ID[20];
    char target_gene_ID[20];
    int checked;
} R_table;

Relation_table *Pos_Coexp;

int num_of_pos_edge = 0;
int num_of_neg_edge = 0;

float pos_cutoff_LD;
float pos_cutoff_TD;

float neg_cutoff_LD;
float neg_cutoff_TD;

float pos_no_cutoff = 0.5;
float neg_no_cutoff = -0.5;

char seed_TF_ID[20];  //seed_TF_ID for testing: "Zm00001d041056" (ZmARF2) for C1_+C2+

int done = 0;

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
    	TF_exp_table[i].level = -1;
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
    	gene_exp_table[i].level = -1;
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

void Read_seed_genes (char *input) {
    FILE *fptr = fopen(input, "r");
    char GID[20];
    
    num_of_seeds = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
        num_of_seeds++;
    }

    rewind(fptr);
    //table initialization
    gene_seed_table = new Time_Course[num_of_seeds];

    int index = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
        strcpy(gene_seed_table[index].gene_ID, GID);
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

void node_pair_generator_LD_or_TD(int opt) {

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

    if(opt == 0) { //considering both positive and negative correlations  C1+C2+
        for(int i=0; i<num_of_TFs; i++) {
            for(int j=0; j<num_of_genes; j++) {
                if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                    R_LD = r_calculator(i,j,0);
                    R_TD = r_calculator(i,j,1);
    		    
                    if(R_LD >= pos_cutoff_LD && R_TD >= pos_cutoff_TD) { // C1+C2+
                        num_of_pos_edge++;
                    }
                }
            }
        }
    } else if(opt == 1) { //considering positive correlation in LD only  C1+C20
        for(int i=0; i<num_of_TFs; i++) {
            for(int j=0; j<num_of_genes; j++) {
                if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                    R_LD = r_calculator(i,j,0);
                    R_TD = r_calculator(i,j,1);
   
                    if(R_LD >= pos_cutoff_LD && R_TD < pos_no_cutoff && R_TD >= neg_no_cutoff) { // C1+C20
                        num_of_pos_edge++;
                    }
                }
            }
        }
    } else if(opt == 2) { //considering positive correlation in TD only C10C2+
        for(int i=0; i<num_of_TFs; i++) {
            for(int j=0; j<num_of_genes; j++) {
                if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                    R_LD = r_calculator(i,j,0);
                    R_TD = r_calculator(i,j,1);
                    
                    if(R_TD >= pos_cutoff_TD && R_LD < pos_no_cutoff && R_LD >= neg_no_cutoff) { // C10C2+
                        num_of_pos_edge++;
                    }
                }
            }
        }
    } else if(opt == 3) { //considering positive correlation in LD only  C1+C2-
        for(int i=0; i<num_of_TFs; i++) {
            for(int j=0; j<num_of_genes; j++) {
                if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                    R_LD = r_calculator(i,j,0);
                    R_TD = r_calculator(i,j,1);
   
                    if(R_LD >= pos_cutoff_LD && R_TD < neg_cutoff_TD) { // C1+C2-
                        num_of_pos_edge++;
                    }
                }
            }
        }
    } else if(opt == 4) { //considering positive correlation in TD only C1-C2+
        for(int i=0; i<num_of_TFs; i++) {
            for(int j=0; j<num_of_genes; j++) {
                if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                    R_LD = r_calculator(i,j,0);
                    R_TD = r_calculator(i,j,1);
                    
                    if(R_TD >= pos_cutoff_TD && R_LD < neg_cutoff_LD) { // C1-C2+
                        num_of_pos_edge++;
                    }
                }
            }
        }
	}    
    
     

    Pos_Coexp = new R_table[num_of_pos_edge];
    
    for (int i=0; i<num_of_pos_edge; i++) {
        Pos_Coexp[i].checked = 0;
    }
    
    int index_pos = 0;
    int index_neg = 0;
    for(int i=0; i<num_of_TFs; i++) {
        for(int j=0; j<num_of_genes; j++) {
            if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                R_LD = r_calculator(i,j,0);
                R_TD = r_calculator(i,j,1);

                if(opt == 0) { //C1+C2+
                    if(R_LD >= pos_cutoff_LD && R_TD >= pos_cutoff_TD) {
                        strcpy(Pos_Coexp[index_pos].query_gene_ID, TF_exp_table[i].gene_ID);
                        strcpy(Pos_Coexp[index_pos].target_gene_ID, gene_exp_table[j].gene_ID);
                        index_pos++;
                    }
                } else if(opt == 1) { //C1+C20
                    if(R_LD >= pos_cutoff_LD && R_TD < pos_no_cutoff && R_TD >= neg_no_cutoff) {
                        strcpy(Pos_Coexp[index_pos].query_gene_ID, TF_exp_table[i].gene_ID);
                        strcpy(Pos_Coexp[index_pos].target_gene_ID, gene_exp_table[j].gene_ID);
                        index_pos++;
                    }
                } else if(opt == 2) { //C10C2+
                    if(R_TD >= pos_cutoff_TD && R_LD < pos_no_cutoff && R_LD >= neg_no_cutoff) {
                        strcpy(Pos_Coexp[index_pos].query_gene_ID, TF_exp_table[i].gene_ID);
                        strcpy(Pos_Coexp[index_pos].target_gene_ID, gene_exp_table[j].gene_ID);
                        index_pos++;
                    }
                } else if(opt == 3) { //C1+C2-
                    if(R_LD >= pos_cutoff_LD && R_TD < neg_cutoff_TD) {
                        strcpy(Pos_Coexp[index_pos].query_gene_ID, TF_exp_table[i].gene_ID);
                        strcpy(Pos_Coexp[index_pos].target_gene_ID, gene_exp_table[j].gene_ID);
                        index_pos++;
                    }
                } else if(opt == 4) { //C1-C2+
                    if(R_TD >= pos_cutoff_TD && R_LD < neg_cutoff_LD) {
                        strcpy(Pos_Coexp[index_pos].query_gene_ID, TF_exp_table[i].gene_ID);
                        strcpy(Pos_Coexp[index_pos].target_gene_ID, gene_exp_table[j].gene_ID);
                        index_pos++;
                    }
                }
            }
        }
    }

}

void set_neighbor_pos(char *check_ID, int L) {
    
    for (int i=0; i<num_of_pos_edge; i++) {
        if (Pos_Coexp[i].checked == 0) {
            if (strcmp(check_ID, Pos_Coexp[i].query_gene_ID) == 0) {
                Pos_Coexp[i].checked = 1;
                for (int j=0; j<num_of_TFs; j++) {
                    if(strcmp(Pos_Coexp[i].target_gene_ID, TF_exp_table[j].gene_ID) == 0) {
                        if (TF_exp_table[j].level == -1) {
                            TF_exp_table[j].level = L;
                            done = 0;
                        }
                    }
                }
            }
        }
    }
    
    for (int i=0; i<num_of_pos_edge; i++) {
        if (Pos_Coexp[i].checked == 0) {
            if (strcmp(check_ID, Pos_Coexp[i].target_gene_ID) == 0) {
                Pos_Coexp[i].checked = 1;
                for (int j=0; j<num_of_TFs; j++) {
                    if(strcmp(Pos_Coexp[i].query_gene_ID, TF_exp_table[j].gene_ID) == 0) {
                        if (TF_exp_table[j].level == -1) {
                            TF_exp_table[j].level = L;
                            done = 0;
                        }
                    }
                }
            }
        }
    }
}

void level_assignment() {
    
    int level = 1;
    int num_of_found_seeds = 0;
    int seed_found_in_GCN = 0;

    for (int i=0; i<num_of_TFs; i++) {
        for (int j=0; j<num_of_seeds; j++) {
            if(strcmp(gene_seed_table[j].gene_ID, TF_exp_table[i].gene_ID) == 0) {
                    TF_exp_table[i].level = level;
                num_of_found_seeds++;
            }
        }
    }
    
    if(num_of_found_seeds == num_of_seeds) {
        seed_found_in_GCN = 1;
    }
    
    if (seed_found_in_GCN == 0) {
        printf("\nSome ID cannot be found in the GCN. Please check the see ID list again!\n\n");
    } else {

        for (int i=0; i<num_of_TFs; i++) {
            for (int j=0; j<num_of_seeds; j++) {
                if(strcmp(gene_seed_table[j].gene_ID, TF_exp_table[i].gene_ID) == 0) {
                    set_neighbor_pos(TF_exp_table[i].gene_ID, level+1);
                }
            }
        }
        
        while (done == 0) {
            level++;
            done = 1;
            for (int i=0; i<num_of_TFs; i++) {
                if(TF_exp_table[i].level == level) {
                    set_neighbor_pos(TF_exp_table[i].gene_ID, level+1);
                }
            }
        }
    }
}

void function_three () {
    FILE *fout3;
    
    fout3 = fopen("TF_level.csv","w");
    
    fprintf(fout3, "TF gene ID,level in GCN\n");
    for(int i=0; i<num_of_TFs; i++) {
        if(TF_exp_table[i].level >= 0) {
            fprintf(fout3, "%s, %d\n", TF_exp_table[i].gene_ID, TF_exp_table[i].level);
        }
    }
    
    fclose(fout3);
}
 
int main(int argc, char* argv[]) {
    char input_file1[100];   //TF gene list
    char input_file2[100];   //All gene list
    char input_file3[100];   //seed gene list
    int coex_type; //coexpression type
    
    if (argc != 11) {
        printf("\nUsage: TO-GCN #Cond1_samples #Cond2_samples file_of_TF_genes file_of_all_genes Cutoff_pos_C1 Cutoff_pos_C2 Cutoff_neg_C1 Cutoff_neg_C2 Seed_ID_list coexpression_type\n");
        printf("coexpression_type = 0: C1+C2+\n");
        printf("coexpression_type = 1: C1+C20\n");
        printf("coexpression_type = 2: C10+C2+\n");
        printf("coexpression_type = 3: C1+C2-\n");
        printf("coexpression_type = 4: C1-C2+\n\n");
        
    } else {

        num_of_point_LD = atoi(argv[1]);
        num_of_point_TD = atoi(argv[2]);
        
        strcpy(input_file1, argv[3]);
        strcpy(input_file2, argv[4]);
        
        pos_cutoff_LD = atof(argv[5]);
        pos_cutoff_TD = atof(argv[6]);
        
        neg_cutoff_LD = atof(argv[7]);
        neg_cutoff_TD = atof(argv[8]);      
        
        strcpy(input_file3, argv[9]);
        coex_type = atoi(argv[10]);

        FILE *fptr1 = fopen(input_file1, "r");
        FILE *fptr2 = fopen(input_file2, "r");
        FILE *fptr3 = fopen(input_file3, "r");
        
        if(fptr1 == NULL || fptr2 == NULL || fptr3 == NULL) {
            
            printf("\nCan't find the input file. Please check the inupt file again!\n\n");
            
        } else {
            
            fclose (fptr1);
            fclose (fptr2);
            fclose (fptr3);
        
            Read_Time_Course_Data_TFs(input_file1);
            Read_Time_Course_Data_genes(input_file2);
            Read_seed_genes(input_file3);

            printf("NO. of TFs: %d\n", num_of_TFs);
            printf("NO. of Genes: %d\n", num_of_genes);
            printf("No. of samples under Cond.1: %d\n", num_of_point_LD);
            printf("No. of samples under Cond.2: %d\n", num_of_point_TD);
            printf("Cutoffs for (Pos_C1, Pos_C2, Neg_C1, Neg_C2): (%1.2lf, %1.2lf, %1.2lf, %1.2lf)\n\n", pos_cutoff_LD, pos_cutoff_TD, neg_cutoff_LD, neg_cutoff_TD);
            printf("Assigning levels for TFs in GCN by Breadth-First-Search (BFS) method......\n");
        
            node_pair_generator_LD_or_TD(coex_type); //0: for C1+C2+; 1: for C1+C20; 2: for C10C2+ 3: for C1+C2-; 4: for C1-C2+
            level_assignment();
            function_three();
        
            printf("Done!\n");
        }
    }
    
    return 0;	
}
