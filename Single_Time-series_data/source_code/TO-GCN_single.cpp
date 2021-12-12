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
    double pcc;
    int checked;
} R_table;

Relation_table *Pos_Coexp;

int num_of_pos_edge = 0;
//int num_of_neg_edge = 0;

float pos_cutoff_LD;

int done = 0;

void Read_Time_Course_Data_TFs (char *input) {
	FILE *fptr = fopen(input, "r");
	char GID[20];
	double LDE[num_of_point_LD];

	num_of_TFs = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
		}

		num_of_TFs++;
	}		

	rewind(fptr);
	//table initialization
	
	TF_exp_table = new Time_Course[num_of_TFs];
    
    for(int i=0; i<num_of_TFs; i++) {
        TF_exp_table[i].LD_exp = new double[num_of_point_LD];
    }
    
	for(int i=0; i<num_of_TFs; i++) {
    	TF_exp_table[i].level = -1;
    	for(int k=0; k<10; k++) {
    		TF_exp_table[i].network[k] = 0;
    	}
        for(int j=0; j<num_of_point_LD; j++) {
            TF_exp_table[i].LD_exp[j] = 0;
        }
    }
    

	int index = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		strcpy(TF_exp_table[index].gene_ID, GID);
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
			TF_exp_table[index].LD_exp[i] = LDE[i];
		}

		index++;
	}		
	
    fclose(fptr);    
}

void Read_Time_Course_Data_genes (char *input) {
	FILE *fptr = fopen(input, "r");
	char GID[20];
	double LDE[num_of_point_LD];
	
	num_of_genes = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
		}

		num_of_genes++;
	}	

	rewind(fptr);
	//table initialization
	gene_exp_table = new Time_Course[num_of_genes];
    
    for(int i=0; i<num_of_genes; i++) {
        gene_exp_table[i].LD_exp = new double[num_of_point_LD];
    }
    
	for(int i=0; i<num_of_genes; i++) {
    	gene_exp_table[i].level = -1;
    	for(int k=0; k<10; k++) {
    		gene_exp_table[i].network[k] = 0;
    	}
        for(int j=0; j<num_of_point_LD; j++) {
            gene_exp_table[i].LD_exp[j] = 0;
        }
    }	

	int index = 0;
    while(fscanf(fptr,"%s", GID) != EOF) {
		strcpy(gene_exp_table[index].gene_ID, GID);
		for(int i=0; i<num_of_point_LD; i++) {
			fscanf(fptr,"\t%lf", &LDE[i]);
			gene_exp_table[index].LD_exp[i] = LDE[i];
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

double r_calculator (int x, int y) {
    double N_LD = num_of_point_LD;
    double R, SUM_XY, SUM_X, SUM_Y, SUM_X2, SUM_Y2;
    double temp_XY, temp_X, temp_Y, temp_X2, temp_Y2;
    
    temp_XY = temp_X = temp_Y = temp_X2 = temp_Y2 = 0;

    for(int i=0; i<N_LD; i++) {
        temp_XY = temp_XY + (TF_exp_table[x].LD_exp[i] * gene_exp_table[y].LD_exp[i]);
        temp_X = temp_X + TF_exp_table[x].LD_exp[i];
        temp_Y = temp_Y + gene_exp_table[y].LD_exp[i];
        temp_X2 = temp_X2 + (TF_exp_table[x].LD_exp[i] * TF_exp_table[x].LD_exp[i]);
        temp_Y2 = temp_Y2 + (gene_exp_table[y].LD_exp[i] * gene_exp_table[y].LD_exp[i]);
    }
    
    SUM_XY = temp_XY;
    SUM_X = temp_X;
    SUM_Y = temp_Y;
    SUM_X2 = temp_X2;
    SUM_Y2 = temp_Y2;
    
    if (SUM_X == 0 || SUM_Y == 0) {
        R = 0;
    } else {
        R = (SUM_XY-(SUM_X)*(SUM_Y)/N_LD)/(sqrt((SUM_X2-(SUM_X)*(SUM_X)/N_LD)*(SUM_Y2-(SUM_Y)*(SUM_Y)/N_LD)));
    }
    return R;
}

void node_pair_generator_LD() {

    double R_LD;

    //for(int i=0; i<num_of_TFs; i++) {
    //    for(int j=0; j<num_of_genes; j++) {
	//		for(int k=0; k<10; k++) {
	//			TF_exp_table[i].network[k] = 0;
	//			gene_exp_table[j].network[k] = 0;
	//		}
	//	}
	//}

        for(int i=0; i<num_of_TFs; i++) {
            for(int j=0; j<num_of_genes; j++) {
                if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                    R_LD = r_calculator(i,j);
    		    
                    if(R_LD >= pos_cutoff_LD) {
                        num_of_pos_edge++;
                    }
                }
            }
        }
    

    Pos_Coexp = new R_table[num_of_pos_edge];
    
    for (int i=0; i<num_of_pos_edge; i++) {
        Pos_Coexp[i].checked = 0;
    }
    
    int index_pos = 0;
    //int index_neg = 0;
    for(int i=0; i<num_of_TFs; i++) {
        for(int j=0; j<num_of_genes; j++) {
            if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                R_LD = r_calculator(i,j);

                    if(R_LD >= pos_cutoff_LD) {
                        strcpy(Pos_Coexp[index_pos].query_gene_ID, TF_exp_table[i].gene_ID);
                        strcpy(Pos_Coexp[index_pos].target_gene_ID, gene_exp_table[j].gene_ID);
                        Pos_Coexp[index_pos].pcc = R_LD;
                        index_pos++;
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
                //Pos_Coexp[i].checked = 1;
                Pos_Coexp[i].checked = 2;
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
        printf("\nCan't find seed ID in the GCN. Please check the see ID again!\n\n");
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
    
    fout3 = fopen("Node_level.csv","w");
    
    fprintf(fout3, "TF gene ID,level in GCN\n");
    for(int i=0; i<num_of_TFs; i++) {
        if(TF_exp_table[i].level >= 0) {
		if(TF_exp_table[i].level == 1) {
            		fprintf(fout3, "%s,%d\n", TF_exp_table[i].gene_ID, TF_exp_table[i].level);
		} else {
			fprintf(fout3, "%s,%d\n", TF_exp_table[i].gene_ID, TF_exp_table[i].level - 1);
		}
        }
    }
    
    fclose(fout3);
}

void print_relation_table () {
    FILE *fout1;
    
    fout1 = fopen("Node_relation.csv","w");
    
    fprintf(fout1, "node_1 ID, node_2 ID, PCC_value\n");
    
    for(int i=0; i<num_of_pos_edge; i++) {
        if(Pos_Coexp[i].checked == 1) {
            fprintf(fout1, "%s,%s,%lf\n", Pos_Coexp[i].query_gene_ID, Pos_Coexp[i].target_gene_ID, Pos_Coexp[i].pcc);
        }
    }
    
    fclose(fout1);
}

int main(int argc, char* argv[]) {
    char input_file1[100];   //TF gene list
    char input_file2[100];   //All gene list
    char input_file3[100];   //seed gene list
    
    if (argc != 5) {
        printf("\nUsage: TO-GCN_single No_of_features expression_table_file seed_table_file Cutoff\n");
        
    } else {

        num_of_point_LD = atoi(argv[1]);
        
        strcpy(input_file1, argv[2]);
        strcpy(input_file2, argv[2]);
        strcpy(input_file3, argv[3]);
        
        pos_cutoff_LD = atof(argv[4]);

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

            printf("Number of nodes: %d\n", num_of_TFs);
            printf("Number of features: %d\n", num_of_point_LD);
            printf("Number of seeds: %d\n", num_of_seeds);
            printf("Cutoffs for correlation: %1.3lf\n\n", pos_cutoff_LD);
            printf("Assigning levels for nodes in GCN by Breadth-First-Search (BFS) method......\n");
        
            node_pair_generator_LD();
            level_assignment();
            function_three();
            print_relation_table();
        
            printf("Done!\n");
        }
    }
    
    return 0;	
}
