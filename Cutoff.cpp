#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h> 

int num_of_genes;
int num_of_TFs;

long long histogram_LD[201];
long long histogram_TD[201];

int num_of_point_LD;
int num_of_point_TD;

char input_file_TF[100];   //TF gene list
char input_file_gene[100];   //All gene list

typedef struct Time_Course{
        char gene_ID[20];
        double *LD_exp;
        double *TD_exp;
        int network[10];
}Time_Course;

Time_Course *gene_exp_table;
Time_Course *TF_exp_table;

void Read_Time_Course_Data_TFs (char *input) {
	FILE *fptr = fopen(input, "r");
	char GID[20];
	double LDE[num_of_point_LD];
	double TDE[num_of_point_TD];

	num_of_TFs = 0;
	//while(fscanf(fptr,"%s", &GID) != EOF) {
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
	//while(fscanf(fptr,"%s", &GID) != EOF) {
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
	//while(fscanf(fptr,"%s", &GID) != EOF) {
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
	//while(fscanf(fptr,"%s", &GID) != EOF) {
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

void histogram_calculation (double x, int z) {
    int y;
    y = (int)((x + 1) * 100);
    
    if(z == 0) {
        histogram_LD[y]++;
    } else {
        histogram_TD[y]++;
    }
}

void node_pair_generator_LD_or_TD() {

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

    for(int i=0; i<num_of_TFs; i++) {
        for(int j=0; j<num_of_genes; j++) {
            if(strcmp(TF_exp_table[i].gene_ID, gene_exp_table[j].gene_ID) != 0) {
                R_LD = r_calculator(i,j,0);
                R_TD = r_calculator(i,j,1);
                //fprintf(fout9,"%f\t%f\n",R_LD, R_TD);
                histogram_calculation(R_LD,0);
                histogram_calculation(R_TD,1);
            }
        }
    }
}

void function_one () {
    FILE *fout1;
    double bin[201];
    double PDF_LD[201];
    double PDF_TD[201];
    double CDF_asc_LD[201];
    double CDF_asc_TD[201];
    double CDF_desc_LD[201];
    double CDF_desc_TD[201];

    long long sum_LD = 0;
    long long sum_TD = 0;
    
    bin[0] = -1;
    
    for(int i=1; i<=200; i++) {
        bin[i] = bin[i-1] + 0.01;
    }
    
    for(int i=0; i<=200; i++) {
        sum_LD = sum_LD + histogram_LD[i];
        sum_TD = sum_TD + histogram_TD[i];
    }
    
    //printf("%lld, %lld\n", sum_LD, sum_TD);
    
    for(int i=0; i<=200; i++) {
        PDF_LD[i] = (double) histogram_LD[i] / (double) sum_LD;
        PDF_TD[i] = (double) histogram_TD[i] / (double) sum_TD;
    }
    
    CDF_asc_LD[0] = PDF_LD[0];
    CDF_asc_TD[0] = PDF_TD[0];
    
    for(int i=1; i<=200; i++) {
        CDF_asc_LD[i] = CDF_asc_LD[i-1] + PDF_LD[i];
        CDF_asc_TD[i] = CDF_asc_TD[i-1] + PDF_TD[i];
    }
    
    CDF_desc_LD[200] = PDF_LD[200];
    CDF_desc_TD[200] = PDF_TD[200];
    
    for(int i=199; i>=1; i--) {
        CDF_desc_LD[i] = CDF_desc_LD[i+1] + PDF_LD[i];
        CDF_desc_TD[i] = CDF_desc_TD[i+1] + PDF_TD[i];
    }
    
    fout1 = fopen("PCC_histogram.tsv","w");
    
    /*
    fprintf(fout1, "BIN\t# in LD\t# in TD\tPDF in LD\tPDF in TD\tCDF_asc in LD\tCDF_asc in TD\n");
    for(int i=0; i<=200; i++) {
        fprintf(fout1, "%lf\t%lld\t%lld\t%lf\t%lf\t%lf\t%lf\n", bin[i], histogram_LD[i], histogram_TD[i], PDF_LD[i], PDF_TD[i], CDF_asc_LD[i], CDF_asc_TD[i]);
    }
    */

    fprintf(fout1, "BIN\t# in LD\t# in TD\tPDF in LD\tPDF in TD\tCDF_asc in LD\tCDF_asc in TD\tCDF_desc in LD\tCDF_desc in TD\n");
    for(int i=0; i<=200; i++) {
        fprintf(fout1, "%lf\t%lld\t%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", bin[i], histogram_LD[i], histogram_TD[i], PDF_LD[i], PDF_TD[i], CDF_asc_LD[i], CDF_asc_TD[i], CDF_desc_LD[i], CDF_desc_TD[i]);
    }
    
    fclose(fout1);
    
    int pos_cutoff_LD, pos_cutoff_TD;
    int neg_cutoff_LD, neg_cutoff_TD;
    int found_LD = 0;
    int found_TD = 0;
    
    for (int i=200; i>=1; i--) {
        if(found_LD == 0) {
            if(CDF_desc_LD[i] > 0.05) {
                pos_cutoff_LD = i;
                found_LD = 1;
            }
        }
        
        if(found_TD == 0) {
            if(CDF_desc_TD[i] > 0.05) {
                pos_cutoff_TD = i;
                found_TD = 1;
            }
        }
    }
    
    found_LD = 0;
    found_TD = 0;
    for (int i=0; i<=200; i++) {
        if(found_LD == 0) {
            if(CDF_asc_LD[i] > 0.05) {
                neg_cutoff_LD = i;
                found_LD = 1;
            }
        }
        
        if(found_TD == 0) {
            if(CDF_asc_TD[i] > 0.05) {
                neg_cutoff_TD = i;
                found_TD = 1;
            }
        }
    }
    
    pos_cutoff_LD = pos_cutoff_LD + 2;
    pos_cutoff_TD = pos_cutoff_TD + 2;
    
    printf("Cutoff of postive coexpression for Cond. 1 = %1.2lf\n", bin[pos_cutoff_LD]);
    printf("Cutoff of postive coexpression for Cond. 2 = %1.2lf\n", bin[pos_cutoff_TD]);
    printf("Cutoff of negative coexpression for Cond. 1 = %1.2lf\n", bin[neg_cutoff_LD]);
    printf("Cutoff of negative coexpression for Cond. 2 = %1.2lf\n\n", bin[neg_cutoff_TD]);
    printf("You can just copy and paste the following command to get the GCNs or change parameters by yourself:\n");
    printf("GCN %d %d %s %s %1.2lf %1.2lf %1.2lf %1.2f\n", num_of_point_LD, num_of_point_TD, input_file_TF, input_file_gene, bin[pos_cutoff_LD], bin[pos_cutoff_TD], bin[neg_cutoff_LD], bin[neg_cutoff_TD]);
}

 
int main(int argc, char* argv[]) {
    char input_file1[100];   //TF gene list
    char input_file2[100];   //All gene list
    
    if (argc != 5) {
        printf("\nUsage: Cutoff #Cond1_samples #Cond2_samples file_of_TF_genes file_of_all_genes\n\n");
    } else {
        
        for (int i=0; i<=200; i++) {
            histogram_LD[i] = 0;
            histogram_TD[i] = 0;
        }
        
        num_of_point_LD = atoi(argv[1]);
        num_of_point_TD = atoi(argv[2]);
        
        strcpy(input_file1, argv[3]);
        strcpy(input_file2, argv[4]);
        strcpy(input_file_TF, argv[3]);
        strcpy(input_file_gene, argv[4]);
        
        Read_Time_Course_Data_TFs(input_file1);
        Read_Time_Course_Data_genes(input_file2);


        printf("No. of TFs: %d\n", num_of_TFs);
        printf("No. of Genes: %d\n", num_of_genes);
        printf("No. of samples under Cond. 1: %d\n", num_of_point_LD);
        printf("No. of samples under Cond. 2: %d\n\n", num_of_point_TD);
        
        printf("Waiting for histogram generation and cutoff values calculation......\n\n");
        node_pair_generator_LD_or_TD();
        function_one();
        
    }
    
    return 0;	
}
