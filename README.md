# TO-GCN
Pipeline of time-ordered gene coexpression network (TO-GCN) construction from three-dimensional (gene expression, condition, and time) data

The pipeline of constructing a time-ordered gene coexpression network (TO-GCN) contains three steps: (1) Determining the cutoff values, (2) constructing eight GCNs for different types, and (3) determining time-ordered level for nodes in interesting GCN.

## Prepare the gene expression data

Before going to the pipeline, we need to prepare two lists (TF genes and all genes) with RPKM values at different sample points under two conditions, the three-dimensional data. In addition to the data files, we also prepare some information of data like the number of sample under condition 1 (n1) and condition (n2).

In the example folder, there are two data files from the study of "A Comparative Transcriptomics Method to Infer Time-ordered Gene Coexpression Networks and applications". The data should be a Tab-separated values (.tsv) format that contains m rows and n columns, where m is the number of genes (TF genes or all genes) and n represents the summation (n1 + n2) of sample number under condition 1 (n1) and condition 2 (n2).The gene ID is listed in the first column. For each gene, the RPKM values of each sample point under condition 1 and condition 2 are listed from the second to (n1+1)-th columns and from (n1+2)-th to (n1+n2+1)-th columns, respectively. In the example data of TFs_1718.tsv, there are 1718 rows for 1718 TF genes and 27 columns for one gene ID, 13 samples of condition 1, and 13 samples of condition 2.

## Run the program of pipeline

As mentioned above, there are three step for the pipeline. Therefore, we provided program for each step: (1) Cutoff, (2) GCN, and (3) TO-GCN. You can directly run the program by downloading the corresponding binary codes for different system platforms, Linux, MacOS, or Windows. You can also download the C++ source code (.cpp) and compile to executable one by yourself. For compiling source codes by yourself, you can use the following commands:

> g++ Cutoff.cpp -o Cutoff

> g++ GCN.cpp -o GCN

> g++ TO-GCN.cpp -o TO-GCN

### (1) Determining the cutoff values

We need cutoff values of Pearson’s Correlation Coefficient (PCC) under two conditions for constructing the GCN. Our method is to calculate all the PCC values for each TF-gene pair for each condition. With these PCC values, we generate the distribution of probability density function (PDF) and cumulative density function (CDF). According to the CDF, we can give you a suggested cutoff value whose p < 0.05 for each condition. To run the Cutoff program, you have to give four parameters: number of samples under condition 1, number of samples under condition 2, data file of TF genes, and data file of all genes. Here is the example of our study:

> Cutoff 13 13 example_data/TFs_1718.tsv example_data/ All_genes_25489.tsv

### (2) Constructing eight GCNs for different types

In this step, we want to construct eight coexpression types of GCN under two conditions (C1 and C2): C1+C2+, C1+C20, C1+C2–, C10C2+ C1–C2+, C1–C2–, C1–C20, and C10C2–. The output file of each GCN is listed in comma-separated value (.csv) format. The five columns represent the TF gene ID, coexpression type, gene ID, PCC under condition 1, PCC under condition 2. You can import these gene pair into the network generation tool, like Cytoscape, to get the visualization of the GCN. To run the GCN program, you have to give four more parameters of positive and negative cutoff under two conditions. Here is the example:

> GCN 13 13 example_data/TFs_1718.tsv example_data/TFs_1718.tsv 0.84 0.84 -0.75 -0.75

### (3) Determining time-ordered level in the interesting GCN

The final step is to determine the time-order (level) of nodes in the GCN. The time-order is assigned with a seed node by the breadth-first search (BFS) algorithm. In most case, we will select a gene that highly expressed in the first time point and lowly expressed in the following time points as a seed node. In our study, we select a gene with ID, Zm00001d041056, and run the TO_GCN program to assign the time-order (level) of nodes in C1+C2+ GCN. Therefore, we need two more parameters, the seed node gene ID and the coexpression type (0, 1, or 2) where 0, 1, and 2 represent the C1+C2+, C1+C20, and C10C2+, respectively. The level of each node and GCN can be both imported into the Cytoscape. 

> TO-GCN 13 13 example_data/TFs_1718.tsv example_data/TFs_1718.tsv 0.84 0.84 Zm00001d041056 0
