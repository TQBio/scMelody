# scMelody
Here, we propose scMelody (Single-Cell Methylation data Ensemble cLustering based On the Diversity and separabilitY of cell-to-cell similarity patterns), which utilizes an ensemble clustering framework to link multiple methylation similarity patterns of cells for identifying cell types.

## The overview of scMelody

![image](https://github.com/TQBio/scMelody/blob/main/scMelody/pictures/scMelody-main.png)

# Tutorial
Overview: all source code for implementing scMelody is included in this project, source_code file contains all the functions that will be used, after loading these functions, scMelody can be run following the tutorials. Implementation of scMelody requires both R and python running environments, note the configuration in different environments.

## Step_1 Extract methylation patterns of cells through multiple pairwise similarity measures

This step is implemented in the R environment and is mainly used to calculate cell-to-cell similarity matrices. Input your single-cell methylation dataset in the specified format. The recommended format is an R list, where each element is a dataframe containing methylation information for a single cell, organized as follows:

    Chr      location    methylation_state
   
    chr1      131009             1
    
    chr1      131031             0
    
    ...        ...              ...
    
    chrY     56770241            1
    

Run the following instruction to get the similarity matrix by parallel calculation.

    source('./PiarwiseSimilarity.R')
    
    cl <- makeCluster(8) # the number of CPU cores for computing
    
    clusterExport(cl,"mat_F",envir = environment())
    
    results <- parLapply(cl,1:length(data_list),get_res,data_list) #data_list:single-cell methylation profiles organized as above
        
    stopCluster(cl)
    
    res_cor <- do.call(cbind,results)


## Step_2 Perform spectral clustering to generate initial results

This step is implemented in the python environment and is mainly used to produce spectral clustering assignments for different similarity matrices and estimate the optimal number of clusters. Note that scMelody also provides source code of implementing spectral clustering in R, see file Spectral_clusR.r in the scMelody.

Load all functions in "Individual_clust.py" and your pre-computed similarity matrices files in python environment; 

    K_clust = find_kcluster(k_min,k_max,Cosine_Mat,Hamming_Mat,Pearson_Mat) 
    
    #k_min,k_max represents the minimum and maximum number of possible clusters respectively;
    
    #Cosine_Mat,Hamming_Mat,Pearson_Mat represents the Cosine, Hamming and Pearson similarity matrices respectively;

After getting the optimal cluster number k_opt from K_clust, implement spectral clustering and calculate corresponding weights;

    sc_cosine = sc_pre(Cosine_Mat, k_opt) # spectral clustering assignments of Cosine similarity measure.
    
    sc_hamming = sc_pre(Cosine_Mat, k_opt) # spectral clustering assignments of Hamming similarity measure.
    
    sc_pearson = sc_pre(Cosine_Mat, k_opt) # spectral clustering assignments of Pearson similarity measure.
 
    weight_SIL = weight_sil(Cosine_Mat,Hamming_Mat,Pearson_Mat,sc_cosine, sc_hamming, sc_pearson) 
    
    # get weights for spectral clustering results defined by silhouette coefficient.
    
    weight_PNMI = weight_PNMI(sc_cosine, sc_hamming, sc_pearson)  
    
    # get weights for spectral clustering results defined by pairwise NMI.
    
Save the spectral clustering results of and the corresponding weights in the dataframe as a CSV file. 

    spc = merge_clus(sc_cosine, sc_hamming, sc_pearson,weight_SIL,weight_PNMI)
    
    spc.to_csv('spc.csv')
    
Note: You can download the original similarity matrices of dataset GSE87197 from the "Demo_data" file to reproduce the result or familiarize yourself with the algorithm flow.

## Step_3 Ensemble clustering

This step is implemented in the R environment and is mainly used to calculate the co-association matrix and the finall hierarchical clustering results.
  
     df_clus <- read.csv(file = 'spc.csv',header = T) #load spectral clustering results and the corresponding weights
     
     source('./Ensemble.R')
     
     co_associationMat <- Get_CAMat(df_clus) #calculate the co-association matrix
     
     pre_cluster <- cutree(hclust(co_associationMat),k_opt) #output the final cluster assignments 
        
## 4 Performance evaluation

This step is used to evaluate clustering performance between the truth and the prediction. 

Calculate ARI and V-measure in python. Load the corresponding functions in "Perforamce_evaluation.py" in python.

     ARI_score = ARI(true_label, pre_cluster)
     
     Vm _score = V_measure(true_label, pre_cluster)
     
## 5 Visualization of the co-association matrix
     
Using hierarchical clustering heatmap to visualize the weighted consensus matrix. Here, the co-association matrix obtained by scMelody for real data GSE87197 is provided as a in Demo_data file. Load the corresponding functions in "Visualization.r" in R environment.

Using the co-association matrix for real dataset GSE87197 as an example.

     cell_label <- read.csv(file = 'celllabels_GSE87197.csv')
     
     co-associationMat <- read.csv(file = 'coassMat_GSE87197.csv', header=F)
     
     HC_heatmapGSE87197 <- hc_heatmap(co-associationMat)

![image](https://github.com/TQBio/scMelody/blob/main/scMelody/pictures/Heatmap_of_gse87197.png)

For the generation of synthetic datasets, please refer to https://github.com/andreaskapou/Melissa

If you have any questions in use, please feel free to give me feedback. You can use the following contact information to communicate with me.

Email: tqglowing@std.uestc.edu.cn
