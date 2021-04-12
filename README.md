# scMelody
scMelody is a clustering tool for single-cell methylation data analysis. It is an ensemble clustering approach to evaluate epigenetic heterogeneity
based on genome-wide methylation similarity between cell pairs. It adopts an ensemble clustering framework to aggregate individual clustering results from spectral clustering on multiple intercellular methylation similarity matrices, and yields a co-association matrix for hierarchical clustering using an adaptive weight assignment strategy.

# Tutorial
Overview: all source code for implementing scMelody is included in this project, source_code file contains all the functions that will be used, after loading these functions, scMelody can be run following the tutorials. Implementation of scMelody requires both R and python running environments, note the configuration in different environments.

## 1 Extract methylation patterns of cells through multiple pairwise similarity measures

This step is implemented in the R environment and is mainly used to calculate cell-to-cell similarity matrices. Input your single-cell methylation dataset in the specified format. The recommended format is an R list, where each element is a dataframe containing methylation information for a single cell, organized as follows:

    Chr      location    methylation_state
   
    chr1      131009             1
    
    chr1      131031             0
    
    ...        ...              ...
    
    chrY     56770241            1
    

Run the following instruction to get the similarity matrix by parallel calculation.

    source('./PiarwiseSimilarity.R')
    
    cl <- makeCluster(8) 
    
    clusterExport(cl,"mat_F",envir = environment())
    
    results <- parLapply(cl,1:length(data_list),get_res,data_list) #data_list:single-cell methylation profiles organized as above
        
    stopCluster(cl)
    
    res_cor <- do.call(cbind,results)


## 2 Perform spectral clustering to generate initial results

This step is implemented in the python environment and is mainly used to produce clustering assignments for different similarity matrices and number of clusters. Here, the subsequent flow is demonstrated using the pre-calculated similarity matrices from the gse87197 dataset as an example.

Load all functions in "Individual_clust.py" and similarity matrices files {Cosine_Mat.csv,Hamming_Mat.csv,Pearson_Mat.csv} in python environment; 

    K_clust = find_kcluster(k_min,k_max,Cosine_Mat,Hamming_Mat,Pearson_Mat) #k_min,k_max represents the minimum and maximum number of possible clusters;

After getting the optimal cluster number k_opt from K_clust, implement spectral clustering and calculate corresponding weights;

    sc_cosine = sc_pre(Cosine_Mat, k_opt) # spectral clustering assignments of Cosine similarity measure.
    
    sc_hamming = sc_pre(Cosine_Mat, k_opt) # spectral clustering assignments of Hamming similarity measure.
    
    sc_pearson = sc_pre(Cosine_Mat, k_opt) # spectral clustering assignments of Pearson similarity measure.
 
    weight_SIL = weight_sil(Cosine_Mat,Hamming_Mat,Pearson_Mat,sc_cosine, sc_hamming, sc_pearson) # get weights for spectral clustering results defined by silhouette coefficient.
    
    weight_PNMI = weight_PNMI(sc_cosine, sc_hamming, sc_pearson)  # get weights for spectral clustering results defined by pairwise NMI.
    
Save the spectral clustering results of and the corresponding weights in the dataframe las CSV files. 

    spc = merge_clus(sc_cosine, sc_hamming, sc_pearson,weight_SIL,weight_PNMI)
    
    spc.to_csv('spc.csv')
    
## 3 Ensemble clustering

This step is implemented in the R environment and is mainly used to calculate the co-association matrix and the finall hierarchical clustering results.
  
     df_clus <- read.csv(file = 'spc.csv',header = T) #load spectral clustering results and the corresponding weights
     
     source('./Ensemble.R')
     
     co_associationMat <- Get_CAMat(df_clus) #calculate the co-association matrix
     
     pre_cluster <- cutree(hclust(co_associationMat),k_opt) #output the final cluster assignments 
        
## 4 Performance evaluation

This step is used to evaluate clustering performance between the truth and the prediction. It also includes visualization process for UMAP and hierarchical clustering heatmap.

Calculate ARI and V-measure in python. Load the functions in "results_analysis.txt".

     ARI1 = ARI(true_label, pre_cluster)
     
     V_measure1 = V_measure(true_label, pre_cluster)
     
Using UMAP and hierarchical clustering for visualization. Load the functions in "results_analysis.txt".

     UMAP_GSE87197 <- UMAP_plot(co_associationMat,true_label,"GSE87197")
     
     HC_heatmapGSE87197 <- hc_heatmap(co_associationMat)


If you have any questions in use, please feel free to give me feedback. You can use the following contact information to communicate with me.

Email: tqglowing@std.uestc.edu.cn
