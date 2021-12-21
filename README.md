# scMelody
Here, we present scMelody, which utilizes an enhanced consensus-based clustering model to reconstruct cell-to-cell methylation similarity patterns and identifies cell subpopulations with the leveraged information from multiple basic similarity measures.

## The overview of scMelody

![image](https://github.com/TQBio/scMelody/blob/main/scMelody/pictures/Fig1.png)

# Tutorial
Overview: all source code for implementing scMelody is included in this project, source_code file contains all the functions that will be used, after loading these functions, scMelody can be run following the tutorials. Implementation of scMelody requires both R and python running environments, note the configuration in different environments.

## 1 Extract methylation patterns of cells through multiple pairwise similarity measures

This step is implemented in the R environment and is mainly used to calculate cell-to-cell similarity matrices. Input your single-cell methylation dataset in the specified format. The recommended format is an R list, where each element is a dataframe containing methylation information for a single cell, organized as follows:

    Chr      location    methylation_state
   
    chr1      131009             1
    
    chr1      131031             0
    
    ...        ...              ...
    
    chrY     56770241            1
    

We provide code that can run multithreaded in R environment to calculate the basic similarity matrix between cells. Run the following instruction to get the similarity matrix by parallel calculation.

    source('./PiarwiseSimilarity.R')
    
    cl <- makeCluster(8) # the number of CPU cores for computing
    
    clusterExport(cl,"mat_F",envir = environment())
    
    results <- parLapply(cl,1:length(data_list),get_res,data_list) #data_list:single-cell methylation profiles organized as above
        
    stopCluster(cl)
    
    res_cor <- do.call(cbind,results)


## 2 Perform spectral clustering to generate initial results

This step is implemented in the python environment and is mainly used to produce the estimated the optimal number of clusters and the resulting cell partitions. Note that scMelody also provides source code of implementing spectral clustering in R, see file Spectral_clusR.r in the scMelody.

Load all functions in "Step_2.py" and your pre-computed similarity matrices files in python environment; 

    k_sp = find_kspcluster(k_max,Cosine_Mat,Hamming_Mat,Pearson_Mat) ### the optimal number of clusters for spectral clustering 
    
    res_mat = scM(Cosine_Mat,Hamming_Mat,Pearson_Mat,k_sp) ### the reconstructed similarity matrix 
    
    k_opt =  find_kopcluster(k_max,res_mat) ### the optimal number of clusters for the resulting cell partitions
    
    C = hc_pre(res_mat,k_opt) ### the resulting cell partitions
    
Note: You can download the original similarity matrices of dataset GSE87197 from the "Demo_data" file to reproduce the result and familiarize yourself with the algorithm flow.

       
## 3 Performance evaluation

This step is used to evaluate clustering performance between the truth and the prediction. 

Calculate ARI and V-measure in python. Load the corresponding functions in "Perforamce_evaluation.py" in python.

     ARI_score = ARI(true_label, pre_cluster)
     
     Vm _score = V_measure(true_label, pre_cluster)


## 4 Visualization of the reconstructed similarity matrix
     
Using hierarchical clustering heatmap to visualize the reconstructed similarity matrix.

Taking the real dataset GSE87197 as an example.

     cell_label <- read.csv(file = 'celllabels_GSE87197.csv')
        
     HC_heatmapGSE87197 <- hc_heatmap(res_mat)

![image](https://github.com/TQBio/scMelody/blob/main/scMelody/pictures/Heatmap_Farlik2016.png)


## 5 Note

(1) For the generation of synthetic datasets, please refer to https://github.com/andreaskapou/Melissa

(2) If you have any questions in use, please feel free to give me feedback. You can use the following contact information to communicate with me.

    Email: tqglowing@std.uestc.edu.cn
