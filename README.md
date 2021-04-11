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
    

Run "source()"
