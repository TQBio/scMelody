library(pheatmap)
library(diceR)
library(Nbclust)
library(diceR)
library(stats)

Get_CAMat <- function(df_clus){
    weight_sil <- clus[c(nrow(clus)-1),]
    weight_PNMI <- clus[c(nrow(clus)),]
    sc_res <-  clus[-c(nrow(clus)-1,nrow(clus)),]
    cons_mat <- 0.5 * (consensus_matrix(sc_res,weight=weight_sil)+consensus_matrix(sc_res,weight=weight_PNMI))
    return(cons_mat)
}
