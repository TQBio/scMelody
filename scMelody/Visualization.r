
###Using hierarchical clustering heatmap to visualize the weighted consensus matrix in R###

###Taking the dataset of GSE87197 for example####

library(pheatmap)

cell_bd  <- c(rep('CLP',times=21),rep('CMP',times=19),rep('GMP',times=22),
              rep('HSC',times=18),rep('MLP0',times=24),rep('MPP',times=18))

inx_1 <- c(paste('GSM23244',43:63,sep =''))
inx_1 <- c(paste(inx_1,"CLP",1:21,sep ='_'))
inx_2 <- c(paste('GSM2324',497:515,sep =''))
inx_2 <- c(paste(inx_2,"CMP",1:19,sep ='_'))
inx_3 <- c(paste('GSM2324',549:570,sep =''))
inx_3 <- c(paste(inx_3,"GMP",1:22,sep ='_'))
inx_4 <- c(paste('GSM2324',625:642,sep =''))
inx_4 <- c(paste(inx_4,"HSC",1:18,sep ='_'))
inx_5 <- c(paste('GSM2324',815:838,sep =''))
inx_5 <- c(paste(inx_5,"MLP0",1:24,sep ='_'))
inx_6 <- c(paste('GSM2324',977:994,sep =''))
inx_6 <- c(paste(inx_6,"MPP",1:18,sep ='_'))
inx_bd<- c(inx_1,inx_2,inx_3,inx_4,inx_5,inx_6)

hc_heatmap <- function(CA){
   return(pheatmap(as.matrix(CA),cutree_rows = NA,cutree_cols=NA,
                    cellwidth = 8, cellheight = 8))
}

bd_df <- data.frame(cons_Mat)
rownames(bd_df) <- inx_bd
colnames(bd_df) <- inx_bd
annotation <- data.frame(cell_bd,pre_cluster) 
colnames(annotation) <- c("true clusters","inferred clusters")
rownames(annotation) <- inx_bd
tiff(filename = "Heatmap_of_gse87197.tiff", res = 600, width = 4500, height = 4000,compression = 'lzw')
pheatmap(as.matrix(bd_df),cutree_rows = NA,cutree_cols=NA,
         annotation_col = annotation,cellwidth = 3, cellheight = 3, fontsize =4)
dev.off()
