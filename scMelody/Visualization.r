
###Using UMAP and hierarchical clustering for visualization in R###

library(umap)
library(pheatmap)


UMAP_plot <- function(CA,cell_label){
     em_bd <- umap(1-CA,input='dist',n_components=2)$layout
     df_plot <- data.frame(cell_label,em_bd[,1],em_bd[,2])
     colnames(df_plot) <- c("Data","Dimension1","Dimension2")
     df_plot$Data<- factor(df_plot$Data)
     plot_umap <- ggplot(df_plot, aes(x = Dimension1, y = Dimension2,color=Data,shape=Data)) +
                                        geom_point(size=3)+
            		        scale_color_brewer(palette = "Set2")+
            		        scale_shape_manual(values = c(0,1,2,5,6,7))
            		        theme_set(theme_bw())
    return(plot_umap)
   }

hc_heatmap <- function(CA){
   return(pheatmap(as.matrix(CA),cutree_rows = NA,cutree_cols=NA,
                    cellwidth = 8, cellheight = 8))
}