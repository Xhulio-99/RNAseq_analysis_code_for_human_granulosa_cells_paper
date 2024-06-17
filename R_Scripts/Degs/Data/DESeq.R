library(DESeq2)
library(dplyr)
library(readxl)
library(genefilter)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)

PADJ=0.05
FC_THR=1
PADJ_HEAT=1e-15

get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

#prj_list <- c("PRJNA200696","PRJNA263182","PRJNA216966","PRJNA649934")
#prj_list <- c("PRJNA200696")
prj_list <- c("PRJNA200696","PRJNA216966","PRJNA649934")
#pdf(file ="plots4.pdf",width = 10,height = 17)
dir.create(file.path("plots"), showWarnings = FALSE)
prj_list <- c("PRJNA649934")
for (prj in prj_list){
  gene_count_matrix <- read.csv(paste("gene_count_matrix_",prj,".csv",sep = ""), row.names=1)
  rownames(gene_count_matrix) <- sapply(strsplit(rownames(gene_count_matrix),
                                                 split = "|",fixed = TRUE),`[`, 1)
  
  count_data <- as.matrix(gene_count_matrix)
  #gene_count_matrix <- gene_count_matrix %>% filter(rowSums(.) > 40)
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  
  #Apply it to my matrix:
  filter1 <-  apply( count_data, 1, fun)
  count_data <- count_data[filter1,]
  coldata <- as.data.frame(read_excel("phenotype_data.xlsx", sheet = prj))
  
  #coldata$source_name <- factor(x = coldata$source_name,levels = c('cumulus_granulosa_cells','mural_granulosa_cells'))
  all(coldata$Sample %in% colnames(count_data))
  all(coldata$Sample == colnames(count_data))
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = coldata,
                                design = ~ source_name)

  dd <- DESeq(dds)
  
  source_names <- unique(coldata$source_name)
  res <- results(dd,c("source_name",source_names[1],source_names[2]))
  
  resSig <- subset(res, padj < PADJ & abs(log2FoldChange)>=FC_THR)
  #sort and show by log2 fold change value:
  resSig <- resSig[order(abs(resSig$log2FoldChange),decreasing=T), ]
  
  write.csv(as.data.frame(resSig), 
            file=paste("DEG_",prj,".csv",sep = ""))
  
# VOLCANO PLOT  
  pdf(file =paste("plots/prova_volcano_",prj,".pdf",sep = ""),height = 10,width = 9)
  print(EnhancedVolcano(resSig,lab = resSig@rownames,
                  pCutoff = 1e-05,
                  pointSize = 4,
                  labSize = 6,
                  labFace = 'bold',
                  colAlpha = 0.8,
                  #colCustom = keyvals,
                  x = 'log2FoldChange',
                  y = 'padj',
                  ylab = bquote(~-Log[10] ~ italic(P-adj)),
                  legendLabels = c("NS", expression(Log[2] ~ FC), "p-adj", expression(p-adj ~ and
                                                                                        ~ log[2] ~ FC)),
                  #xlim = c(-3,3),
                  #ylim = c(0,25),
                  title = paste("DEGS",prj,source_names[1],"vs",source_names[2]),
                  gridlines.minor = F,
                  drawConnectors = T,
                  widthConnectors = 0.6,
                  lengthConnectors = unit(0.013, "npc"),
                  legendPosition = "right",
                  subtitle = "",
                  caption = "",
                  max.overlaps = 10,
                  boxedLabels = T))
  dev.off()
  # PRE-PROCESS DATI PER PHEATMAP
  
  deseq2VST <- vst(dds)                # estimate dispersion trend and apply a variance stabilizing transformation (find a simple function ƒ to apply to values x in a data set to create new values y = ƒ(x) such that the variability of the values y is not related to their mean value)
  deseq2VST <- assay(deseq2VST)
  idx <- rownames(res)[ which(res$padj < PADJ_HEAT & abs(res$log2FoldChange)>=FC_THR) ]
  rownames(coldata) <- coldata$Sample
  
  ann_colors <- list (source_name=c("purple1","limegreen"))
  ann_colors <- setNames(ann_colors$source_name,source_names)
  ann_colors <- list(source_name=ann_colors)
  
  ann_col=data.frame(row.names = coldata$Sample,source_name=factor(coldata$source_name))
  if(nrow(deseq2VST[ idx, ])!=0)
  {
    heat_map <- pheatmap(deseq2VST[ idx, ],scale = "row",
                         clustering_distance_rows = "correlation", 
                         clustering_method = "complete",
                         annotation_col = ann_col,
                         annotation_colors = ann_colors,
                         angle_col = "0",
                         main = prj,
                         fontsize = 130,
                         fontsize_col = 40,
                         fontsize_row = 40)
    
    plot_dims <- get_plot_dims(heat_map)
    pdf(file =paste("plots/prova_heatmap_",prj,".pdf",sep = ""),height = plot_dims$height*25, 
        width = plot_dims$width*15)
    print(heat_map)
    dev.off()
  } else {
    print(paste("the number of rows for the heatmap input matrix is ",nrow(deseq2VST[ idx, ])))
    print("table of DEGs is empty or it is not a dataframe: no heatmap will be produced")
  }
}
dev.off()









#PCA on PRJNA263182 dataset...
pca_res <- prcomp(t(gene_count_matrix))
autoplot(pca_res,data=t(gene_count_matrix), label = TRUE ,label.size=2)



