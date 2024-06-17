library(DESeq2)
library(readxl)
library(ggfortify)
library(genefilter)
library(corrr)
library(ggcorrplot)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(enrichR)
library(multiMiR)
library(clusterProfiler)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
library(dplyr)
library(tidyverse)

PADJ=0.05
FC_THR=1
PADJ_HEAT=0.05
get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}
setwd("data")
PRJNA200699=read.csv("PRJNA200699_miRNAs_expressed_all_samples_20_03_2023.csv",sep = "\t")
PRJNA417973=read.csv("PRJNA417973_miRNAs_expressed_all_samples_20_03_2023.csv",sep = "\t")
setwd("..")
PRJNA200699[c(1,5:length(PRJNA200699))]
str_detect(colnames(PRJNA200699),pattern = ".norm")
mirnas_list <- list("PRJNA200699"=PRJNA200699,"PRJNA417973"=PRJNA417973)
clean_df <- function(df){
  df <- df %>% filter(total>50) %>% rename("miRNA"=X.miRNA)
  #non faccio la somma perchè in mirdeep2 la stessa read viene allineata a più precursori,
  #quindi la funzione mantiene solo il miRNA con conta max tra i miRNA che sono stati
  #allineati su precursori diversi:
  #miDeep2 runs bowtie with -a (report all alignments) and --best --strata, which means the same 
  #read will be mapped to multiple precursor miRNAs
  #df <- df %>% group_by(miRNA) %>% summarise_each(funs(sum)) <- per fare la somma
  df <- df %>% group_by(miRNA) %>% slice_max(total,with_ties = FALSE)
  df <- df[c(1,5:length(df))]
  df <- df[!str_detect(colnames(df),pattern = ".norm")]
  return(df)
}
mirnas_list <- lapply(mirnas_list,clean_df)
mirnas_list[["PRJNA200699"]] %>% filter(miRNA=="hsa-miR-149-5p")

demir <- list()
result_list <- list()
dds_mat_list <- list()
for (prj in names(mirnas_list)){
  prj1 <- mirnas_list[[prj]]
  count_matrix <- as.matrix(prj1[2:length(prj1)])
  rownames(count_matrix) <- prj1$miRNA
  coldata <- as.data.frame(read_excel("data/mir_phenotype_data.xlsx", sheet = prj))
  
  #coldata$source_name <- factor(x = coldata$source_name,levels = c('cumulus_granulosa_cells','mural_granulosa_cells'))
  all(coldata$Sample %in% colnames(count_matrix))
  all(coldata$Sample == colnames(count_matrix))
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = coldata,
                                design = ~ source_name)
  dds_mat_list[[prj]] <- dds
  dd <- DESeq(dds)
  #vstd <- vst(dds)
  rld <- rlog(dds)
  plotPCA(rld, intgroup = "source_name", ntop = 500) + ggtitle(paste('RLog Transformation ',prj,sep = "")) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste("pca_",prj,".png",sep = ""))
  source_names <- unique(coldata$source_name)
  res <- results(dd,c("source_name",source_names[1],source_names[2]))
  resSig <- subset(res, padj < PADJ & abs(log2FoldChange)>=FC_THR)
  resSig <- resSig[order(abs(resSig$log2FoldChange),decreasing=T), ]
  result_list[[prj]] <- res
  demir[[prj]] <- as.data.frame(resSig)
  write.csv(as.data.frame(resSig), 
            file=paste("deMIR_",prj,".csv",sep = ""))
}


demir[["PRJNA200699"]]
str(demir)
#View(demir[["PRJNA200699"]])
for (prj in names(result_list)){
  pheno_data <- as.data.frame(read_excel("data/mir_phenotype_data.xlsx", sheet = prj))
  groups <- unique(pheno_data$source_name)
  EnhancedVolcano(result_list[[prj]],lab = result_list[[prj]]@rownames,
                  pCutoff = 5e-02,
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
                  title = paste("DEMirs",prj,groups[1],"vs",groups[2]),
                  gridlines.minor = F,
                  drawConnectors = T,
                  widthConnectors = 0.6,
                  lengthConnectors = unit(0.013, "npc"),
                  legendPosition = "right",
                  subtitle = "",
                  caption = "",
                  max.overlaps = 7,
                  boxedLabels = T)
  ggsave(paste("volcano_",prj,".jpg",sep = ""),scale = 2,width = 5,height = 6)
  
  deseq2VST <- vst(dds_mat_list[[prj]],nsub=nrow(dds_mat_list[[prj]]))             # estimate dispersion trend and apply a variance stabilizing transformation (find a simple function ƒ to apply to values x in a data set to create new values y = ƒ(x) such that the variability of the values y is not related to their mean value)
  deseq2VST <- assay(deseq2VST)
  idx <- rownames(result_list[[prj]])[ which(result_list[[prj]]$padj < PADJ_HEAT & abs(result_list[[prj]]$log2FoldChange)>=FC_THR) ]
  rownames(pheno_data) <- pheno_data$Sample
  ann_colors <- list (source_name=c("purple1","limegreen"))
  ann_colors <- setNames(ann_colors$source_name,groups)
  ann_colors <- list(source_name=ann_colors)
  
  ann_col=data.frame(row.names = pheno_data$Sample,source_name=factor(pheno_data$source_name))
  
  mirna_to_order_heatmap <- intersect(rownames(demir[[1]]),rownames(demir[[2]]))
  
  ##cambiato ordine della heatmap
  ##per ripristinare basta rimettere al posto di ordered_ -> deseq2VST[idx, ]
  if(nrow(deseq2VST[idx, ])!=0)
  {
    ordered_ <- as.data.frame(deseq2VST[idx, ] , stringsAsFactors = FALSE)
    ordered_ <- ordered_ %>% mutate(order_col=ifelse(rownames(ordered_) %in% 
                                                       mirna_to_order_heatmap, 1, 2)) %>% 
      arrange(order_col, rownames(ordered_)) %>% select(-order_col)
    ordered_ <- as.matrix(ordered_)
    heat_map <- pheatmap(ordered_,scale = "row",
                      cluster_rows = F, 
                     clustering_method = "complete",
                     annotation_col = ann_col,
                     annotation_colors = ann_colors,
                     angle_col = "0",
                     main = prj)
                     #fontsize = 130,
                     #fontsize_col = 40,
                     #fontsize_row = 40)
    plot_dims <- get_plot_dims(heat_map)
    base_file_name=paste0("heatmap_",prj,"_no_cluster_row")
    pdf(file =paste0(base_file_name,".pdf"),height = plot_dims$height*3.5, 
        width = plot_dims$width*1.5)
    print(heat_map)
    dev.off()
    ## per convertire da command line con imagemagick pdf in jpg
    system2("convert",paste0("-density 300 ",base_file_name,".pdf",
                   " -quality 100 ",base_file_name,".jpg"))
  }
}

library(cooltools)





scale(count_matrix)
corr_matrix <- cor(scale(count_matrix))
ggcorrplot(corr_matrix)




df_tot <- inner_join(mirnas_list$PRJNA200699,mirnas_list$PRJNA417973,by="miRNA",suffix=c("_PRJNA200699","_PRJNA417973"))
coldata_tot <- as.data.frame(read_excel("data/mir_phenotype_data.xlsx", sheet = "tot"))
count_matrix_tot <- as.matrix(df_tot[2:length(df_tot)])
rownames(count_matrix_tot) <- df_tot$miRNA
dds_tot <- DESeqDataSetFromMatrix(countData = count_matrix_tot,
                                  colData = coldata_tot,
                                  design = ~ source_name + bioproject)


