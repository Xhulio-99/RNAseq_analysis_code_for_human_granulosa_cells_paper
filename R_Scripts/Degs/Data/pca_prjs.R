prj_list <- c("PRJNA200696","PRJNA216966","PRJNA649934")

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
  # Regularized log transformation for different analysis (clustering, heatmaps, etc)
  
  dd <- DESeq(dds)
  rld <- rlogTransformation(dd)
  pca <- plotPCA(rld,intgroup=c("source_name"))
  ggsave(paste0("pca_",prj,".png"))
}

