library(clusterProfiler)
library(org.Hs.eg.db)

prj_list <- c("PRJNA200696","PRJNA216966","PRJNA649934")


for (prj in prj_list){
  
  degs <- read.csv(paste("DEG_",prj,".csv",sep=""),row.names = 1)
  all <- rownames(degs)
  all_ids <- bitr(all, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ###get gene_list
  
  gene_count_matrix <- read.csv(paste("gene_count_matrix_",prj,".csv",sep = ""), row.names=1)
  rownames(gene_count_matrix) <- sapply(strsplit(rownames(gene_count_matrix),
                                                 split = "|",fixed = TRUE),`[`, 1)
  
  count_data <- as.matrix(gene_count_matrix)
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  filter1 <-  apply( count_data, 1, fun)
  count_data <- count_data[filter1,]
  gene_list <- rownames(count_data)
  
  kk_all <- enrichKEGG(gene         = all_ids$ENTREZID,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)
  kk_all <- setReadable(kk_all, 'org.Hs.eg.db', 'ENTREZID')
  write.csv(as.data.frame(kk_all), 
            file=paste("kegg/all_kegg_",prj,".csv",sep = ""))
  pdf(file=paste("plots/kegg_plot",prj,".pdf"),height = 12)
  print(dotplot(kk_all,showCategory=30,title=paste("Kegg",prj)))
  dev.off()
  
}





#prj <- "PRJNA200696"
for (prj in prj_list){
  
  degs <- read.csv(paste("DEG_",prj,".csv",sep=""),row.names = 1)
  up_regulated <- rownames(degs %>% filter(log2FoldChange>0))
  down_regulated <- rownames(degs %>% filter(log2FoldChange<0))
  all <- rownames(degs)
  
  up_ids <- bitr(up_regulated, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  down_ids <- bitr(down_regulated, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  all_ids <- bitr(all, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ###get gene_list
  
  gene_count_matrix <- read.csv(paste("gene_count_matrix_",prj,".csv",sep = ""), row.names=1)
  rownames(gene_count_matrix) <- sapply(strsplit(rownames(gene_count_matrix),
                                                 split = "|",fixed = TRUE),`[`, 1)
  
  count_data <- as.matrix(gene_count_matrix)
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  filter1 <-  apply( count_data, 1, fun)
  count_data <- count_data[filter1,]
  gene_list <- rownames(count_data)
  
  
  
  kk_up <- enrichKEGG(gene         = up_ids$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
  write.csv(as.data.frame(kk_up), 
            file=paste("kegg/up_kegg_",prj,".csv",sep = ""))
  

  kk_down <- enrichKEGG(gene         = down_ids$ENTREZID,
                      organism     = 'hsa',
                      pvalueCutoff = 0.05)
  write.csv(as.data.frame(kk_down), 
            file=paste("kegg/down_kegg_",prj,".csv",sep = ""))
  
  
  kk_all <- enrichKEGG(gene         = all_ids$ENTREZID,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  write.csv(as.data.frame(kk_all), 
            file=paste("kegg/all_kegg_",prj,".csv",sep = ""))
  
}

browseKEGG(kk_up,"hsa04914")
