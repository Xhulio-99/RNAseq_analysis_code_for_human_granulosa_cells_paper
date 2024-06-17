library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(qpdf)
library(genefilter)
HS <- org.Hs.eg.db

get_gos <- function(lista,tipo_lista){
  go_list <- list()
  for (prj in names(lista)){
    count_data <- as.matrix(gene_counts[[prj]])
    fun <- kOverA(round(dim(count_data)[2]*80/100),5)
    filter1 <-  apply(count_data, 1, fun)
    count_data <- count_data[filter1,]
    gene_list <- rownames(count_data)
    lsv_type = unique(lista[[prj]]$gene_name)
    ego <- enrichGO(gene       = lsv_type,
                    #universe      = gene_list,
                    OrgDb         = HS,
                    keyType       ="SYMBOL",
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    
    go_list[[prj]] <- simplify(ego)
  }
  return(go_list)
}


go_lsv <- get_gos(lsv,"ALL")
go_lsv_tf <- get_gos(lsv_tf,"TF")
go_lsv_deg <- get_gos(lsv_filtered,"DEG")
#go_no_filter <- get_gos(lsv_all,"NO_FILTER")

plot_gos <- function(gos_list,tipo_lista){
  for (go_prj in names(gos_list)){
  #go_prj <- "PRJNA216966"
  ego <- gos_list[[go_prj]]
  write.csv(as.data.frame(filter(ego@result)),
            file=paste("go_enrichment/all_ego_bp_",go_prj,"_",tipo_lista,".csv",sep = ""))

    bp <- ego@result$ONTOLOGY=="BP"
    cc <- ego@result$ONTOLOGY=="CC"
    mf <- ego@result$ONTOLOGY=="MF"
  
  
  file_name1=paste("go_enrichment/go_plot_",go_prj,"_",tipo_lista,".pdf",sep = "")
  pdf(file=file_name1,height = 12,width = 9)
  
  if (sum(bp)>0){
    print(dotplot(filter(ego,bp),showCategory=30,title=paste("Biological Process",go_prj)))
  }
  if (sum(cc)>0){
    print(dotplot(filter(ego,cc),showCategory=30,title=paste("Cellular Component",go_prj)))
  }
  if (sum(mf)>0){
    print(dotplot(filter(ego,mf),showCategory=30,title=paste("Molecular Function",go_prj)))
  }
  dev.off()
  
  file_name2=paste("go_enrichment/go_plot_temp_",go_prj,"_",tipo_lista,".pdf",sep = "")
  pdf(file=file_name2,height = 20,width = 26)
  # degss <- read.csv(paste("DEG_",go_prj,".csv",sep=""),row.names = 1)
  # fc <- degss$log2FoldChange
  # names(fc) <- as.character(rownames(degss))
  # fc = sort(fc, decreasing = TRUE)
  if (sum(bp)>0){
    edox <- setReadable(filter(ego,bp), 'org.Hs.eg.db', 'ENTREZID')
    print(cnetplot(edox, circular = TRUE, colorEdge = TRUE,title=go_prj)+
          theme(text = element_text(size = 28),
                legend.text=element_text(size=30),
                legend.key.size = unit(1,"cm"))
    )
  }
  dev.off()
  
  file_out=paste("go_enrichment/go_plots_",go_prj,"_",tipo_lista,".pdf",sep = "")
  pdf_combine(input=c(file_name1,file_name2),output = file_out)
  unlink(file_name2)
  unlink(file_name1)
  
  }
}

plot_gos(go_lsv,"ALL")
plot_gos(go_lsv_tf,"TF")
plot_gos(go_lsv_deg,"DEG")
plot_gos(go_no_filter,"NO_FILTER")


library(plotrix)
pie3D(data, mar = rep(1.75, 4),
      col = hcl.colors(length(data), "Spectral"),
      labels = data,
      explode = 0.2)

kegg_list <- list()
for (prj in names(lsv)){
  lsv_ids <- bitr(lsv[[prj]]$gene_name, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  count_data <- as.matrix(gene_counts[[prj]])
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  filter1 <-  apply(count_data, 1, fun)
  count_data <- count_data[filter1,]
  gene_list <- rownames(count_data)
  gene_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ekk <- enrichKEGG(gene =lsv_ids$ENTREZID,
                    #universe = gene_ids$ENTREZID,
                    keyType = "kegg",
                    organism = 'hsa',
                    #pvalueCutoff = 0.05,
                    use_internal_data=T)
  kegg_list[[prj]] <- ekk
}
for (prj in names(kegg_list)){
  if (sum(kegg_list[[prj]]@result$p.adjust<0.05)>0){
    print(dotplot(kegg_list[[prj]],title=paste("Enriched KEGGs",prj)))
    ggsave(paste("go_enrichment/kegg_",prj,".jpg",sep = ""),scale = 1.5,height = 8)
  }
}
