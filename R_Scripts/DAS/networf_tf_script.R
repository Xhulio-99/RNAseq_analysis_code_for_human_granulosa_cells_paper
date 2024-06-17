library(ComplexHeatmap)
library(RColorBrewer)
library(igraph)
library(qgraph)
library(colourvalues)

tf_genes <- read.delim("data/trrust_rawdata.human.tsv",header = FALSE,col.names = c("TF","Gene","Activity","ID"))
tf_degs <- list()
for (prj in names(lsv_tf)){
  tf_degs[[prj]] <- tf_genes %>% filter(TF %in% lsv_tf[[prj]]$gene_name,Gene %in% rownames(degs[[prj]]))
  tf_degs[[prj]]["Gene_FC"] <- sapply(tf_degs[[prj]]$Gene,FUN = function(x) degs[[prj]][rownames(degs[[prj]])==x,]$log2FoldChange)
  
  edge_list <- as.matrix(tf_degs[[prj]][,1:2])
  network <- graph_from_edgelist(edge_list,directed = TRUE)
  V(network)$type <- V(network)$name %in% edge_list[,1]
  V(network)$shape <- ifelse(V(network)$type ,"rectangle", "circle")
  V(network)$size <- ifelse(V(network)$type ,16, 8)
  #colori PRJ
  tf_color <- colorRampPalette(brewer.pal(9, "Set1"))(sum(V(network)$type))
  V(network)$color[V(network)$type] <- tf_color
  #colori geni
  genes = V(network)$name[!V(network)$type]
  fcs = sapply(genes,FUN = function(x) degs[[prj]][rownames(degs[[prj]])==x,]$log2FoldChange)
  heat_pal <- grDevices::colorRamp(c("blue","yellow","red"))(1:100/100)
  V(network)$color[!V(network)$type] <- color_values(fcs,palette = heat_pal)
  #V(network)$color[!V(network)$type] <- colorRampPalette(brewer.pal(11, "RdYlBu"))(diff(range(fcs)))
  
  coled <- ifelse(tf_degs[[prj]]$Activity=="Activation","green",ifelse(tf_degs[[prj]]$Activity=="Repression","red","grey"))
  
  e <- get.edgelist(network,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(network),
                                         area=30*(vcount(network)^2),
                                         repulse.rad=(vcount(network)^2))
  #l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(network))
  #plot
  LEN_LEGEND = 30
  lgd_= rep(NA,length(LEN_LEGEND))
  lgd_[c(1,LEN_LEGEND/2,LEN_LEGEND)] = c(max(fcs),median(fcs),min(fcs))
  pdf(file = paste("network_",prj,".pdf",sep = ""), width = 50, height = 50)
  plot(network,layout=l, vertex.label.font = 2,
       vertex.label.color= c("black"), vertex.label.cex= 5, 
       edge.curved=0, edge.width=8,edge.color=coled)
  legend("bottomright",
         legend = lgd_,
         fill = colorRampPalette(colors = c("red","yellow","blue"))(LEN_LEGEND),
         border = NA,
         y.intersp = 0.5,
         cex = 4, text.font = 2)
  dev.off()
}


go_list_tf_targets <- list()
for (prj in names(tf_degs)){
  count_data <- as.matrix(gene_counts[[prj]])
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  filter1 <-  apply(count_data, 1, fun)
  count_data <- count_data[filter1,]
  gene_list <- rownames(count_data)
  tf_targets = tf_degs[[prj]]$Gene
  ego <- enrichGO(gene       = tf_targets,
                  universe      = gene_list,
                  OrgDb         = HS,
                  keyType       ="SYMBOL",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  go_list_tf_targets[[prj]] <-ego
}



go_list_tf_targets$PRJNA200696@result <- go_list_tf_targets$PRJNA200696@result %>% filter(!Description %in% c("response to lipopolysaccharide",
                                                                                                             "response to molecule of bacterial origin",
                                                                                                             "response to bacterium",
                                                                                                             "learning"))
#go_list_tf_targets$PRJNA216966@result <- go_list_tf_targets$PRJNA216966@result %>% filter(!Description %in% c(""))
#go_list_tf_targets$PRJNA649934@result <- go_list_tf_targets$PRJNA649934@result %>% filter(!Description %in% c(""))

plot_gos(go_list_tf_targets,"TF_TARGETS")

#tf_degs$PRJNA216966 %>% dplyr::filter(Gene=="ESR1")

library(plotrix)
library(dplyr)
library(gsubfn)
library(cowplot)
library(stringr)
library(stringi)
library(ggplot2)

choose_high_fc <- function(genes,prj){
  genes <- paste(stri_split(genes,fixed = "/")[[1]])
  if (prj=="PRJNA649934"){
    res <-  tf_degs[[prj]] %>% arrange((Gene_FC)) %>% filter(Gene %in% genes) %>% distinct(Gene)
    res <- (unname(unlist(res)[1:2]))
    return (paste(res,collapse = "\n"))
  }
  res <-  tf_degs[[prj]] %>% arrange(desc(Gene_FC)) %>% filter(Gene %in% genes) %>% distinct(Gene)
  res <- (unname(unlist(res)[1:2]))
  return (paste(res,collapse = "\n"))
  
}

labs <- sapply(data$geneID,FUN = choose_high_fc,"PRJNA649934")
labs

paste(stri_split(data$geneID,fixed = "/")[[1]])
GENE_TO_TEST <- c("CYP51A1",
                  "TNC",
                  "EPHB1",
                  "FBXO32",
                  "ADAMTS16",
                  "GRIN2A",
                  "CLDN11",
                  "CMAS")


for (prj in names(go_list_tf_targets)){
  data <- go_list_tf_targets[[prj]]@result
  data <- data %>% mutate(GeneRatio_numeric=as.numeric(gsubfn("(\\d+)/(\\d+)", ~ as.numeric(x)/as.numeric(y),GeneRatio))) %>% arrange(desc(GeneRatio_numeric)) %>% slice_head(n=10)
  print(data)
  #labs <- sapply(data$geneID,FUN = function(x){paste(stri_split(x,fixed = "/")[[1]][1:2],collapse = "\n")})
  labs <- sapply(data$geneID,FUN = choose_high_fc,prj)
  print(labs)
  pdf(file = paste("go_enrichment/FC_order_tf_targets_GO_pie3D_",prj,".pdf",sep=""),width = 9,height = 7)
  pie3D(data$GeneRatio_numeric, mar = rep(1.75, 4),
        col = hcl.colors(length(data$GeneRatio_numeric), "Spectral"),
        labels = labs,
        explode = 0.2)
  par(xpd=TRUE)
  legend(1,0.7,legend=data$Description,cex=0.8,xjust = 0.9,yjust = 4.5,
         fill = hcl.colors(length(data$GeneRatio_numeric), "Spectral"),ncol = 2)
  dev.off()
  
}
tf_degs[["PRJNA200696"]]

kegg_list_tf_targets <- list()
for (prj in names(tf_degs)){
  tf_degs_ids <- bitr(tf_degs[[prj]]$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  count_data <- as.matrix(gene_counts[[prj]])
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  filter1 <-  apply(count_data, 1, fun)
  count_data <- count_data[filter1,]
  gene_list <- rownames(count_data)
  gene_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ekk <- enrichKEGG(gene =tf_degs_ids$ENTREZID,
                    universe = gene_ids$ENTREZID,
                    keyType = "kegg",
                    organism = 'hsa',
                    #pvalueCutoff = 0.05,
                    use_internal_data=T)
  kegg_list_tf_targets[[prj]] <- ekk
}
for (prj in names(kegg_list_tf_targets)){
  if (sum(kegg_list_tf_targets[[prj]]@result$p.adjust<0.05)>0){
    print(dotplot(kegg_list_tf_targets[[prj]],title=paste("Enriched KEGGs",prj)))
    ggsave(paste("go_enrichment/kegg_tf_targets_",prj,".jpg",sep = ""),scale = 1.5,height = 8)
  }
}



