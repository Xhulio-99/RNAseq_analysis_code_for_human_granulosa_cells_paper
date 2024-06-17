library(rlist)
library(dplyr)

intersect(rownames(demir[[1]]),rownames(demir[[2]]))

library(multiMiR)
multimir_result <- list()
for (prj in names(demir)){
multimir_result[[prj]] <- get_multimir(org     = 'hsa',
                                       mirna = rownames(demir[[prj]]),
                                       table   = 'validated',
                                       summary = FALSE)
}
all_target <- list()
for (prj in names(multimir_result)){
  all_target[[prj]] <- unique(multimir_result[[prj]]@data$target_symbol)
}



mir_target_df_list <- list()
for (prj in names(multimir_result)){
  mir_target_df <- multimir_result[[prj]]@data %>% select(database,mature_mirna_id,target_symbol,target_entrez,target_ensembl,pubmed_id) %>%
    group_by(mature_mirna_id,target_entrez) %>% slice_head(n = 1) %>% ungroup() %>% 
      filter(target_symbol!="")
  mir_target_df <- mir_target_df %>% group_by(target_symbol) %>% mutate(occ=n()) %>% ungroup()
  mir_target_df_list[[prj]] <- mir_target_df
}


go_list <- list()
for (prj in names(all_target)){
  occm5 <-  unique((mir_target_df_list[[prj]] %>% filter(occ>5))$target_symbol)
  ego <- enrichGO(gene       = occm5,
                  universe      = all_target[[prj]],
                  OrgDb         = hs,
                  keyType       ="SYMBOL",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  go_list[[prj]] <- clusterProfiler::simplify(ego)
}


for (prj in names(go_list)){
  go <- go_list[[prj]]
  bp <- go@result$ONTOLOGY=="BP"
  cc <- go@result$ONTOLOGY=="CC"
  mf <- go@result$ONTOLOGY=="MF"
  
  dotplot(filter(go,bp),showCategory=20,title=paste("Biological Process",prj))
  ggsave(paste("bp_",prj,".jpg",sep = ""),scale = 1.5,height = 8)
  
  dotplot(filter(go,cc),showCategory=20,title=paste("Cellular Component",prj))
  ggsave(paste("cc_",prj,".jpg",sep = ""),scale = 1.5,height = 8)
  
  dotplot(filter(go,mf),showCategory=20,title=paste("Molecular Function",prj))
  ggsave(paste("mf_",prj,".jpg",sep = ""),scale = 1.5,height = 8)
  
}

all_target_occ5 <- list()
for(prj in names(mir_target_df_list)){
  all_target_occ5[[prj]] <- unique((mir_target_df_list[[prj]] %>% filter(occ>5))$target_symbol)
}

ccg <- compareCluster(all_target_occ5,fun = "enrichGO",
               OrgDb         = org.Hs.eg.db,
               keyType       ="SYMBOL",
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05,
               readable      = TRUE)

dotplot(ccg)
ggsave("compare_cluster.jpg",width = 7)
# library(createKEGGdb)
# create_kegg_db(species = "hsa")
#DOWNLOAD FROM "https://www.bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz"
#install.packages("KEGG.db_3.2.4.tar.gz", repos=NULL,type="source")
library(KEGG.db)
kegg_list <- list()
for (prj in names(all_target)){
  all_ids <- bitr(all_target[[prj]], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  occm5 <- unique(mir_target_df_list[[prj]] %>% filter(occ>5))$target_symbol
  occm5_ids <- bitr(occm5, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  ekk <- enrichKEGG(gene =occm5_ids$ENTREZID,
                    universe = all_ids$ENTREZID,
                    keyType = "kegg",
                  organism = 'hsa',
                  pvalueCutoff = 0.05,
                  use_internal_data=T)
  kegg_list[[prj]] <- ekk
}

for (prj in names(kegg_list)){
  kegg <- kegg_list[[prj]]
  dotplot(kegg,showCategory=20,title=paste("Enriched Keggs",prj))
  ggsave(paste("kegg_",prj,".jpg",sep = ""),scale = 1.5,height = 8)
  
}


mirna_target_list=list()
for (prj in names(multimir_result)){
  mir_target_df <- multimir_result[[prj]]@data %>% select(database,mature_mirna_id,target_symbol,target_entrez,target_ensembl,pubmed_id) %>%
    group_by(mature_mirna_id,target_entrez) %>% slice_head(n = 1) %>% ungroup()
  mir_target_df <- mir_target_df %>% group_by(target_symbol) %>% mutate(occ=n()) %>% ungroup() %>% 
    filter(target_symbol!="")
  mirnas <- unique(mir_target_df$mature_mirna_id)
  targets_per_mirna <- list()
  for (mir in mirnas){
    targets <- unique(mir_target_df[mir_target_df$mature_mirna_id==mir,]$target_symbol)
    targets_per_mirna[[mir]] <- targets
  }
  mirna_target_list[[prj]] <- targets_per_mirna
}

library(rlist)

library(igraph)
library(qgraph)
library(ComplexHeatmap)
library(RColorBrewer)
create_network <- function(inputlist,geni_filter,fc_df,file_name,title,legend_title){
  adj_matrix <- list_to_matrix(inputlist)
  adj_matrix <- adj_matrix[rownames(adj_matrix) %in% geni_filter,]
  adj_matrix <- adj_matrix[,colSums(adj_matrix)>0]
  network <- graph_from_incidence_matrix(adj_matrix)
  V(network)$shape <- ifelse(V(network)$type ,"sphere", "circle")
  #V(network)$shape[V(network)$type & V(network)$name %in% mirna_8_occ] <- "sphere"
  V(network)$size <- ifelse(V(network)$type ,9, 6)
  #V(network)$size[V(network)$type & V(network)$name %in% mirna_8_occ] <- 15
  #colori MIR
  v_proj_names <- V(network)$name[V(network)$type]
  fc <- fc_df[v_proj_names,]$log2FoldChange
  mir_color <- colorRampPalette(brewer.pal(n = 10, name = "RdYlBu"))(length(fc))
  V(network)$color[V(network)$type] <- mir_color
  #colori geni
  V(network)$color[!V(network)$type] <- "white"
  
  edge_list <- get.edgelist(network)
  names(mir_color) <- v_proj_names
  coled <- c()
  col_cex <- ifelse(V(network)$type,"white","black")
  cex_size <- ifelse(V(network)$type,6,5)
  for (i in edge_list[,2]){
    coled <- append(coled,mir_color[[i]])
  }
  
  e <- get.edgelist(network,names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(network),
                                         area=30*(vcount(network)^2),
                                         repulse.rad=(vcount(network)^2))
  
  lgd_= rep(NA,length(fc)/2)
  lgd_[c(1,length(fc)/4,length(fc)/2)] = sapply(c(max(fc),median(fc),min(fc)),FUN = function(x){format(round(x, 2), nsmall = 2)})
  
  pdf(file = paste0(file_name,".pdf"), width = 200, height = 230)
  plot(network,layout=l, vertex.label.font = 2,
       vertex.label.color= col_cex, vertex.label.cex= cex_size, 
       edge.curved=0, edge.width=8,edge.color=coled)
  title(main=title,cex.main=15,line=-2)
  legend("bottomright",
        title = legend_title,
         legend = lgd_,
         fill = colorRampPalette(brewer.pal(n = 10, name = "RdYlBu"))(length(fc)/2),
         border = NA,
         y.intersp = 0.5,
         cex = 10, text.font = 2)
  dev.off()
}

TARGET_OCC=15
for (prj in names(demir)){
  geni_filter_network <- unique(mir_target_df_list[[prj]] %>% filter(occ>TARGET_OCC))$target_symbol
  df_fc <- demir[[prj]]
  create_network(mirna_target_list[[prj]],geni_filter_network,df_fc,
                file_name = paste(prj,"_network_target_mir",sep = ""),
                title = paste("miRNA-Target Network",prj),
                legend_title = "miRNA log2FC (ccs vs mgcs)")

}
prj <- "PRJNA200699"
df_fc <- demir[[prj]]

degs_PRJNA200696 <- read.csv("data/DEG_PRJNA200696.csv")
degs_PRJNA200696 <- degs_PRJNA200696 %>% rename("gene"=X)
down_regulated_degs <- degs_PRJNA200696[degs_PRJNA200696$log2FoldChange<0,]$gene
up_regulated_mirna <- rownames(demir[["PRJNA200699"]][demir[["PRJNA200699"]]$log2FoldChange>0,])
PRJNA200699_mirna_target_list <- mirna_target_list[["PRJNA200699"]][up_regulated_mirna]
PRJNA200699_mirna_target_list <- sapply(PRJNA200699_mirna_target_list,function(x){x[x%in%down_regulated_degs]})
PRJNA200699_mirna_target_list <- PRJNA200699_mirna_target_list[!is.na(names(PRJNA200699_mirna_target_list))]

create_network(PRJNA200699_mirna_target_list,down_regulated_degs,df_fc,
               file_name = paste(prj,"_network_target_mir_filter-down-mRNA",sep = ""),
               title = paste("CCs Up-regulated miRNA-Target Network",prj),
               legend_title = "miRNA log2FC (ccs vs mgcs)")


degs_PRJNA200696 <- read.csv("data/DEG_PRJNA200696.csv")
degs_PRJNA200696 <- degs_PRJNA200696 %>% rename("gene"=X)
down_regulated_degs_MGCs <- degs_PRJNA200696[degs_PRJNA200696$log2FoldChange>0,]$gene
up_regulated_mirna_MGCs <- rownames(demir[["PRJNA200699"]][demir[["PRJNA200699"]]$log2FoldChange<0,])
#PRJNA200699_mirna_target_list <- mirna_target_list[["PRJNA200699"]][up_regulated_mirna_MGCs]
common_miRNAS <- intersect(dems[[1]], dems[[2]])
PRJNA200699_mirna_target_list <- mirna_target_list[["PRJNA200699"]][common_miRNAS]
PRJNA200699_mirna_target_list <- sapply(PRJNA200699_mirna_target_list,function(x){x[x%in%down_regulated_degs_MGCs]})
PRJNA200699_mirna_target_list <- PRJNA200699_mirna_target_list[!is.na(names(PRJNA200699_mirna_target_list))]


base_file_name = paste0(prj,"_network_target_common_miRNAs")
create_network(PRJNA200699_mirna_target_list,down_regulated_degs_MGCs,df_fc,
               file_name = base_file_name,
               title = paste("Common mirna miRNA-Target Network",prj),
               legend_title = "miRNA log2FC (ccs vs mgcs)")

system2("convert",args = paste0("-density 300 ",base_file_name,".pdf",
                         " -quality 100 ",base_file_name,".jpg"))
system2()
