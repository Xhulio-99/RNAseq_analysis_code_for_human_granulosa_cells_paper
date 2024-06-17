library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(genefilter)
library(RColorBrewer)
library(qpdf)
library(ggplot2)
library(dplyr)
library(scatterpie)

prj_list <- c("PRJNA200696","PRJNA216966","PRJNA649934")

# FUNCTION TO GET GOS FROM PRJ_LIST
get_gos <- function(prj_list){
  go_list <- list()
  for (prj in prj_list){
    degs <- rownames(read.csv(paste("DEG_",prj,".csv",sep=""),row.names = 1))
    ###get gene_list
    
    gene_count_matrix <- read.csv(paste("gene_count_matrix_",prj,".csv",sep = ""), row.names=1)
    rownames(gene_count_matrix) <- sapply(strsplit(rownames(gene_count_matrix),
                                                   split = "|",fixed = TRUE),`[`, 1)
    
    count_data <- as.matrix(gene_count_matrix)
    fun <- kOverA(round(dim(count_data)[2]*80/100),5)
    filter1 <-  apply( count_data, 1, fun)
    count_data <- count_data[filter1,]
    gene_list <- rownames(count_data)
    
    
    ego <- enrichGO(gene       = degs,
                    universe      = gene_list,
                    OrgDb         = org.Hs.eg.db,
                    keyType       ="SYMBOL",
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    go_list <- append(go_list,ego)
  }
  names(go_list) <- prj_list
  return(go_list)
}
gos_list <- get_gos(prj_list)
gos_list[["PRJNA649934"]]
# LOOP TO PLOT GO_ENRICHED_LIST
for (go_prj in names(go_list)){
  #go_prj <- "PRJNA216966"
  ego <- gos_list[[go_prj]]
  bp <- ego@result$ONTOLOGY=="BP"
  cc <- ego@result$ONTOLOGY=="CC"
  mf <- ego@result$ONTOLOGY=="MF"
  
  write.csv(as.data.frame(filter(ego,bp)),
             file=paste("go_enrichment/all_ego_bp_",go_prj,".csv",sep = ""))
  write.csv(as.data.frame(filter(ego,cc)),
             file=paste("go_enrichment/all_ego_cc_",go_prj,".csv",sep = ""))
  write.csv(as.data.frame(filter(ego,mf)),
             file=paste("go_enrichment/all_ego_mf_",go_prj,".csv",sep = ""))
  file_name1=paste("plots/go_plot_",go_prj,".pdf",sep = "")
  pdf(file=file_name1,height = 12,width = 9)

  print(dotplot(filter(ego,bp),showCategory=30,title=paste("Biological Process",go_prj)))
  print(dotplot(filter(ego,cc),showCategory=30,title=paste("Cellular Component",go_prj)))
  print(dotplot(filter(ego,mf),showCategory=30,title=paste("Molecular Function",go_prj)))
  dev.off()

  file_name2=paste("plots/go_plot_temp_",go_prj,".pdf",sep = "")
  pdf(file=file_name2,height = 20,width = 26)
  degss <- read.csv(paste("DEG_",go_prj,".csv",sep=""),row.names = 1)
  fc <- degss$log2FoldChange
  names(fc) <- as.character(rownames(degss))
  fc = sort(fc, decreasing = TRUE)
  edox <- setReadable(filter(ego,bp), 'org.Hs.eg.db', 'ENTREZID')
  
#   selected <- c("nuclear division",
# "circulatory system process",
# "regulation of system process",
# "cell cycle phase transition",
# "cell cycle G2/M phase transition")
  
  print(cnetplot(edox, foldChange=fc, circular = TRUE, colorEdge = TRUE,title=go_prj)+
          theme(text = element_text(size = 28),
                legend.text=element_text(size=30),
                legend.key.size = unit(1,"cm"))
        )
  dev.off()
  
  file_out=paste("plots/go_plots_",go_prj,".pdf")
  pdf_combine(input=c(file_name1,file_name2),output = file_out)
  unlink(file_name2)
  unlink(file_name1)
  
}

# LOOP TO CLUSTER COMPARE BIOPROJECTS
degs_list=list()
for (prj in prj_list){
  # get degs from csv
  deg <- read.csv(paste("DEG_",prj,".csv",sep=""))[[1]]
  # convert symbol to entrez_id
  deg <- bitr(deg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  degs_list[[prj]] <- deg$ENTREZID
}
ck <- compareCluster(geneCluster = degs_list, fun = enrichKEGG,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

#pathways_intersection <- ck@compareClusterResult[["Description"]][duplicated(ck@compareClusterResult[["Description"]])]
#df_ck <- ck@compareClusterResult
df_ck <-ck@compareClusterResult %>% group_by(Description) %>% 
  mutate(path_occ=n()) %>% arrange(desc(path_occ),qvalue)

ck@compareClusterResult <-ck@compareClusterResult %>% group_by(Description) %>% 
    mutate(path_occ=n()) %>% arrange(desc(path_occ),qvalue)

head(ck)
pdf("plots/common_kegg.pdf",width = 15,height = 12)
print(dotplot(ck,showCategory=10,title="Common enriched Keggs"))
print(cnetplot(ck,title="Common enriched Keggs"))
dev.off()

pdf("plots/common_kegg.pdf",width = 15,height = 12)
print(dotplot(ck))
print(cnetplot(ck))
dev.off()



cg <- compareCluster(geneCluster = degs_list, fun = "enrichGO",
                     OrgDb         = org.Hs.eg.db,
                     keyType       ="ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
cg <- setReadable(cg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

df_cg <- as.data.frame(cg@compareClusterResult)
df_cg <-df_cg %>% group_by(Description) %>% 
  mutate(path_occ=n()) %>% arrange(desc(path_occ),qvalue)
cg2 <- cg
cg2@compareClusterResult <- cg2@compareClusterResult %>% group_by(Description) %>% 
  mutate(path_occ=n()) %>% arrange(desc(path_occ),qvalue)

### prova uniione cg
l <- list(cg, cg2, cg)
names(l) <- c("cg","cg2","cg")
keys <- unique(unlist(lapply(l, names)))
setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, keys))), keys)

cg_cg <- do.call(append,list(cg,cg))

### togliere legenda gene number
prova <- cnetplot(cg,title="Common enriched Gos")
prova$layers

plot$layers[c(4,5,6,7)] <- NULL
prova$layers
prova
prova <- prova %>% remove_geom("segment") %>% remove_geom("text",idx = 2) %>% 
  remove_geom("arc",idx=3) %>% remove_geom("text",idx = 1)

pdf("plots/common_go_prove_legend.pdf",width = 20,height = 17)
#print(dotplot(cg,showCategory=10,font.size=15,title="Common enriched Gos"))
#print(cnetplot(cg,title="Common enriched Gos")
print(prova)
dev.off()

library(ggedit)


layersList(prova)
prova$layers[[4]]

prova <- remove_geom(prova, "geom_segment")
prova %>% remove_geom("arc") %>% remove_geom("arc")

str(prova, max.level = 2, size = FALSE)

go_intersection <- cg@compareClusterResult[["Description"]][duplicated(ck@compareClusterResult[["Description"]]) & ck@compareClusterResult[["qvalue"]]<1e-03]
pdf("plots/common_go_intersection.pdf",width = 20,height = 17)
print(dotplot(cg,showCategory=go_intersection,title="Intersection and qval < 1e-03"))
print(cnetplot(cg,showCategory=go_intersection,title="Intersection and qval < 1e-03"))
dev.off()


cg_all <- compareCluster(geneCluster = degs_list, fun = "enrichGO",
                     OrgDb         = org.Hs.eg.db,
                     keyType       ="ENTREZID",
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)


head(cg_all,20)
cg_all <- setReadable(cg_all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
pdf("plots/common_go_all.pdf",width = 20,height = 17)
print(dotplot(cg_all,font.size=15))
print(cnetplot(cg_all))
dev.off()



for (prj in prj_list){
  degs <- read.csv(paste("DEG_",prj,".csv",sep=""),row.names = 1)
  up_regulated <- rownames(degs %>% filter(log2FoldChange>0))
  down_regulated <- rownames(degs %>% filter(log2FoldChange<0))
  all <- rownames(degs)
  regs <- list("up"=up_regulated,"down"=down_regulated,"all"=all)
  ###get gene_list
  
  gene_count_matrix <- read.csv(paste("gene_count_matrix_",prj,".csv",sep = ""), row.names=1)
  rownames(gene_count_matrix) <- sapply(strsplit(rownames(gene_count_matrix),
                                                 split = "|",fixed = TRUE),`[`, 1)
  
  count_data <- as.matrix(gene_count_matrix)
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  filter1 <-  apply( count_data, 1, fun)
  count_data <- count_data[filter1,]
  gene_list <- rownames(count_data)
  
  reg_names <- names(regs)
  for (reg in reg_names){
    ego <- enrichGO(gene       = regs[[reg]],
                    universe      = gene_list,
                    OrgDb         = org.Hs.eg.db,
                    keyType       ="SYMBOL",
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    bp <- ego@result$ONTOLOGY=="BP"
    cc <- ego@result$ONTOLOGY=="CC"
    mf <- ego@result$ONTOLOGY=="MF"
    
    write.csv(as.data.frame(filter(ego,bp)),
              file=paste("go_enrichment/",reg,"_ego_bp_",prj,".csv",sep = ""))
    write.csv(as.data.frame(filter(ego,cc)),
              file=paste("go_enrichment/",reg,"_ego_cc_",prj,".csv",sep = ""))
    write.csv(as.data.frame(filter(ego,mf)),
              file=paste("go_enrichment/",reg,"_ego_mf_",prj,".csv",sep = ""))
    
  }
  
}
### plot all only

for (prj in prj_list){
  degs <- rownames(read.csv(paste("DEG_",prj,".csv",sep=""),row.names = 1))
  ###get gene_list
  
  gene_count_matrix <- read.csv(paste("gene_count_matrix_",prj,".csv",sep = ""), row.names=1)
  rownames(gene_count_matrix) <- sapply(strsplit(rownames(gene_count_matrix),
                                                 split = "|",fixed = TRUE),`[`, 1)
  
  count_data <- as.matrix(gene_count_matrix)
  fun <- kOverA(round(dim(count_data)[2]*80/100),5)
  filter1 <-  apply( count_data, 1, fun)
  count_data <- count_data[filter1,]
  gene_list <- rownames(count_data)
  
  
  
  ego <- enrichGO(gene       = degs,
                  universe      = gene_list,
                  OrgDb         = org.Hs.eg.db,
                  keyType       ="SYMBOL",
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  bp <- ego@result$ONTOLOGY=="BP"
  cc <- ego@result$ONTOLOGY=="CC"
  mf <- ego@result$ONTOLOGY=="MF"
  
  pdf(file=paste("plots/go_plot",prj,".pdf"))
  print(dotplot(filter(ego,bp),showCategory=30,title="Biological Process"))
  print(dotplot(filter(ego,cc),showCategory=30,title="Cellular Component"))
  print(dotplot(filter(ego,mf),showCategory=30,title="Molecular Function"))
  dev.off()
  
}

dotplot(ego)





