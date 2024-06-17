#Compare cluster di tutti i DEGS(per prj)
#Selezionare dal compare cluster i geni up e down regolati nel sigolo progetto
#Mettere al posto dei dot un grafico a torta con % di up regolati nei diversi tipi cellulari

library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(genefilter)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(stringr)

degs_list_up_down = list(
  #DEG cumulus GV vs M2
  up_reg_PRJNA216966_GV=(read.csv("DEG_PRJNA216966.csv") %>% filter(log2FoldChange>0))[[1]],
  down_reg_PRJNA216966_GV_UP_M2=(read.csv("DEG_PRJNA216966.csv") %>% filter(log2FoldChange<0))[[1]],
  
  #DEG oocyte vs cumulus
  up_reg_PRJNA649934_CC=(read.csv("DEG_PRJNA649934.csv") %>% filter(log2FoldChange<0))[[1]],
  down_reg_PRJNA649934_CC=(read.csv("DEG_PRJNA649934.csv") %>% filter(log2FoldChange>0))[[1]],
  
  #DEG cumulus vs mural
  up_reg_PRJNA200696_CC=(read.csv("DEG_PRJNA200696.csv") %>% filter(log2FoldChange>0))[[1]],
  down_reg_PRJNA200696_CC=(read.csv("DEG_PRJNA200696.csv") %>% filter(log2FoldChange<0))[[1]]
)

degs_list = list(
  #DEG cumulus GV vs M2
  PRJNA216966=(read.csv("DEG_PRJNA216966.csv"))[[1]],
  
  #DEG oocyte vs cumulus
  PRJNA649934=(read.csv("DEG_PRJNA649934.csv"))[[1]],
  
  #DEG cumulus vs mural
  PRJNA200696=(read.csv("DEG_PRJNA200696.csv"))[[1]]
)
hs <- org.Hs.eg.db
convert_to_ENTREZID <- function(genes){
  r=AnnotationDbi::select(hs,genes, columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")
  return(r$ENTREZID)
}
degs_list=lapply(degs_list,convert_to_ENTREZID)
compare_cluster <- compareCluster(geneCluster = degs_list, fun = "enrichGO",
                                      OrgDb         = hs,
                                      keyType       ="ENTREZID",
                                      readable      = TRUE)
write.csv(compare_cluster,file = "compare_cluster_all_prjs.csv")
compare_cluster@compareClusterResult$Description
p <- dotplot(compare_cluster)
#Get the actual data from the dotplot (selects the GOs to represent)
df_compare_cluster <- p$data
df_compare_cluster$Cluster
df_compare_cluster <- df_compare_cluster %>% mutate(ClusterName=str_split_fixed(Cluster,"\n",2)[,1])
names(degs_list_up_down)

#Split the string in vector of genes(enriched for a GO)
#Get the percentage of genes(on the lenght of the enriched genes vector) in the named degs list element 
get_perc <- function(name,gene_list){
  genes <- str_split(gene_list,pattern = "/")[[1]]
  return((sum(genes %in% degs_list_up_down[[name]])/length(genes))*100)
}

#This function gets the list of vector of genes(in string format separated by "/")
#Appends to cc or others vector the percentage of genes in up or down regulated degs vector
#for each bioproject
get_cc_percent <- function(prjs,gene_list){
  cc <- c()
  others <- c()
  gv <- c()
  m2 <- c()
  print("CIAOO")
  for( i in seq_along(prjs)){
    
  switch (prjs[i],
    "PRJNA216966" = {gv <- c(gv,get_perc("up_reg_PRJNA216966_GV",gene_list[i]))
                     m2 <- c(m2,get_perc("down_reg_PRJNA216966_GV_UP_M2",gene_list[i]))
                     cc <- c(cc,0)
                     others <- c(others,0)
                     
                     # print(str_split(gene_list[i],pattern = "/")[[1]])
                     print(gv)
                     print(m2)
                     },
    "PRJNA649934" = {cc <- c(cc,get_perc("up_reg_PRJNA649934_CC",gene_list[i]))
                     others <- c(others,get_perc("down_reg_PRJNA649934_CC",gene_list[i]))
                     gv <- c(gv,0)
                     m2 <- c(m2,0)
                     
                     },
    "PRJNA200696" = {cc <- c(cc,get_perc("up_reg_PRJNA200696_CC",gene_list[i]))
                     others <- c(others,get_perc("down_reg_PRJNA200696_CC",gene_list[i]))
                     gv <- c(gv,0)
                     m2 <- c(m2,0)
                     
                     }
    )
  }
  # print(cc)
  res_l = list(cc=cc,others=others,gv=gv,m2=m2)
  print(res_l)
  return(res_l)
}
df_compare_cluster <- p$data
df_compare_cluster$Cluster
df_compare_cluster <- df_compare_cluster %>% mutate(ClusterName=str_split_fixed(Cluster,"\n",2)[,1])
names(degs_list_up_down)
colnames(df_compare_cluster)
df_compare_cluster <- df_compare_cluster %>% mutate(CC=get_cc_percent(ClusterName,geneID)[[1]],Others=get_cc_percent(ClusterName,geneID)[[2]],
                                                    CC_GV=get_cc_percent(ClusterName,geneID)[[3]],CC_M2=get_cc_percent(ClusterName,geneID)[[4]])
# library(ggforce)
library(scatterpie)

# Plot
fig <- function(width, heigth){
  options(repr.plot.width = width, repr.plot.height = heigth)
}
options(repr.plot.width = 30, repr.plot.height =10)
dots <- ggplot() +
  geom_point(data = df_compare_cluster, aes(x = Cluster, y = Description),color = "transparent") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  ) +
  labs(title = "Go enrichment analysis",
       x = "Bioprojects",
       y = "GOs",
       fill = "Cell Type")
dots
#To get the coords of the dots
# options(repr.plot.width = 20, repr.plot.height =3)
dd <- ggplot_build(dots)
#Getting x and y coords and creating index column
x_y <- dd$data[[1]] %>% select(x,y) %>% mutate(index=rownames(dd$data[[1]]))
#Resetting df index for the full join
row.names(df_compare_cluster) <- NULL
pie_data <- x_y %>% full_join(df_compare_cluster %>% mutate(index=rownames(df_compare_cluster)),by = "index") %>%  select(x,y,CC,Others,CC_GV,CC_M2)
dots+
  geom_scatterpie(aes(x=x, y = y),data = pie_data, cols = c("CC","Others","CC_GV","CC_M2"),pie_scale = 7)
ggsave("Compare cluster PIE_CHARTS.jpg",width = 6,height = 8)

library(settings)
reset(options)
