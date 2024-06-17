PRJNA200696 <- as.matrix(unique(read.table("DAS_1.txt", header = T)))

PRJNA649934 <- as.matrix(unique(read.table("DAS_649934.txt", header = T)))

PRJNA216966 <- as.matrix(unique(read.table("DAS_216966.txt", header = T)))

#convert from SYMBOL to ENTREZ:
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
gcSample <- lsv_all
lsv_all$PRJNA200696
#Mural-Cumulus

names(gcSample) <- c("PRE-OV(M2)_CCs-vs-MCs","CCs_GV-vs-M2","GV_Oocyte-Vs-CCs")

bitr(gcSample[[1]]$gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = hs)[,2]
gcSample <- lapply(gcSample, function(df) bitr(df$gene_name,fromType = "SYMBOL",
                                   toType = "ENTREZID",OrgDb = hs)[,2])
names(gcSample)
#Genero le comparazioni tra cluster di GO e plotto solo quelle condivise da piu configurazioni:
library(clusterProfiler)
ck <- compareCluster(geneCluster = gcSample,  fun = "enrichGO", OrgDb = "org.Hs.eg.db")
ck <- setReadable(ck, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
head(ck) 
write.csv(ck,"GO-DAS-all-info-11-10-23.csv") 
library(clusterProfiler)
library(createKEGGdb)

createKEGGdb::create_kegg_db("hsa")
install.packages("KEGG.db_1.0.tar.gz", repos=NULL,type="source")
library(KEGG.db)

#KEGG:
kg <- compareCluster(geneCluster = gcSample,  fun = "enrichKEGG", organism="hsa")
kg <- setReadable(kg, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
head(kg) 
write.csv(kg,"KEGG-DAS-11-10-23.csv")
clusterProfiler:::kegg_list("hsa")
library(dplyr)
#filter(x, p.adjust < .05) 

#parte grafica:
library(ggplot2)
names(gcSample) <- c("PRE-OV(M2)_CCs-vs-MCs","CCs_GV-vs-M2","GV_Oocyte-Vs-CCs")
mycolors <- c("CCs_GV-vs-M2"= "#619dff", "PRE-OV(M2)_CCs-vs-MCs" = "#f9766e", "GV_Oocyte-Vs-CCs" = "#01ba38")
cn<- cnetplot(ck, cex_label_gene=2,cex_gene=3, cex_label_category=3, showCategory = 10)
cn+scale_color_manual(mycolors)
cn$layers[c(4,5,6,7)] <- NULL
cn
ggsave("cnetplot-DAS-GO-10categories.jpeg",scale=5)

dp1<-dotplot(ck)
dp1 + theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
ggsave("dotplot-DAS-GO-Fig2B.jpeg")
  
dp2<-dotplot(kg)
cn2 <- cnetplot(kg,cex_label_gene=1,cex_gene=3, cex_label_category=2)
cn2$layers[c(4,5,6,7)] <- NULL
cn2

#cowplot::plot_grid(dp1, dp2, labels=LETTERS[1:2])


