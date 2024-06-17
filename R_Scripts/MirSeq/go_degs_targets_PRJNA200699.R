library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(genefilter)
library(RColorBrewer)
library(qpdf)
library(ggplot2)
library(scatterpie)
library(dplyr)
hs <- org.Hs.eg.db
convert_to_ENTREZID <- function(genes){
  r=AnnotationDbi::select(hs,genes, columns = c("ENTREZID", "SYMBOL"),
                          keytype = "SYMBOL")
  return(r$ENTREZID)
}
mir_target_df_list <- list()
prj="PRJNA200699"
mir_target_df <- multimir_result[[prj]]@data %>% select(database,mature_mirna_id,target_symbol,target_entrez,target_ensembl,pubmed_id) %>%
  group_by(mature_mirna_id,target_entrez) %>% slice_head(n = 1) %>% ungroup() %>% 
  filter(target_symbol!="")
mir_target_df <- mir_target_df %>% group_by(target_symbol) %>% mutate(occ=n()) %>% ungroup()
target_occ5 <- unique((mir_target_df %>% filter(occ>5))$target_symbol)
PRJNA200699_target <- target_occ5
degs_PRJNA200696

degs_target <- list(targets=PRJNA200699_target,degs=degs_PRJNA200696$gene)
degs_target=lapply(degs_target,convert_to_ENTREZID)
cc_degs_target <- compareCluster(degs_target,fun = "enrichGO",
                      OrgDb         = org.Hs.eg.db,
                      keyType       ="ENTREZID",
                      readable      = TRUE)
cc_degs_target_simply <- simplify(cc_degs_target)
cc_degs_target_simply_df <- cc_degs_target_simply@compareClusterResult
View(cc_degs_target_simply_df)
duplicated(cc_degs_target_simply_df$ID)
df = cc_degs_target_simply_df
# Identify duplicated rows (based on all columns)
duplicates <- df[duplicated(df$ID) | duplicated(df$ID, fromLast = TRUE), ]

# Identify unique rows based on specific columns
uniques <- df[!(duplicated(df$ID) | duplicated(df$ID, fromLast = TRUE)), ]

# Combine duplicated rows first followed by unique rows
cc_degs_target_simply@compareClusterResult <- bind_rows(duplicates, uniques)
cc_degs_target_simply@compareClusterResult <- duplicates
View(cc_degs_target_simply@compareClusterResult)
dotplot(cc_degs_target_simply,showCategory=30)
ggsave("commonGO_targets_degs_PRJNA200699.jpg",scale=1.5)
cc_degs_target@compareClusterResult

targets_go <- enrichGO(gene       = degs_target[["targets"]],
                       OrgDb         = hs,
                       keyType       ="ENTREZID",
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
targets_go_df <- targets_go@result
degs_go <- enrichGO(gene       = degs_target[["degs"]],
                       OrgDb         = hs,
                       keyType       ="ENTREZID",
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)
degs_go_df <- degs_go@result
targets_degs_go_df <- targets_go_df %>% filter(ID %in% degs_go_df$ID)
colnames(targets_degs_go_df)
library(ggplot2)
ggplot() +
  geom_point(data = targets_degs_go_df, aes(x = Cluster, y = Description),fill=p.adjust) +
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
