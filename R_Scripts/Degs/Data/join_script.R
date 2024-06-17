library(DESeq2)
library(dplyr)
library(readxl)
library(tidyverse)
library(UpSetR)
library(gtools)
#find degs intersection

prj_list <- c("PRJNA200696","PRJNA216966","PRJNA649934")
degs <- read.csv(paste("DEG_",prj_list[1],".csv",sep=""))
intersection <- degs[[1]]
degs_dfs <- list()
for (prj in prj_list){
  degs <- read.csv(paste("DEG_",prj,".csv",sep=""))
  degs_dfs[[prj]] <- degs[,1]
  #degs_dfs[[prj]] <- degs %>% rename(gene=X) %>% select(gene,log2FoldChange,padj)
  #degs_dfs[[prj]] <- degs %>% rename(gene=X) %>% select(gene)
}

up_set <- upset(fromList(degs_dfs))
intersections <- up_set$New_data
temp <- unlist(degs_dfs, use.names = FALSE)
rownames(intersections) <- temp[ !duplicated(x1) ]
write.csv(intersections,"degs_intersection.csv")

perm <- permutations(3,2,prj_list)
perm <- as.data.frame(perm)
apply(perm,1,FUN = function(x) print(cbind(intersections[x[1]],intersections[x[2]])))


pdf("upset.pdf")
print(up_set)
dev.off()


common_elements <- Reduce(full_join, degs_dfs);

common_degs <- degs_dfs %>%
                  bind_rows(.id='df') %>%
                  filter(gene %in% unique(gene[duplicated(gene)])) %>%
                  mutate(df2=df) %>% 
                  spread(df,df2)
#prova
#"GABRA5" %in% degs_dfs[["PRJNA649934"]][,1]



join <- degs_dfs %>% reduce(inner_join, by='X')

join2 <- inner_join(degs_dfs[[prj_list[1]]],degs_dfs[[prj_list[2]]],by = "X",
                   suffix = c("", prj_list[2])) %>% 
                  inner_join(degs_dfs[[prj_list[3]]],by = "X",suffix=c(prj_list[1],prj_list[3]))
write_csv(join2,"DEG_intersection.csv")


gene_count_pubblic <- read_excel("PRJNA200696_mRNA_counts.xls", 
                                      sheet = "Filtered")

gene_count_mat <- read.csv("gene_count_matrix_PRJNA200696.csv")
gene_count_mat <- read.csv(paste("gene_count_matrix_PRJNA200696.csv"))
gene_count_mat$gene_id <- sapply(strsplit(gene_count_mat$gene_id,
                                               split = "|",fixed = TRUE),`[`, 1)


join_gene_count <- inner_join(gene_count_pubblic,gene_count_mat,by="gene_id",
                              suffix=c("pubblic","reanalyzed"))


