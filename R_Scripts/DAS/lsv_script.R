library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
#library(qpdf)
library(genefilter)
library(ComplexHeatmap)
library(RColorBrewer)
library(igraph)
library(qgraph)
library(colourvalues)
library(readr)
library(data.table)
library(splitstackshape)
library(janitor)
library(stringr)
library(pheatmap)
library(dplyr)
library(tidyverse)

parse_voila_delta_tsv = function(file_path){
  ####this function takes the file path to a voila tsv and parses it using splitstackshape and then combines the lsv_id and junction coordinates to make a unique
  ###lsv_id_junction for each lsv to get a PSI per splice event in long format
  #read in the file using data table and janitor to get clean names
  if(file.exists(file_path)){
    voila_org =  as.data.table(clean_names(fread(file_path)))
    # look for columns with conditional psi
    condition_psi_cols = colnames(voila_org)[grep("_mean_psi",colnames(voila_org))]
    
    
    columns_to_split = c("mean_dpsi_per_lsv_junction",
                         "probability_changing",
                         "probability_non_changing",
                         "junctions_coords",
                         "de_novo_junctions",
                         "ir_coords",condition_psi_cols)
    
    # for some reason....some of this columns didn't show up with me running majiq on the old bam files, different GFF or genome or something is possible route cuase sofor now
    voila_melt = cSplit(voila_org, splitCols = columns_to_split, sep = ';', direction = "long")
    #the splitting and melting can produce non unique rows, so remove those
    voila_melt = unique(voila_melt)
    
    #also can result in some rows that have NA in their psi columns
    voila_melt = voila_melt[!is.na(probability_changing)]
    #now to make a unique identifier cmbined the lsv_id and junction coordinates
    voila_melt[,lsv_junc := paste(lsv_id, junctions_coords,sep = "_")]
    # split the coords to have a start and an exon end
    # same with the junctions and the retained introns
    voila_melt[,c("junc_start","junc_end") := tstrsplit(junctions_coords, "-")]
    voila_melt[,c("ir_start","ir_end") := tstrsplit(ir_coords, "-")]
    voila_melt[,paste_into_igv_junction := paste0(seqid, ":",junc_start, "-",junc_end)]
    # I find this column irritating
    voila_melt$ucsc_lsv_link = NULL
    # remove these as they're now redundant
    voila_melt$junctions_coords = NULL
    voila_melt$ir_coords = NULL
    voila_melt$exons_coords = NULL
    # helpful columns
    voila_melt[,junc_dist := abs(as.numeric(junc_start) - as.numeric(junc_end))]
    
    
    
    return(voila_melt)
  } else{
    print(file_path)
    return("File not found")
  }
}

print(getwd())

prj_list <- c("PRJNA200696","PRJNA216966","PRJNA649934")
human_tf <- read.delim("_TF.txt",sep = "\t")

degs <- list()
lsv <- list()
lsv_filtered <- list()
lsv_tf <- list()
lsv_all <- list()
nrow(lsv_all$PRJNA200696)
#check if gene is lsv
lsv_all$PRJNA649934[lsv_all$PRJNA649934$gene_name=="BCAS2",]
colnames(lsv_all$PRJNA649934)

colnames(lsv_all$PRJNA216966)
for (prj in prj_list){
  degs[[prj]] <- read.csv(paste("../DEG_",prj,".csv",sep=""),row.names = 1)
  #carico lsv dei geni con almeno un dpsi per junction >|0.2|
  ##lsv[[prj]] <- read_tsv(paste("data/",prj,"_voila_filtered.tsv",sep = ""),comment = '#',skip_empty_rows = T)
  lsv[[prj]] <- as_tibble(parse_voila_delta_tsv(paste("data/",prj,"_voila.tsv",sep = "")))
  lsv_all[[prj]] <- lsv[[prj]]
  lsv[[prj]] <- lsv[[prj]] %>% filter(abs(mean_dpsi_per_lsv_junction) > 0.2)
  #filtro per lsv di TF
  lsv_tf[[prj]] <- lsv[[prj]] %>% filter(gene_name %in% human_tf$Symbol)
  #filtro i geni con lsv prendendo quelli significativamente espressi differenzialmente
  lsv_filtered[[prj]] <- lsv[[prj]][lsv[[prj]]$gene_name %in% rownames(degs[[prj]]),]
  #adding FC column
  lsv_filtered[[prj]]["FC"] <- sapply(lsv_filtered[[prj]]$gene_name,FUN = function(x) degs[[prj]][rownames(degs[[prj]])==x,]$log2FoldChange)
  #reordering columns in df
  lsv_filtered[[prj]] <- lsv_filtered[[prj]][, colnames(lsv_filtered[[prj]])[c(1:2,ncol(lsv_filtered[[prj]]),3:ncol(lsv_filtered[[prj]])-1)]]
}

View(lsv_all$PRJNA200696)

library(readxl)
library(openxlsx)
write.xlsx(lsv_all, "DAS_all_bioprojects.xlsx")

heatmaps_from_list <- function(list,tipo_lista){
  for(prj in names(list)){
    heat_matrix <- list[[prj]] %>% select(c(colnames(list[[prj]])[grep("_mean_psi",colnames(list[[prj]]))],"lsv_junc")) 
    heat_matrix <- as.matrix(heat_matrix)
    lsv_junc <- heat_matrix[,3]
    #cambia l'lsv con il gene name
    #lsv_junc <- sapply(lsv_junc, function(x){str_replace(x,pattern = (str_split(x,pattern = ":")[[1]][1]),list[[prj]][list[[prj]]$gene_id==str_split(x,pattern = ":")[[1]][1],]$gene_name[1])})
    heat_matrix <- heat_matrix[,1:2]
    heat_matrix <- apply(heat_matrix,2,as.numeric)
    rownames(heat_matrix) <- lsv_junc
    heat_map <- pheatmap(heat_matrix,
                         #annotation_col = ann_col,
                         #annotation_colors = ann_colors,
                         angle_col = "0",
                         main = paste(prj,"LSV",tipo_lista,sep=" "),
                         fontsize = 4,
                         fontsize_col = 4,
                         fontsize_row = 4)

    pdf(file =paste(prj,"_heatmap_",tipo_lista,".pdf",sep = ""),height = 15)
    print(heat_map)
    dev.off()
  }
}

heatmaps_from_list(lsv,"ALL")
heatmaps_from_list(lsv_tf,"TF")
heatmaps_from_list(lsv_filtered,"DEG")

#lsv contiene i DAS con dPSI>0.2 (dPSI tra sample di stesso dataset)
#join dei due progetti per vedere lsv in comune tra MURAL e CumulusCell in diverse fasi di maturazione
join_lsv <- inner_join(x = lsv$PRJNA200696,y = lsv$PRJNA216966,by="gene_name",suffix=c("_20","_21"))
join_lsv <- join_lsv %>% select(c("gene_name","mean_dpsi_per_lsv_junction_20","lsv_id_20","mean_dpsi_per_lsv_junction_21","lsv_id_21","lsv_junc_20","lsv_junc_21"))
#filtro per gli lsv comnuni, applicando doppio controllo su junction e lsv_id
lsv_common <- join_lsv %>%filter(lsv_junc_20==lsv_junc_21,lsv_id_20==lsv_id_21)
View(lsv_common)

#heatmap join
join_mat <- inner_join(x = lsv$PRJNA200696,y = lsv$PRJNA216966,by="lsv_junc")
join_mat <- join_mat %>% select(c("lsv_junc",colnames(join_mat)[grep("_mean_psi",colnames(join_mat))]))
join_mat <- as.matrix(join_mat)
tmp_rn <- join_mat[,1]
tmp_rn <- sapply(tmp_rn, function(x){str_replace(x,pattern = (str_split(x,pattern = ":")[[1]][1]),lsv$PRJNA216966[lsv$PRJNA216966$gene_id==str_split(x,pattern = ":")[[1]][1],]$gene_name[1])})
join_mat <- join_mat[,-1]
join_mat <- apply(join_mat,2,as.numeric)
rownames(join_mat) <- tmp_rn
pdf(file =paste("JOIN_heatmap_.pdf",sep = ""))
pheatmap(join_mat,
             #annotation_col = ann_col,
             #annotation_colors = ann_colors,
             angle_col = "0",
             main = "Join HEATMAP",
             fontsize = 4,
             fontsize_col = 4,
             fontsize_row = 4)
dev.off()





gene_counts <- list()
gene_counts_filtered <- list()
for (prj in prj_list){
  gene_count_matrix <- read.csv(paste("../gene_count_matrix_",prj,".csv",sep = ""), row.names=1)
  rownames(gene_count_matrix) <- sapply(strsplit(rownames(gene_count_matrix),
                                                 split = "|",fixed = TRUE),`[`, 1)
  gene_counts[[prj]] <- gene_count_matrix
  gene_counts_filtered[[prj]] <- gene_counts[[prj]][rownames(gene_counts[[prj]]) %in% lsv_filtered[[prj]]$gene_name,]
}





parse_try <- parse_voila_delta_tsv("data/PRJNA649934_voila_filtered.tsv")
parse_try <- parse_try %>% filter(gene_id %in% human_tf$Ensembl, abs(mean_dpsi_per_lsv_junction) > 0.2)


View(read_tsv("data/het_tsv_voila.tsv",comment = '#',skip_empty_rows = T))

