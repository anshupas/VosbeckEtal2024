# For recreating the UMAP, you can directly use the R Seurat object at:
# data/Lietal2023/lietalSeuratObj.RData
# Once you load the object into R, 
# replace the cellType slot in the Seurat object
# with the numbers in the 2nd column of the attached csv file after matching the cell Ids.

library(Seurat)
library(tidyverse)

load("data/Lietal2023/lietalSeuratObj.RData")

seo

DimPlot(seo, reduction = "umap")



df = read.table("data/Lietal2023/Hongzhe_clustering_and_2dUMAP_louvain11_merged_95_renamed.csv", 
                sep = ",", row.names = 1)

df$cell_id = rownames(df)

df_names = read_tsv("data/Lietal2023/hongzhe_celltype_names.tsv")
df_names = df_names %>% 
    select(cluster_id, cellType, Genes)

df = df %>% 
    left_join(df_names, by = c("louvain11_merged.95_renamed" = "cluster_id"))



df = df %>%
    arrange(match(cell_id, rownames(seo@meta.data)))

any(df$cell_id != rownames(seo@meta.data))

setdiff(df$cell_id, rownames(seo@meta.data))

seo@meta.data$cellType = as.character(df$cellType)


pdf("umap_with_hongzhe_clustering.pdf", width=12, height=9)
DimPlot(seo, reduction = "umap", group.by = "cellType", label = TRUE, repel = T)
dev.off()



df_anno = df %>% 
    select(cellType, louvain11_merged.95_renamed) %>% 
    distinct()



df %>% 
    mutate(louvain11_merged.95_renamed = factor(louvain11_merged.95_renamed, 
                                                levels = sort(unique(df$louvain11_merged.95_renamed)))) %>% 
    #mutate(louvain11_merged.95_renamed = as.character(louvain11_merged.95_renamed)) %>% 
    ggplot(mapping = aes(x = umap2d_X, y = umap2d_Y, color = louvain11_merged.95_renamed)) +
    geom_point(size = 0.1) +
    theme_classic() +
    xlab("UMAP_1") +
    ylab("UMAP_2") + 
    guides(color = guide_legend(nrow = 5,
                                override.aes = list(size = 2))) +
    scale_color_discrete(labels = setNames(paste(df_anno$louvain11_merged.95_renamed, 
                                                 df_anno$cellType), 
                                           df_anno$louvain11_merged.95_renamed)) +
    labs(color = "Hongzhe clustering") + 
    theme(legend.position="bottom",
          legend.text  = element_text(size = 8),
          legend.title=element_blank())
ggsave("umap_with_hongzhe_clustering.pdf", width = 14, height = 7)




msc_cells = rownames(seo@meta.data[seo@meta.data$orig.ident == "MSC",])  


df

df_anno = df %>% 
    #filter(cell_id %in% msc_cells) %>% 
    filter(cellType == "Stromal cells") %>% 
    select(cellType, louvain11_merged.95_renamed) %>% 
    distinct()



df %>% 
    #filter(cell_id %in% msc_cells) %>% 
    filter(cellType == "Stromal cells") %>% 
    mutate(louvain11_merged.95_renamed = factor(louvain11_merged.95_renamed, 
                                                levels = sort(unique(df$louvain11_merged.95_renamed)))) %>% 
    # mutate(louvain11_merged.95_renamed = as.character(louvain11_merged.95_renamed)) %>% 
    ggplot(mapping = aes(x = umap2d_X, y = umap2d_Y, color = louvain11_merged.95_renamed)) +
    geom_point(size = 0.1) +
    theme_classic() +
    xlab("UMAP_1") +
    ylab("UMAP_2") + 
    guides(color = guide_legend(nrow = 3,
                                override.aes = list(size = 2))) +
    scale_color_discrete(labels = setNames(paste(df_anno$louvain11_merged.95_renamed, 
                                                 df_anno$cellType), 
                                           df_anno$louvain11_merged.95_renamed)) +
    labs(color = "Hongzhe clustering") + 
    theme(legend.position="bottom",
          legend.text  = element_text(size = 8),
          legend.title=element_blank())
ggsave("umap_with_hongzhe_clustering_stromal.pdf", width = 5, height = 7)



## get percentage of NRP2 and NCAM1 positive cells in 
## 29 and/or 37

df_s = df[df$louvain11_merged.95_renamed %in% c(29, 39),]

counts = seo@assays$RNA@counts[c('Nrp2', 'Ncam1'), ]

counts = as.data.frame(counts)

counts$gene = rownames(counts)

counts = counts %>% 
    pivot_longer(cols = -gene, names_to = 'cell_id') %>% 
    pivot_wider(id_cols = c(cell_id), 
                names_from = gene, 
                values_from = value)


df_s
dim(df_s)
dim(counts)


df_s = df_s %>% 
    left_join(counts, by = "cell_id")


df_s %>% 
    group_by(louvain11_merged.95_renamed, cellType) %>% 
    summarise(no_nrp2 = sum(Nrp2 > 0), 
              no_ncam1 = sum(Ncam1 > 0),
              no_both = sum(Ncam1 > 0 & Nrp2 > 0),
              no_cells = n())




df = df %>% 
    left_join(counts, by = "cell_id")


df_genes = df %>% 
    group_by(louvain11_merged.95_renamed, cellType) %>% 
    summarise(no_nrp2 = sum(Nrp2 > 0), 
              no_ncam1 = sum(Ncam1 > 0),
              no_both = sum(Ncam1 > 0 & Nrp2 > 0),
              no_cells = n())



df_genes$percent_nrp2 = round(df_genes$no_nrp2/df_genes$no_cells,3)
df_genes$percent_ncam1 = round(df_genes$no_ncam1/df_genes$no_cells,3)
df_genes$percent_both = round(df_genes$no_both/df_genes$no_cells,3)

df_genes %>% 
    arrange(-percent_nrp2) %>% 
    rename(cluster_id = louvain11_merged.95_renamed) %>% 
    write_tsv("percent_nrp2_and_ncam1_in_hongzhe_clustering.tsv")
















