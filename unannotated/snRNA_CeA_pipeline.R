---
title: "PARC HP/LP snMultiome scratch work"
author: "Justin Anderson"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
    code_folding: hide
    df_print: kable
    toc: true
    toc_depth: 3
    number_sections: true
    highlight: tango
editor_options: 
  chunk_output_type: console
---
This is where I am prototyping a workflow/pipeline for the single nuclei multome data for the CeA in the HP/LP mice as part of the PARC P001. The ultimate goal here is to recreate something similar to what Kip/Rita made, in order to verify those results and finish the paper ASAP.

This is the final version of the pipeline, where we are not longer removed doublets called on the cell level

#Setup

Setup some default options and methods

## Defaults
```{r setup, include = FALSE}
knitr::opts_chunk$set(error = 0,
                      warning = FALSE,
                      tidy = TRUE)
# Default kable table
kableDefault <- function(table, ...) {
  knitr::kable(table, ...) %>%
    kable_styling(full_width = T,
                  fixed_thead = T,
                  bootstrap_options = c("hover","striped"))
}

theme_JQA <- function(x) {
  theme_bw(base_size = 12) +
    theme(panel.border = element_blank(),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
}

#40GB for the workers (we will only use 4).
options(future.globals.maxSize = 10 * 1024 ^ 3)

#rng
saved_seed <- 23587
set.seed(saved_seed)
```

Load in the necessary packages

## Libraries
```{r}
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/groups/rosalind/andejust/projects/snRNA/Rdata")

#Tidyverse
pacman::p_load(
  plyr,
  tidyverse,
  dplyr,
  kableExtra,
  ggplot2
  )

#Single Cell
pacman::p_load(
  Signac,
  Seurat
)

#Futurevese
pacman::p_load(future,
               future.callr)

+#Gene Annotation
pacman::p_load(
  biomaRt,
  org.Mm.eg.db
  )

#Visualization
pacman::p_load(
  ggpubr,
  ggstatsplot,
  patchwork,
  colorspace,
  ComplexHeatmap,
  viridis,
  circlize,
  clustree
)

#Export
pacman::p_load(
  openxlsx,
  readxl
)
```

## Data
We are doing the simple example first, just working with the snRNA data
https://satijalab.org/seurat/articles/pbmc3k_tutorial

https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html suggests we start with the raw filtered matrix or that we actually import every sample ourselves (so that we can do our own filtering)

### Sample info
```{r}
sampleInfo <- read_delim("/home/groups/parc/primeseq/TamaraPR/20/@files/cellranger_analysis/sample.names",
                       col_names = F)
sampleInfo <- sampleInfo[,c(3,5)]

colnames(sampleInfo) <- c("name", "id")

sampleInfo$genotype <- str_sub(sampleInfo$name, 1,1)
sampleInfo$sex <- str_sub(sampleInfo$name, 5,5)
sampleInfo$name <- str_sub(sampleInfo$name, 1,5)

#we can add in the lane manually.. it lvies at the above path.. 20/@files/raw/Reports/html/HF52FDSX3/MRL211207TP/all/all/laneBarcode.html

sampleInfo <- sampleInfo %>%
  mutate(GEX_lane = if_else(str_detect(sampleInfo$name, "100"), 2, 3))
```

### 10x Samples
```{r}
#~13 seconds per sample

#This is the OHSU data, we are switching to Kip's which is aligned to the newer mouse reference genome
# raw_data <- list()
# for(j in 1:16) {
#   cat(sprintf("Working on sample HS_%d.\n", j))
#   tmp_time <- Sys.time()
#   raw_data[[paste0("HS_",j)]] <- Read10X_h5(
#     paste0("/home/groups/rosalind/andejust/projects/snRNA/OHSU/HS_",
#            j,
#            "/outs/raw_feature_bc_matrix.h5")
#   )
#   cat(sprintf("    Sample HS_%d is done (%s elapsed).\n", j, format(Sys.time()-tmp_time)))
# }
# rm(tmp_time,j)
# gc()


raw_data <- list()

for(j in sampleInfo$name) {
  cat(sprintf("Working on sample %s.\n", j))
  tmp_time <- Sys.time()
  raw_data[[j]] <- Read10X_h5(
    paste0("/home/groups/rosalind/andejust/projects/snRNA/WAKE/Cell_Ranger_Sample_Outputs/",
           j,
           "_outs/outs/raw_feature_bc_matrix.h5")
  )
  cat(sprintf("    Sample %s is done (%s elapsed).\n", j, format(Sys.time()-tmp_time)))
}
rm(tmp_time,j)
gc()
```

We set a low cutoff for the first pass
```{r}
#~30 seconds per sample
GE_data <- list()
for(j in sampleInfo$name) {
  cat(sprintf("Working on sample %s.\n", j))
  tmp_time <- Sys.time()
  tmp_GE <- CreateSeuratObject(
    counts = raw_data[[j]]$`Gene Expression`,
    project = j)
  
  #Drop cells with low numbers of genes early here
  GE_data[[j]] <- subset(x = tmp_GE, 
                         subset = 
                           (nCount_RNA > 250))
  
  rm(tmp_GE)
  cat(sprintf("    Sample %s is done (%s elapsed).\n", j, format(Sys.time()-tmp_time)))
}
rm(j)
gc()
```

Combine into one object - This takes awhile (mostly the join layers one)
```{r}
seurat_GE <- merge(GE_data[[1]],
              GE_data[2:16],
              add.cell.ids = names(GE_data))

seurat_GE <- JoinLayers(seurat_GE)

#We can get rid of the raw data now
rm(raw_data)
```

#### Meta QC
We are going to build up the meta data file and work through some more formal QC, compared to what we did above.

```{r}
GE_meta <- seurat_GE@meta.data

#Mitochondrial contamination - should be ~0 for snRNA data. Wang et al. use 0.05%
GE_meta$percent.mt <- PercentageFeatureSet(seurat_GE, 
                                           pattern="^mt-")

GE_meta$percent.rp <- PercentageFeatureSet(seurat_GE, 
                                           pattern="^Rp[ls]")
```

Babraham's lab also investigates the highest expressed gene in each cell.

```{r}
#There isn't a method for which.max for sparse matrices, it forces it to dense which takes like
#20 gigs of allocation and also forever.  This is as fast at small numbers of columns, and respects
#sparse memory allocation. Because of the allocation issue it is MUCH faster at large numbers of
#columns (10-30x).

sparseWhich.max <- function(m) {
  #allocate the matrix, 2 by ncols
  tmp_data <- array(NA, ncol(m))
  for(j in 1:ncol(m)) {
    if((m@p[j+1]-m@p[j])>0)
    tmp_data[j] <- m@p[j]+which.max(m@x[(m@p[j]+1):m@p[j+1]])
  }
  
  return(data_frame(
    max = m@x[tmp_data],
    gene = rownames(m)[m@i[tmp_data]+1],
    cell = colnames(m)
  ))
}


#First lets verify that this new method works
# test2 <- as.matrix(counts)
# 
# 
# test2 <- data_frame(max = unname(sparseMatrixStats::colMaxs(counts)),
#                     gene = rownames(counts)[apply(counts,2,which.max)],
#                     cell = colnames(counts))
# test2$max <- unname(test2$max)
# all.equal(test1,test2)
#Hoorah

# test <- counts[,1:2000]
# microbenchmark::microbenchmark(
#   sparseWhich(test),
#   data_frame(max = unname(sparseMatrixStats::colMaxs(test)),
#                     gene = rownames(test)[apply(test,2,which.max)],
#                     cell = colnames(test)),
#   times = 10,
#   check = 'equal'
# )
#This is already 20x faster at 2000 columns
#rm(test,test1,test2,test3)
```

Okay lets lookinto the highest expressed gene for a bit
```{r}
counts <- GetAssayData(object = seurat_GE, layer = "counts")
#Rows = genes, columns = cells

tmp_HGE <- sparseWhich.max(counts)


tmp_HGE %>%
  dplyr::count(gene) %>%
  mutate(prop = n/sum(n)) %>%
  arrange(desc(n))

#Okay, as they said MALAT1 is the highest expressed gene in 83% of cells, so we need to drop it.

tmp_HGE <- sparseWhich.max(counts[-which(rownames(counts) == "Malat1"),])
tmp_HGE %>%
  dplyr::count(gene) %>%
  mutate(prop = n/sum(n)) %>%
  arrange(desc(n))
```

```{r}
tmp_cuttoffs <- read_xlsx("/home/groups/rosalind/andejust/projects/snRNA/WAKE/snRNA_cutoffs.xlsx") %>%
  mutate(across(starts_with("n"), function(x) {round(10^x)}))
```

```{r}
GE_meta$cells <- rownames(GE_meta)
GE_meta$sample <- GE_meta$orig.ident

GE_meta <- left_join(GE_meta,
                  sampleInfo,
                  by = c("sample" = "name"))

GE_meta <- left_join(GE_meta,
                     tmp_HGE,
                     by = c("cells" = "cell"))

GE_meta <- GE_meta %>%
  mutate(percent.hge = max/nCount_RNA) 

GE_meta <- left_join(GE_meta,
                  tmp_cuttoffs) %>%
  dplyr::mutate(class = dplyr::case_when(
    ((nCount_RNA > nCount_split) & (nFeature_RNA > nFeature_split)) ~ "A",
    ((nCount_RNA > nCount_low) & (nFeature_RNA > nFeature_low)) ~ "B",
    .default = "C"
  ))
                     
rm(tmp_cuttoffs,tmp_HGE)
```

Apply the cutoffs
```{r}
keep <- GE_meta %>%
  filter(percent.mt < 0.5,
         percent.hge < 10/100,
         log10(nFeature_RNA)/log10(nCount_RNA) > complexity,
         nCount_RNA > nCount_low,
         nFeature_RNA > nFeature_low,
         nCount_RNA < nCount_high,
         nFeature_RNA < nFeature_high)

GE_meta$keep = GE_meta$cells %in% keep$cells
rownames(GE_meta) <- GE_meta$cells

table(GE_meta$keep)

seurat_GE@meta.data <- GE_meta

rm(keep)
```

### Subset/Filtering Cells
```{r}
filtered_GE <- subset(x = seurat_GE,
                      subset = (keep==T))

filtered_GE
```

### Filtering Genes

Gene level filtering... we dont want genes only expressed in a tiny number of cells.. we have 69318  cells.. and 33589 genes. Lets say we have a gene to have more than 2 counts in over 10 cells (ucdavis criteria)

```{r}
counts <- GetAssayData(object = filtered_GE, layer = "counts")

nonzeros <- Matrix::rowSums(counts>2)

nonzeros <- data_frame(gene = rownames(counts),
                       ncells = nonzeros)

nonzeros <- nonzeros %>%
arrange(ncells) %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) 

sum(nonzeros$ncells > 10)

nonzeros %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) %>%
  ggplot() +
  geom_point(aes(y = idx, x = ncells)) +
  theme_JQA() +
    coord_cartesian(xlim=c(0,10000))

nonzeros[c(
  which.max(nonzeros$idx>10),
  which.max(nonzeros$idx>20),
  which.max(nonzeros$idx>30),
  which.max(nonzeros$idx>40),
  which.max(nonzeros$idx>50),
  which.max(nonzeros$idx>60),
  which.max(nonzeros$idx>70),
  which.max(nonzeros$idx>80),
  which.max(nonzeros$idx>90)), ]


#This leaves us with ~ 14k genes which is great  
to_drop <- nonzeros$ncells < 10

sum(!to_drop)
#that leaves 14k genes, which is fine I think
rownames(filtered_GE@meta.data)

filtered_GE <- subset(x = filtered_GE,
                      features = nonzeros$gene[!to_drop])
filtered_GE 
rm(counts,nonzeros,to_drop)
```

Final is 69318 cells and 14582 genes

### Dropping IEGs and Sex Genes
```{r}
counts <- GetAssayData(object = filtered_GE, layer = "counts")

to_drop <- rownames(counts) %in% c("Btg2", "Jun", "Egr4", "Fosb",
                                   "Junb", "Gadd45g", "Fos", "Arc",
                                   "Nr4a1", "Npas4", "Coq10b", "Tns1", "Per2",
                                   "Ptgs2", "Rnd3", "Tnfaip6", "Srxn1",
                                   "Tiparp", "Ccnl1", "Mcl1", "Dnajb5", "Nr4a3",
                                   "Fosl2", "Nptx2", "Rasl11a", "Mest", "Sertad1",
                                   "Egr2", "Midn", "Gadd45b", "Dusp6", "Irs2", "Plat",
                                   "Ier2", "Rrad", "Tpbg", "Csrnp1", "Peli1", "Per1",
                                   "Kdm6b", "Inhba", "Plk2", "Ifrd1", "Baz1a", "Trib1",
                                   "Pim3", "Lrrk2", "Dusp1", "Cdkn1a", "Pim1", "Sik1",
                                   "Frat2", "Dusp5", "Xist", "Tsix", "Eif2s3y", "Ddx3y",
                                   "Uty", "Kdm5d")

filtered_GE <- subset(x = filtered_GE,
                      features = rownames(counts)[!to_drop])
filtered_GE 
rm(counts, to_drop)
```
14533 features across 69318 samples within 1 assay 

#Overall Analysis

```{r}


rm(seurat_GE, GE_data, GE_meta)
```

```{r}
tmp_file = "filtered_seurat.RData"

if(!file.exists(tmp_file)) {
  save(filtered_GE, file = tmp_file)
} else {
  load(tmp_file)
}
rm(tmp_file)
```

## Calling Clusters
We are dropping HS16 and HS3

### Normalization
```{r}
lognormal_GE <- SplitObject(filtered_GE, split.by = "id")
# lognormal_GE[["HS_16"]] <- NULL
# lognormal_GE[["HS_3"]] <- NULL
```
Dropping those two samples drops us to 53839 samples


```{r}
plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)
lognormal_GE <- future.apply::future_lapply(lognormal_GE,
                                      function(x) {
                                        SCTransform(x,
                                                    assay = "RNA",
                                                    vars.to.regress = "percent.mt",
                                                    verbose = FALSE,
                                                    seed.use = saved_seed,
                                                    )
                                      },
                                      future.seed = saved_seed
) 

lognormal_GE <- future.apply::future_lapply(lognormal_GE,
                                            function(x) {
                                              NormalizeData(x,
                                                            assay = "RNA",
                                              )
                                            },
                                            future.seed = saved_seed
) 

lognormal_GE <- merge(lognormal_GE[[1]],
              lognormal_GE[2:14],
              add.cell.ids = names(lognormal_GE))

lognormal_GE[["RNA"]] <- JoinLayers(lognormal_GE[["RNA"]])

VariableFeatures(lognormal_GE, assay = "SCT") <- rownames(lognormal_GE[["SCT"]]@scale.data)

# lognormal_GE <- NormalizeData(subset(filtered_GE,
#                                      subset = (id != "HS_3") &
#                                        (id != "HS_16")))
# 
# lognormal_GE <- FindVariableFeatures(lognormal_GE)
# 
# lognormal_GE <- ScaleData(lognormal_GE)

plan(sequential)
gc()
```


```{r}
lognormal_GE <- RunPCA(lognormal_GE,
                       assay = "SCT")
```

```{r}
ElbowPlot(lognormal_GE, ndims = 40)
```

```{r}
DimHeatmap(lognormal_GE,dims=1:10, cells=500)
DimHeatmap(lognormal_GE,dims=11:20, cells=500)
DimHeatmap(lognormal_GE,dims=21:30, cells=500)
DimHeatmap(lognormal_GE,dims=31:40, cells=500)
```


```{r}
lognormal_GE <- RunUMAP(lognormal_GE, 
                        dims = 1:40,
                        reduction = "pca",
                        reduction.name = "umap",
                        seed.use = saved_seed)
```

### Clustering
```{r}
lognormal_GE <- FindNeighbors(lognormal_GE,
                              dims = 1:40,
                              reduction = "pca",
                              graph.name = paste0("pca_", c("nn","snn")))
```

```{r}
library(reticulate)
numpy <- reticulate::import("numpy")
pandas <- reticulate::import("pandas")
leidenalg <- reticulate::import("leidenalg")

plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)

lognormal_GE <- FindClusters(lognormal_GE,
                             method  = "igraph",
                             algorithm = 4,
                             random.seed = saved_seed,
                             graph.name = "pca_snn",
                             resolution = c(0.01,0.05,0.1,0.25,0.5,0.8,1,1.5,3,4.5,6,10))

plan(sequential)
gc()
```

#### Save/load normalized 
```{r}
tmp_file = "lognormal_seurat.RData"


if(!file.exists(tmp_file)) {
  save(lognormal_GE, file = tmp_file)
} else {
  load(tmp_file)
}

rm(tmp_file)
```


```{r}
clustree(lognormal_GE, prefix = "pca_snn_res.")
```

```{r}
DimPlot(lognormal_GE,
        reduction = "umap",
        group.by = c("genotype",
                     "id",
                     "class",
                     "GEX_lane",
                     "pca_snn_res.0.01",
                     "pca_snn_res.0.05",
                     "pca_snn_res.0.1",
                     "pca_snn_res.0.25",
                     "pca_snn_res.0.5"),
        label = T)
```

```{r}
plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)
lognormal_GE <- PrepSCTFindMarkers(lognormal_GE)
plan(sequential)
gc()

#Found 14 SCT models. Recorrecting SCT counts using minimum median counts: 5314
```

As before, we will move ahead with 0.5 as a good intermediate choice.. we also do not really care about overclustering here, its just about finding populations of neuronal cells, particlarly GABAergic ones.

#### Cluster Markers

```{r}
tmp_filter <- list()
tmp_cluster <- list()
tmp_filter_id <- list()
```

```{r}
tmp_gaba <- c("Gad1", "Gad2", "Slc32a1")
tmp_glut <- c("Slc17a7", "Slc17a6")
tmp_neuronal <- c(tmp_gaba, tmp_glut)

#what if we use the markers from Rita, which is driven by this data
tmp_non <- list("Gja1_Aqp4_Fgfr3_Aldoc_Mfge8",
                "C1qc_C1qa_C1qb",
                "Vcan_Apoe_Olig2_Cdh19_Sox10",
                "Olig1_C1ql1_Plp1_Mbp_Mog_Mobp_Mag_Cnp_Cldn11") %>%
  sapply(function(x) str_split(x,"_")) %>%
  unlist() %>% unique()

test <- lognormal_GE@assays[["SCT"]]@counts@Dimnames[[1]]

#tmp_non <- unique(tmp)
#Add in some markers given by Wang and Hochgerner
tmp_non <- c(tmp_non, "Ctss", "Opalin", "Pdgfra", "Enpp2", "Otx2os1") %>% unique()
tmp_neuronal <- tmp_neuronal[tmp_neuronal %in% (test)]
tmp_non <- tmp_non[tmp_non %in% (test)]

rm(test)
```

```{r}
for( i in paste0("pca_snn_res.",c (0.01,0.05,0.1,0.25,0.5,0.8,1,1.5,3,4.5,6,10)) )
{
  
 test <- AggregateExpression(lognormal_GE,
                          return.seurat = T,
                          assay = "SCT",
                          group.by= i)$SCT$count

test[test==0] <- NA
  
test2 <- data_frame(
  neuronal = colMeans(test[tmp_neuronal, ], na.rm = T),
  non = colMeans(test[tmp_non, ], na.rm = T),
  gaba = colMeans(test[tmp_gaba, ], na.rm = T),
  glut = colMeans(test[tmp_glut, ], na.rm = T),
  seurat_cluster = 1:(ncol(test))) %>%
  mutate(across(where(is.numeric), function(x) replace_na(x, 0)))

#Hochberg called doublet clusters by looking for ratios greater than 2..

#Clusters that are both neuronal and non
test2 <- test2 %>%
  dplyr::mutate(
    cluster_class = case_when(
      (non > glut & non > gaba) ~ "Non",
      gaba > glut ~ "Gaba",
      .default = "Glut",
    ),
    cluster_class = if_else(non+gaba+glut == 0, "Other", cluster_class)
  )

test2 %>%
  dplyr::count(cluster_class)
#This seems about right..

test2 <- test2 %>%
  rowwise() %>%
  dplyr::mutate(
    worst_cluster_ratio = case_match(
      cluster_class,
      "Gaba" ~ min(gaba/glut, gaba/non),
      "Glut" ~ min(glut/gaba, glut/non),
      "Non" ~ min(non/gaba, non/glut)
    )) %>%
  dplyr::mutate(
    doublet = worst_cluster_ratio < 3
  ) %>% ungroup()
  
  test3 <- left_join(test2 %>%
                       mutate(seurat_cluster = factor(seurat_cluster)),
                     lognormal_GE@meta.data %>%
                       dplyr::select(-starts_with(colnames(test2))),
                     by = c("seurat_cluster" = i)) %>%
    dplyr::filter(doublet == F) %>%
    group_by(id, cluster_class) %>%
    summarize(n = n()) %>%
    ungroup() %>% 
    group_by(id) %>%
    mutate(n = n / sum(n))
  
  tmp_filter_id[[i]] <- test3
  
  tmp_cluster[[i]] <- test2
  
  test2 <- left_join(test2 %>%
                       mutate(seurat_cluster = factor(seurat_cluster)),
                     lognormal_GE@meta.data %>%
                       group_by(!!sym(i), id) %>%
                       summarize(n = n(),
                                 cells = list(cells)) %>%
                       mutate(prop = n / sum(n),
                              n=sum(n),
                              cells = list(unique(unlist(cells)))) %>%
                       slice_max(prop),
                     by = c("seurat_cluster" = i))
  
  
  tmp_filter[[i]] <- test2
}
```

##### Cluster figure
```{r}
test <- lapply(tmp_filter, 
       function(x) {
         x %>%
           dplyr::mutate(cluster_class = if_else(doublet==F, cluster_class, "Doublet")) %>%
           group_by(cluster_class) %>%
           summarize(n = sum(n)) 
       }
)

(test %>% bind_rows(.id = "res") %>%
  ggplot(aes(x = cluster_class,
             y = n,
             group = cluster_class,
             color = res)) +
  geom_violin(scale = "width",
              show.legend = F) +
  geom_boxplot(show.legend = F,
               width = 0.25) +
  ggforce::geom_sina(scale = "width") +
  theme_JQA() +
  scale_color_discrete_diverging()) %>% plotly::ggplotly()


test %>% bind_rows(.id = "res") %>%
    dplyr::mutate(res = (res == "pca_snn_res.0.8")) %>%
  ggplot(aes(x = cluster_class,
             y = n,
             group = cluster_class,
             color = res)) +
  geom_violin(scale = "width",
              show.legend = F) +
  geom_boxplot(show.legend = F,
               width = 0.25) +
  ggforce::geom_sina(scale = "width") +
  theme_JQA() 
```
We are going to use 0.8, it is in the middle for all of em

##### Sample figure
```{r}
tmp_filter_id %>% bind_rows(.id = "res") %>%
  dplyr::mutate(res = (res == "pca_snn_res.0.8")) %>%
  ggplot(aes(x = cluster_class,
             y = n,
             group = cluster_class,
             color = res)) +
  geom_violin(scale = "width",
              show.legend = F) +
  geom_boxplot(show.legend = F,
               width = 0.25) +
  ggforce::geom_sina(scale = "width") +
  theme_JQA() +
  facet_wrap(~id, ncol = 4)
```


```{r}
colnames(lognormal_GE@meta.data)

test2 <- tmp_cluster[["pca_snn_res.0.8"]]

test2 <- left_join(lognormal_GE@meta.data %>%
                     dplyr::select(-starts_with(colnames(test2))),
                   tmp_cluster[["pca_snn_res.0.8"]] %>%
                     mutate(seurat_cluster = factor(seurat_cluster)),
                   by = c("pca_snn_res.0.8" = "seurat_cluster"))

rownames(test2) <- lognormal_GE[["SCT"]]@data@Dimnames[[2]]
test2$cells_old <- test2$cells
test2$cells <- rownames(test2)

lognormal_GE@meta.data <- test2

rm(test, test2, test3, tmp_cluster, tmp_filter, tmp_filter_id)
```

##### Non-Neuronals
We used non-neuronal markers given by Rita and found in Hochgerner and Wang.. but we know the nonneuronal populations are going to be more heterogeneous than the neuronal clusters, so maybe it would be better to look at each of those clusters individually.. so we can build better lists of markers and do a better job of doublet calling, maybe?

```{r}
Idents(object = lognormal_GE) <- "pca_snn_res.0.8"
tmp_meta <- lognormal_GE@meta.data

non_clusters <- tmp_meta %>%
  dplyr::filter(cluster_class != "Gaba",
                cluster_class != "Glut") %>%
  pull(pca_snn_res.0.8) %>% unique %>%
  unfactor() %>% as.numeric

neuron_clusters <- tmp_meta %>%
  dplyr::filter(cluster_class == "Gaba" |
                cluster_class == "Glut") %>%
  pull(pca_snn_res.0.8) %>% unique %>%
  unfactor %>% as.numeric()

gaba_clusters <- tmp_meta %>%
  dplyr::filter(cluster_class == "Gaba") %>%
  pull(pca_snn_res.0.8) %>% unique  %>%
  unfactor %>% as.numeric()

glut_clusters <- tmp_meta %>%
  dplyr::filter(cluster_class == "Glut") %>%
  pull(pca_snn_res.0.8) %>% unique  %>%
  unfactor %>% as.numeric()

plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)

tmp_markers <- future.apply::future_lapply(non_clusters,
                             function(i)
                             {
                               SeuratWrappers::RunPresto(lognormal_GE,
                                                         assay = "SCT",
                                                         ident.1 = i,
                                                         ident.2 = neuron_clusters)
                             },
                             future.seed = saved_seed
                             )

tmp_markers  <- c(tmp_markers,
                  future.apply::future_lapply(list(gaba = gaba_clusters,
                                                  glut = glut_clusters,
                                                   non = non_clusters),
                                              function(i) {
                                                SeuratWrappers::RunPresto(lognormal_GE,
                                                                          assay = "SCT",
                                                                          ident.1 = i)
                                              },
                                              future.seed = saved_seed
                  )
)
plan(sequential)
gc()        

# tmp_markers[["28"]] <- SeuratWrappers::RunPresto(lognormal_GE,
#                                                  assay = "SCT",
#                                                  ident.1 = 28,
#                                                  random.seed = saved_seed)

tmp_markers <- lapply(tmp_markers, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})


names(tmp_markers) <- c(non_clusters, "gaba", "glut", "non")
```

```{r}
overall_markers <- tmp_markers
```

##### DotPlot

Okay, now that we have markers for our stuff, can we find n=3-5 for each type so we can have an even footing for cell-type/cluster-type annotations (as apposed to above where we just grabbed a bunch of non-neuronal things)

```{r}
tmp <- overall_markers %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

#the number of markers we want per cluster
tmp2 <- tmp %>%
  #mutate(cluster = if_else(cluster %in% c("6", "9", "10"), "6", cluster)) %>%
  group_by(cluster) %>%
  slice_max(rank, n = 10) %>%
  #dplyr::filter(rank > 2.8) %>%
  #dplyr::select(-n) %>%
  ungroup()


DotPlot(lognormal_GE,
        features = c(tmp2$gene %>% unique, "Hbb-bt","Ttr"),
        scale = T,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))
```


```{r}
tmp_gaba <- c("Gad1", "Gad2", "Dlx6os1", "Slc32a1")
tmp_glut <-c("Slc17a7", "Slc17a6", "Nrn1")
tmp_astro <- c("Gja1", "S1pr1", "Fgfr3", "Gli2", "Ntsr2")
tmp_micro <- c("Cx3cr1", "Ikzf1", "C1qa", "C1qb", "C1qc")
tmp_olig <- c("Mal", "Mog", "Mag", "Prr5l", "Plp1")
tmp_OPC <- c("Pdgfra", "Cspg4", "Olig1", "Olig2", "C1ql1")
tmp_other <- c("Tmem72", "Otx2os1", "Mecom", "Otx2", "Folr1")
tmp_vascular <- c("Vtn", "Cldn5", "Pdgfrb", "Flt1", "Acta2", "Ly6c1")
tmp_blood <- c("Hbb-bs", "Hbb-bt", "Hba-a2")#, "Ttr")
```

#### Cluster Markers 2

```{r}
i <- "pca_snn_res.0.8"
Idents(lognormal_GE) = i
test <- AggregateExpression(lognormal_GE,
                            return.seurat = T,
                            assay = "SCT",
                            group.by= i)$SCT$count

tmp_neuronal <- c(tmp_gaba, tmp_glut)
tmp_non <- c(tmp_astro, tmp_micro, tmp_olig, tmp_OPC, tmp_other, tmp_vascular)

tmp_neuronal <- tmp_neuronal[tmp_neuronal %in% rownames(test)]
tmp_non <- tmp_non[tmp_non %in% rownames(test)]
tmp_vascular <- tmp_vascular[which(tmp_vascular %in% rownames(test))]
```


```{r}
test[test==0] <- NA
  
test2 <- data_frame(
  neuronal = Matrix::colMeans(test[tmp_neuronal, ], na.rm = T),
  non = Matrix::colMeans(test[tmp_non, ], na.rm = T),
  gaba = Matrix::colMeans(test[tmp_gaba, ], na.rm = T),
  glut = Matrix::colMeans(test[tmp_glut, ], na.rm = T),
  micro = Matrix::colMeans(test[tmp_micro, ], na.rm = T),
  astro = Matrix::colMeans(test[tmp_astro, ], na.rm = T),
  olig = Matrix::colMeans(test[tmp_olig, ], na.rm = T),
  OPC = Matrix::colMeans(test[tmp_OPC, ], na.rm = T),
  imm = Matrix::colMeans(test[tmp_other, ], na.rm = T),
  vascular = Matrix::colMeans(test[tmp_vascular, ], na.rm = T),
  seurat_cluster = 1:(ncol(test))) %>%
  mutate(across(where(is.numeric), function(x) replace_na(x, 0)))

#Hochberg called doublet clusters by looking for ratios greater than 2..

#Clusters that are both neuronal and non
test2 <- test2 %>%
  rowwise() %>%
  dplyr::mutate(
    cluster_class = colnames(test2)[3:10][which.max(c_across(gaba:vascular))],
    cluster_class_broad = if_else(cluster_class %in% c("glut", "gaba"), cluster_class, "non")
  ) %>%
  ungroup()
  

test2 %>%
  dplyr::count(cluster_class_broad)

test2 %>%
  dplyr::count(cluster_class)
#This seems about right..

test2 <- test2 %>%
  dplyr::mutate(
    worst_cluster_ratio = case_match(
      cluster_class,
    "gaba" ~ apply(cbind(gaba/glut, gaba/astro, gaba/micro,
                          gaba/olig, gaba/OPC, gaba/imm,
                          gaba/vascular),
                     1,
                     min, na.rm = T),
      "glut" ~ apply(cbind(glut/gaba, glut/astro, glut/micro,
                           glut/olig, glut/OPC, glut/imm,
                           glut/vascular),
                     1,
                     min, na.rm = T),
      "astro" ~ apply(cbind(astro/gaba, astro/glut, astro/micro,
                            astro/olig, astro/OPC, astro/imm,
                            astro/vascular),
                     1,
                     min, na.rm = T),
      "micro" ~ apply(cbind(micro/gaba, micro/glut, micro/astro,
                            micro/olig, micro/OPC, micro/imm,
                            micro/vascular),
                     1,
                     min, na.rm = T),
      "olig" ~ apply(cbind(olig/gaba, olig/glut, olig/micro,
                           olig/astro, olig/OPC, olig/imm,
                           olig/vascular),
                     1,
                     min, na.rm = T),
      "OPC" ~ apply(cbind(OPC/gaba, OPC/glut, OPC/micro,
                          OPC/olig, OPC/astro, OPC/imm,
                          OPC/vascular),
                     1,
                     min, na.rm = T),
      "imm" ~ apply(cbind(imm/gaba, imm/glut, imm/micro,
                          imm/olig, imm/astro, imm/OPC,
                          imm/vascular),
                     1,
                     min, na.rm = T),
     "vascular" ~ apply(cbind(vascular/gaba, vascular/glut, vascular/micro,
                         vascular/olig, vascular/astro, vascular/OPC,
                         vascular/imm),
                   1,
                   min, na.rm = T),
          .default = NA
    )) %>%
  dplyr::mutate(
    doublet = worst_cluster_ratio < 3
  ) %>% ungroup()

test2 %>% dplyr::count(cluster_class, doublet)

test2 <- left_join(lognormal_GE@meta.data %>%
                     dplyr::select(-starts_with(colnames(test2))),
                   test2 %>%
                     mutate(seurat_cluster = factor(seurat_cluster)),
                   by = c("pca_snn_res.0.8" = "seurat_cluster"))

rownames(test2) <- lognormal_GE[["SCT"]]@data@Dimnames[[2]]
test2$cells_old <- test2$cells
test2$cells <- rownames(test2)

lognormal_GE@meta.data <- test2

rm(test, test2)
```

```{r}
lognormal_GE@meta.data$cluster_class_dotplot <- lognormal_GE@meta.data$cluster_class_broad
lognormal_GE@meta.data$cluster_class_dotplot[lognormal_GE@meta.data$doublet] <- "doublet"
lognormal_GE@meta.data$cluster_class_dotplot[lognormal_GE@meta.data$cluster_class=="imm"] <- "other"

lognormal_GE@meta.data$cluster_class_dotplot <- factor(lognormal_GE@meta.data$cluster_class_dotplot,
                                                       levels = c("gaba", "glut","non",
                                                                  "doublet", "other"))

#table(lognormal_GE@meta.data$cluster_class_dotplot)
#Idents(object = lognormal_GE) <- "cell_class_dotplot"
#Idents(object = lognormal_GE) <- "pca_snn_res.0.5"

tmp_plot_2 <-
  DotPlot(lognormal_GE,
        features = c(tmp_gaba, tmp_glut, tmp_non),
        group.by = "cluster_class_dotplot",
        split.by = "cluster_class_dotplot",
        cols = c("grey50",
                 "#0072B4",
                 "#CC1C2F",
                 "#FDE333",
                 "#188B41")[c(5,2,3,1,4)],
        scale = T,
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  xlab("Marker Gene") +
  ylab("Cluster Cell Type") +
  scale_y_discrete(labels = c("Doublet\n(5 clusters, 1209 cells)",
                              "Glutamatergic\n(11 clusters, 12708 cells)",
                              "GABAergic\n(17 clusters, 25891 cells)",
                              "Non-Neuronal\n(5 clusters, 13931 cells)",
                              "Other \n(1 cluster, 100 cells)")[c(3,2,4,1,5)]
                   )
  
  tmp_plot_2
  
  lognormal_GE@meta.data %>%
  dplyr::select(cluster_class_dotplot, pca_snn_res.0.8) %>%
  distinct() %>% count(cluster_class_dotplot)

lognormal_GE@meta.data %>% count(cluster_class_dotplot) %>%
  mutate(p = n / sum(n))
  
Idents(object = lognormal_GE) <- "pca_snn_res.0.8"

tmp <- lognormal_GE@meta.data %>%
  mutate(test = unfactor(pca_snn_res.0.8)) %>%
  arrange(cluster_class_dotplot, as.numeric(test)) %>% 
  pull(test)

lognormal_GE@meta.data$clusters_dotplot <- factor(lognormal_GE@meta.data$pca_snn_res.0.8,
                                                  levels = unique(tmp))
  
DotPlot(lognormal_GE,
        features = c(tmp_gaba, tmp_glut, tmp_non),
        scale = F,
        cluster.idents = F,
        group.by = "clusters_dotplot",
        split.by = "cluster_class_dotplot",
        cols = c("grey50",
                 "#0072B4",
                 "#CC1C2F",
                 "#FDE333",
                 "#188B41")[c(5,2,3,1,4)]) +theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  xlab("Marker Gene") +
  ylab("Cluster Cell Type") 

DotPlot(lognormal_GE,
        features = c("Slc17a6", "Slc17a7", "Slc17a8",
                     "Slc32a1", "Slc18a2", "Gad1", "Gad2", "Aldh1a1",
                     "Slc6a5",
                     "Slc18a3", "Chat",
                     "Slc6a3", "Th", "Ddc",
                     "Slc6a4", "Tph2",
                     "Slc6a2", "Dbh",
                     "Hdc"),
        scale = F,
        cluster.idents = F,
        group.by = "clusters_dotplot",
        split.by = "cluster_class_dotplot",
        cols = c("grey50",
                 "#0072B4",
                 "#CC1C2F",
                 "#FDE333",
                 "#188B41")[c(5,2,3,1,4)]) +theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  xlab("Marker Gene") +
  ylab("Cluster Cell Type") 


```
I kind of want to dig into the doublets a bit more...




```{r}
DimPlot(
  lognormal_GE,
  cells = (lognormal_GE@meta.data %>% 
             filter(doublet) %>%
             rownames()),
  reduction = "umap",
  split.by = "clusters_dotplot",
  ncol = 8,
  label = T
)
```

Of the doublets, only 37 is actually spatially resolved enough to really judge.

##### UMAP
```{r}
test <- lognormal_GE@reductions$umap@cell.embeddings
test2 <- lognormal_GE@meta.data

all.equal(rownames(test),
          rownames(test2))

test <- bind_cols(test,test2)

test2 <- test %>%
  group_by(pca_snn_res.0.8) %>%
  summarise(n = n(),
            cluster_class_dotplot = unique(cluster_class_dotplot)) %>%
  group_by(cluster_class_dotplot) %>%
  arrange(desc(n)) %>%
  mutate(cluster_color = factor(paste0(str_to_title(cluster_class_dotplot),"_",formatC(row_number(),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  ungroup()

test <- left_join(test, test2 %>% dplyr::select(cluster_color, pca_snn_res.0.8))
  #mutate(cell_color = if_else(doublet == F, cell_color, "Doublet"),
  #       cell_class = if_else(doublet == F, cell_class, "Doublet"))

tmp_plot_1 <- 
  ggplot(test %>%
         group_by(pca_snn_res.0.8) %>%
         slice_sample(prop = 0.25)) +
  geom_point(aes(x = umap_1,
                 y = umap_2,
                 color = cluster_color,
                 shape = cluster_class_dotplot)) +
  ggrepel::geom_text_repel(data = test %>%
                     group_by(pca_snn_res.0.8) %>%
                     summarize(x = median(umap_1),
                               y = median(umap_2),
                               cluster_class_dotplot = unique(cluster_class_dotplot),
                               cluster_color = unique(cluster_color)) %>%
                    dplyr::filter(cluster_class_dotplot != "doublet") %>%
                    group_by(cluster_class_dotplot) %>%
                    slice_min(pca_snn_res.0.8, n = 3) %>%
                    ungroup(),
                   aes(x = x, y = y, label = cluster_color),
                  size = 5,
                  min.segment.length = 0.1,
                  max.overlaps = 15,
                  box.padding = 2,

                  color = "grey25"
                   ) +
  theme_JQA() +
  scale_color_manual(
    values = c(
      hcl.colors(22,palette="Blues3")[1:17],
      hcl.colors(15,palette="Reds3")[1:11],
      hcl.colors(10,palette="Greens3")[1:5],
      rep("grey50", 5),
      "#FDE333"
    ),
  ) +
  scale_shape_manual(
                     values = c(16,16,16,16,16),
                     labels = c("Doublet",
                                "GABAergic",
                                "Glutamatergic",
                                "Other",
                                "Non-Neuronal")
                     ) +
  guides(color = "none",
         shape = guide_legend("Cluster Cell Type",
                              override.aes = list(
                                size = 5,
                                shape = c(15,15,15,15,15),
                                color = c("grey50",
                                          "#0072B4",
                                          "#CC1C2F",
                                          "#FDE333",
                                          "#188B41")
                                )
                              )
  ) +
  xlab("Umap Dimension 1 (arb)") +
  ylab("Umap Dimension 2 (arb)")

tmp_plot_1
```

```{r}
wrap_plots(
  tmp_plot_1,
  tmp_plot_2,
  guides = "collect"
) + plot_annotation(tag_levels = "A")
```


#### Cell Markers

I there is some motivation to instead do this on the cell level... since the clusters arent.. super clean? I.E. some of those doublet clusters are just because of weird clustering (grabbing cells from multiple clusters), rather than being true doublets.  At the very least we can do doublet calling on the cell level.

Dayne and Nihal used scSorter.. we will follow https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html
```{r}
anno <- list(
  gaba = data_frame(Marker = tmp_gaba),
  glut = data_frame(Marker =tmp_glut),
  astro = data_frame(Marker =tmp_astro),
  micro = data_frame(Marker =tmp_micro),
  olig = data_frame(Marker =tmp_olig),
  OPC = data_frame(Marker =tmp_OPC),
  #blood = data_frame(Marker = tmp_blood),
  #vascular = data_frame(Marker = tmp_vascular),
  other = data_frame(Marker = tmp_other)
) %>% bind_rows(.id = "Type") %>%
  mutate(Weight = 2)

topgenes <- head(VariableFeatures(lognormal_GE), 2000)
expr = GetAssayData(lognormal_GE, assay = "SCT")

topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

#At last, we subset the preprocessed expression data. Now, we are ready to run scSorter.

picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]
test <- as.matrix(expr)
coinda 
rts <- scSorter::scSorter(test, as.data.frame(anno))

table(rts$Pred_Type)
```


```{r}
tmp_neuronal <- c(tmp_gaba, tmp_glut)
tmp_non <- c(tmp_astro, tmp_micro, tmp_olig, tmp_OPC, tmp_other, tmp_vascular)

tmp <- unique(c(tmp_neuronal,
                 tmp_non))
test <- GetAssayData(object = lognormal_GE, layer = "count", assay = "SCT")

tmp <- tmp[tmp %in% rownames(test)]
test <- test[tmp, ]

tmp_blood <- tmp_blood[which(tmp_blood %in% tmp)]
tmp_vascular <- tmp_vascular[which(tmp_vascular %in% tmp)]
tmp_non <- tmp_non[which(tmp_non %in% tmp)]

test[test==0] <- NA

test2 <- data_frame(
  neuronal = Matrix::colMeans(test[tmp_neuronal, ], na.rm = T),
  non = Matrix::colMeans(test[tmp_non, ], na.rm = T),
  gaba = Matrix::colMeans(test[tmp_gaba, ], na.rm = T),
  glut = Matrix::colMeans(test[tmp_glut, ], na.rm = T),
  micro = Matrix::colMeans(test[tmp_micro, ], na.rm = T),
  astro = Matrix::colMeans(test[tmp_astro, ], na.rm = T),
  olig = Matrix::colMeans(test[tmp_olig, ], na.rm = T),
  OPC = Matrix::colMeans(test[tmp_OPC, ], na.rm = T),
  imm = Matrix::colMeans(test[tmp_other, ], na.rm = T),
  vascular = Matrix::colMeans(test[tmp_vascular, ], na.rm = T),
  cells = colnames(test)) %>%
  mutate(across(where(is.numeric), function(x) replace_na(x, 0)))

test3 <- pivot_longer(test2,
                       cols = neuronal:vascular,
                       names_to = "cell_type",
                       values_to = "mean")

ggdensity(test3 %>%
            filter(mean > 0, mean < 30),
          facet.by = "cell_type",
          x = "mean")

#Clusters that are both neuronal and non
test2 <- test2 %>%
  dplyr::mutate(
    cell_class = colnames(test2)[3:10][apply(test2[,3:10], 1, which.max)],
    cell_class_broad = if_else(cell_class %in% c("glut", "gaba"), cell_class, "non")
  ) %>%
  dplyr::mutate(
    cell_class = if_else((non+neuronal) == 0, "other", cell_class),
    cell_class_broad = if_else((non+neuronal) == 0, "other", cell_class_broad),
  )

test2 %>%
  dplyr::count(cell_class_broad)

test2 %>%
  dplyr::count(cell_class)
#This seems about right..

test2 <- test2 %>%
  dplyr::mutate(
    worst_cell_ratio = case_match(
      cell_class,
     "gaba" ~ apply(cbind(gaba/glut, gaba/astro, gaba/micro,
                          gaba/olig, gaba/OPC, gaba/imm,
                          gaba/vascular),
                     1,
                     min, na.rm = T),
      "glut" ~ apply(cbind(glut/gaba, glut/astro, glut/micro,
                           glut/olig, glut/OPC, glut/imm,
                           glut/vascular),
                     1,
                     min, na.rm = T),
      "astro" ~ apply(cbind(astro/gaba, astro/glut, astro/micro,
                            astro/olig, astro/OPC, astro/imm,
                            astro/vascular),
                     1,
                     min, na.rm = T),
      "micro" ~ apply(cbind(micro/gaba, micro/glut, micro/astro,
                            micro/olig, micro/OPC, micro/imm,
                            micro/vascular),
                     1,
                     min, na.rm = T),
      "olig" ~ apply(cbind(olig/gaba, olig/glut, olig/micro,
                           olig/astro, olig/OPC, olig/imm,
                           olig/vascular),
                     1,
                     min, na.rm = T),
      "OPC" ~ apply(cbind(OPC/gaba, OPC/glut, OPC/micro,
                          OPC/olig, OPC/astro, OPC/imm,
                          OPC/vascular),
                     1,
                     min, na.rm = T),
      "imm" ~ apply(cbind(imm/gaba, imm/glut, imm/micro,
                          imm/olig, imm/astro, imm/OPC,
                          imm/vascular),
                     1,
                     min, na.rm = T),
     "vascular" ~ apply(cbind(vascular/gaba, vascular/glut, vascular/micro,
                         vascular/olig, vascular/astro, vascular/OPC,
                         vascular/imm),
                   1,
                   min, na.rm = T),
      .default = NA
    )) %>%
  dplyr::mutate(
    doublet = worst_cell_ratio < 2
  )
test2 %>% dplyr::count(doublet)
test2 %>% dplyr::count(cell_class_broad,doublet)
test2 %>% dplyr::count(cell_class,doublet)
#So, a few thousand potential doublets if we count at the cell level - which maybe is problematic.. its probably more like "lower quality cells"

test2 %>% dplyr::filter(doublet == T) %>% View
```


```{r}
test3 <- left_join(lognormal_GE@meta.data %>%
                     mutate(cluster_doublet = doublet) %>%
                     dplyr::select(-neuronal, -non, -gaba, -glut,-micro,-astro,-olig,
                                   -OPC, -imm, -vascular, -cells, -doublet) %>%
                     rownames_to_column(var = "cells"),
                   test2 %>%
                     mutate(cell_doublet = doublet) %>%
                     dplyr::select(-neuronal, -non, -gaba, -glut,
                                   -micro,-astro,-olig,-OPC, -imm,- vascular, -doublet)
                   )
colnames(test3)

rownames(test3) <- rownames(lognormal_GE@meta.data)

lognormal_GE@meta.data <- test3
rm(test3)
```


```{r}
lognormal_GE@meta.data %>%
  dplyr::filter(is.na(cell_doublet)) %>% View

lognormal_GE@meta.data %>% count(cell_class_broad, cell_doublet) %>%
  mutate(p = n / sum(n))

lognormal_GE@meta.data %>% count(cell_doublet) %>%
  mutate(p = n / sum(n))

lognormal_GE@meta.data  %>% dplyr::count(cell_class_broad,cell_doublet) %>%
  ungroup() %>%
  group_by(cell_class_broad) %>%
  dplyr::mutate(p = n / sum(n)*100)

```

```{r}
ggplot(lognormal_GE@meta.data) +
    geom_bar(aes(x=pca_snn_res.0.8, fill=(cell_doublet)), position=position_fill()) + theme_JQA()

ggplot(lognormal_GE@meta.data) +
    geom_bar(aes(x=pca_snn_res.0.8, fill=(cell_class)), position=position_fill()) + theme_JQA()

lognormal_GE@meta.data %>%
  dplyr::select(pca_snn_res.0.8,
                cluster_class, cluster_doublet) %>% distinct %>% View

ggplot(lognormal_GE@meta.data) +
    geom_bar(aes(x=id, fill=(cell_class)), position=position_fill()) + theme_JQA()

ggplot(lognormal_GE@meta.data) +
    geom_bar(aes(fill=id, x=(cell_class)), position=position_fill()) + theme_JQA()

#immature neuronal group is 50% hs_12, but it is in all samples

lognormal_GE@meta.data %>%
  group_by(pca_snn_res.0.5) %>%
  dplyr::count(cell_class) %>%
  mutate(prop = n / sum(n),
         cells = sum(n)) %>%
  slice_max(prop, n=1) %>% 
  arrange(prop) 

lognormal_GE@meta.data %>%
  group_by(pca_snn_res.0.5) %>%
  dplyr::count(cell_doublet) %>%
  mutate(prop = n / sum(n),
        cells = sum(n)) %>%
  dplyr::filter(cell_doublet == T) %>%
  arrange(desc(prop))

```

#### Doublets

```{r}
bp <-  BiocParallel::MulticoreParam(6, RNGseed=saved_seed)

sce <- scDblFinder(
  GetAssayData(object = lognormal_GE, layer = "counts", assay = "SCT"),
  clusters = lognormal_GE@meta.data$pca_snn_res.0.8,
  #clusters  = F,
  samples = lognormal_GE@meta.data$id,
  BPPARAM = bp
  )

tab <- findDoubletClusters(GetAssayData(object = lognormal_GE, layer = "counts", assay = "RNA"),
                           clusters = lognormal_GE@meta.data$pca_snn_res.0.8)

View(as.data.frame(tab))
table(sce$scDblFinder.class)

sce$scDblFinder.mostLikelyOrigin %>% table %>%

table(lognormal_GE@meta.data$cluster_doublet)
```



```{r}
tmp_file = "lognormal_seurat_v3.RData"

if(!file.exists(tmp_file)) {
  save(lognormal_GE, overall_markers, file = tmp_file)
} else {
  load(tmp_file)
}

rm(tmp_file)
```


```{r}
#write.xlsx(overall_markers, file = "OVERALLmarkers_jqa_111224.xlsx")

```


# GABA

We start over from the raw data
```{r}
GABA_GE <- subset(filtered_GE,
                  subset = (id != "HS_3") &
                    (id != "HS_16"))

#Fix the cell names
test4 <- lognormal_GE@meta.data %>%
  dplyr::select(cells, cell_class, cell_class_broad, cell_doublet, cluster_class, cluster_class_broad, cluster_doublet) 
test4$cells <- str_split_i(test4$cells, "[0-9]_", 2)
test4 <- left_join(GABA_GE@meta.data, test4)
rownames(test4) <- test4$cells
GABA_GE@meta.data <- test4

#We keep any cell that was in a GABAergic cluster, or was identified as GABAergic.. but that wasn't also identified as non-neuronal or glutamatergic
GABA_GE <- subset(GABA_GE,
                   subset = 
                       ((cell_class == "gaba") |
                       (cluster_class == "gaba"))  &
                    #(cell_doublet == F) &
                    (cell_class_broad != "non") &
                    (cell_class != "glut"))

GABA_GE

rm(test,test2,test3,test4, tmp, tmp_markers, tmp_meta, tmp_plot_1, tmp_plot_2)
```

We are left with 26181 cells.  We will filter by low expression again, though.

### Cell Filtering
```{r}
counts <- GetAssayData(object = GABA_GE, layer = "counts")

nonzeros <- Matrix::rowSums(counts>2)

nonzeros <- data_frame(gene = rownames(counts),
                       ncells = nonzeros)

nonzeros <- nonzeros %>%
arrange(ncells) %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) 

sum(nonzeros$ncells > 10)

nonzeros %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) %>%
  ggplot() +
  geom_point(aes(y = idx, x = ncells)) +
  theme_JQA() +
    coord_cartesian(xlim=c(0,10000))

nonzeros[c(
  which.max(nonzeros$idx>10),
  which.max(nonzeros$idx>20),
  which.max(nonzeros$idx>30),
  which.max(nonzeros$idx>40),
  which.max(nonzeros$idx>50),
  which.max(nonzeros$idx>60),
  which.max(nonzeros$idx>70),
  which.max(nonzeros$idx>80),
  which.max(nonzeros$idx>90)), ]


#This leaves us with ~ 14k genes which is great  
to_drop <- nonzeros$ncells < 10

sum(!to_drop)
#that leaves 13k genes, which is fine I think

GABA_GE <- subset(x = GABA_GE,
                      features = nonzeros$gene[!to_drop])
GABA_GE 
rm(counts,nonzeros,to_drop)
```

26181 GABA Neurons, 12322 genes

### Normalization
```{r}
GABA_GE <- SplitObject(GABA_GE, split.by = "id")
```

```{r}
plan(sequential)
gc()
plan(multisession, workers = 7)
options(future.globals.maxSize = 80 * 1024 ^ 3)
GABA_GE <- future.apply::future_lapply(GABA_GE,
                                      function(x) {
                                        SCTransform(x,
                                                    assay = "RNA",
                                                    vars.to.regress = "percent.mt",
                                                    verbose = FALSE,
                                                    seed.use = saved_seed,
                                                    )
                                      },
                                      future.seed = saved_seed
) 

GABA_GE <- future.apply::future_lapply(GABA_GE,
                                            function(x) {
                                              NormalizeData(x,
                                                            assay = "RNA",
                                              )
                                            },
                                            future.seed = saved_seed
) 

GABA_GE <- merge(GABA_GE[[1]],
              GABA_GE[2:14],
              add.cell.ids = names(GABA_GE))

GABA_GE[["RNA"]] <- JoinLayers(GABA_GE[["RNA"]])

VariableFeatures(GABA_GE, assay = "SCT") <- rownames(GABA_GE[["SCT"]]@scale.data)

# lognormal_GE <- NormalizeData(subset(filtered_GE,
#                                      subset = (id != "HS_3") &
#                                        (id != "HS_16")))
# 
# lognormal_GE <- FindVariableFeatures(lognormal_GE)
# 
# lognormal_GE <- ScaleData(lognormal_GE)
plan(sequential)
gc()
```

```{r}
GABA_GE <- RunPCA(GABA_GE,
                  assay = "SCT")
```


```{r}
ElbowPlot(GABA_GE,ndims = 40)
```

```{r}
DimHeatmap(GABA_GE,dims=1:10, cells=500)
DimHeatmap(GABA_GE,dims=11:20, cells=500)
DimHeatmap(GABA_GE,dims=21:30, cells=500)
DimHeatmap(GABA_GE,dims=31:40, cells=500)
```

```{r}
GABA_GE <- RunUMAP(GABA_GE, 
                        dims = 1:40,
                        reduction = "pca",
                        reduction.name = "umap",
                        seed.use = saved_seed)
```

### Clustering
```{r}
GABA_GE <- FindNeighbors(GABA_GE,
                              dims = 1:40,
                              reduction = "pca",
                              graph.name = paste0("pca_", c("nn","snn")))
```

```{r}
library(reticulate)
numpy <- reticulate::import("numpy")
pandas <- reticulate::import("pandas")
leidenalg <- reticulate::import("leidenalg")

plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)

GABA_GE <- FindClusters(GABA_GE,
                             method  = "igraph",
                             algorithm = 4,
                             random.seed = saved_seed,
                             graph.name = "pca_snn",
                             resolution = c(0.01,0.05,0.1,0.25,0.5,0.8, 1, 1.5))

plan(sequential)
gc()
```

```{r}
clustree(GABA_GE, prefix = "pca_snn_res.")
```

This is a very pretty one, clean clustering. I like 5 the most, as it breaks up the big cluster. 0.5 is the first break of th second largets cluster.

```{r}
DimPlot(GABA_GE,
        group.by = c("genotype",
                     "id",
                     "class",
                     "GEX_lane",
                     "pca_snn_res.0.01",
                     "pca_snn_res.0.05",
                     "pca_snn_res.0.1",
                     "pca_snn_res.0.25",
                     "pca_snn_res.0.5"),
        alpha = 0.5,
        label = T)
```

```{r}
ggplot(GABA_GE@meta.data) +
    geom_bar(aes(x=pca_snn_res.0.5, fill=(cell_class)), position=position_fill()) + theme_JQA()

```

```{r}
GABA_GE@meta.data %>% dplyr::count(pca_snn_res.0.5)
```
Only one cluster below 100 cells



```{r}
DimPlot(GABA_GE,
        group.by = c("pca_snn_res.0.5"),
        label = T)
```

I would trust 0.01, 0.05, 0.1 and 0.25 fine. I'm actually going to go with 5, because it isolates the cluster in the upper left, we can always combine clusters back together.

```{r}
GABA_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.5, fill=(class)), position=position_fill())+ theme_JQA()
#really only two mixed clusters on cell_class

GABA_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.5, fill=(id)), position=position_fill())+ theme_JQA()
#still some clusters completely dominated by a single sample

GABA_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.5, fill=(genotype)), position=position_fill())+ theme_JQA()

#We have two fewer HP mice now.

GABA_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.5, fill=(sex)), position=position_fill())+ theme_JQA()
```

```{r}
Idents(object = GABA_GE) <- "pca_snn_res.0.5"
```

```{r}
DimPlot(GABA_GE,
        group.by = "pca_snn_res.0.5",
        label = T)
```

```{r}
test <- GABA_GE@meta.data

test2 <- sapply(unique(test$id), function(sample) {
  sapply(as.numeric(levels(test$pca_snn_res.0.5)), function(cluster) {
  fisher.test(test$id == sample,
              test$pca_snn_res.0.5 == cluster,
              alternative = "g")$p.value
  })
})

test3 <- -log10(test2)
test3[test3 == Inf] <- 100
test3[test3 == NA] <- 0
test3[test3 > 25] <- 25
rownames(test3) <- as.numeric(levels(test$pca_snn_res.0.5))

test4 <- data_frame(id = colnames(test3)) %>%
  left_join(sampleInfo)

Heatmap(test3,
        cluster_rows = T,
        row_title_side = "right",
        top_annotation = HeatmapAnnotation(Sex = test4$sex,
                                           Genotype = test4$genotype),
        col = colorRamp2(c(0, 1.3,
                     3.475, 5.650, 7.825, 10, 25), c("black", "white",
                                                 "#00B785","#53CC67","#B2DC3C","#FDE333", "red")))

```
18,7,14,1,21,17 are all driven by HS_13, 8,5 are all driven by a single sample

```{r}
FeaturePlot(
  GABA_GE,
  reduction = "umap",
  features = c("nFeature_RNA",
               "nCount_RNA",
               "nCount_SCT",
               "nFeature_SCT",
               "percent.hge",
               "percent.mt",
               "percent.rp"),
  order = TRUE,
  pt.size = 0.4,
  min.cutoff = 'q10',
  label = T
)
```
None of the clusters really stand out.  Maybe cluster 1 is a little higher in nCount/nFeatures

```{r}
GABA_GE@meta.data %>%
  dplyr::count(pca_snn_res.0.5)

VlnPlot(GABA_GE, features = "nCount_RNA", group.by = "pca_snn_res.0.5", assay = "RNA") /
VlnPlot(GABA_GE, features = "nFeature_RNA", group.by = "pca_snn_res.0.5", assay = "RNA")/
VlnPlot(GABA_GE, features = "nCount_SCT", group.by = "pca_snn_res.0.5", assay = "RNA") /
VlnPlot(GABA_GE, features = "nFeature_SCT", group.by = "pca_snn_res.0.5", assay = "RNA")  
#15 and 21 seem to have double the feature count..

VlnPlot(GABA_GE,
        assay = "SCT",
        layer = "data",
        features = c("Gad1", "Gad2", "Slc32a1", "Slc17a7", "Slc17a6"), group.by = "pca_snn_res.0.5", ncol = 2)
```

```{r}
GABA_GE_orig <- GABA_GE
```

```{r}
plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)
GABA_GE <- PrepSCTFindMarkers(GABA_GE)
plan(sequential)
gc()

#5746
```

## Cluster DE Markers
```{r}
Idents(object = GABA_GE) <- "pca_snn_res.0.5"

plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)

tmp_markers_gaba <- future.apply::future_lapply(1:n_distinct(GABA_GE@meta.data$pca_snn_res.0.5),
                             function(i)
                             {
                               SeuratWrappers::RunPresto(GABA_GE,
                                                         assay = "SCT",
                                                         ident.1 = i)
                             },
                             future.seed = saved_seed
                             )

plan(sequential)
gc()


tmp_markers_gaba <- lapply(tmp_markers_gaba, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})

#write.xlsx(tmp_markers_gaba, file = "GABAmarkers_jqa_103124.xlsx")
```

#### Dendogram
```{r}
tmp_tree <- BuildClusterTree(GABA_GE,
                          assay = "SCT",
                          slot  = "scale.data",
                          dim = 1:40
                          )

PlotClusterTree(tmp_tree, 
                direction = "rightwards")

rm(tmp_tree)
```

```{r}
GABA_GE@meta.data$pca_snn_res.0.5 %>% table
```
23 an 24 are below 10/sample

#### Dotplot doublets
Let's check for doublets against our known types, since we didn't drop any earlier
```{r}
(DotPlot(GABA_GE,
        features = c(tmp_gaba, tmp_glut, tmp_non),
        group.by = c("pca_snn_res.0.5"),
        scale = F,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))) %>%
  plotly::ggplotly()

#18 is gaba+glut
#19 is gaba+opc
#17 is a little noisy, maybe olig

DimPlot(GABA_GE,
        split.by = "pca_snn_res.0.5",
        reduction = "umap",
        ncol = 5)
```
18+19 are looking like fairly robust doublets.. 

#### Dotplot top markers
```{r}
tmp <- tmp_markers_gaba %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

#the number of markers we want per cluster
tmp2 <- tmp %>%
  #mutate(cluster = if_else(cluster %in% c("6", "9", "10"), "6", cluster)) %>%
  group_by(cluster) %>%
  slice_max(rank, n = 10) %>%
  #dplyr::filter(rank > 2.8) %>%
  #dplyr::select(-n) %>%
  ungroup()
  
tmp2 %>% count(cluster) %>% arrange(n)

(DotPlot(GABA_GE,
        features = unique(tmp2$gene),
        group.by = c("pca_snn_res.0.5"),
        scale = T,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))) %>%
  plotly::ggplotly()
```

Drop the obvious doublets and we are going to keep 23 for now, as it looks like it is maybe cholinergic...
```{r}
GABA_GE <- subset(GABA_GE, 
                  subset = (pca_snn_res.0.5 != "19") &
                    (pca_snn_res.0.5 != "18")
                  )
```

We can also merge (7/12, 10/9) as their expression profiles are consistent. I think I will also throw 24 back into 5 since it is only 73 cells.. we could also throw it back into 14 which is its closest pair in the dendrogram... 5 is how it was split when we were doing voerall clustering, though, so lets honor that.

```{r}
tmp <- GABA_GE@meta.data

tmp <- tmp %>%
    mutate(pca_snn_res.0.5 = if_else(pca_snn_res.0.5 %in% c("7", "12"), "7", pca_snn_res.0.5),
           pca_snn_res.0.5 = if_else(pca_snn_res.0.5 %in% c("10", "9"), "9", pca_snn_res.0.5),
           pca_snn_res.0.5 = if_else(pca_snn_res.0.5 %in% c("24"), "5", pca_snn_res.0.5)
           )

tmp2 <- tmp %>%
  count(pca_snn_res.0.5) %>%
  arrange(desc(n)) %>%
  mutate(new = factor(row_number()))

tmp2

tmp <- left_join(tmp, tmp2)
tmp$pca_snn_res.0.5 <- tmp$new
tmp <- tmp %>%
  dplyr::select(-new, -n)
rownames(tmp) <- GABA_GE[["SCT"]]@data@Dimnames[[2]]

GABA_GE@meta.data <- tmp

#Do we need to do this?
Idents(object = GABA_GE) <- "pca_snn_res.0.5"
```

Let's look at that cholinergic population
```{r}
tmp %>%
  filter(pca_snn_res.0.5 == "19") %>%
  count(cell_class)

tmp %>%
  filter(pca_snn_res.0.5 == "19") %>%
  count(cluster_class)
```

They were gaba.

Now lets go re-run the cluster thing, we have to set recorrect_umi to false here, since we subset a SCT stack that we previoulsly ran prepsctmarkers on

## Cluster DE Markers 2
```{r}
plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)

tmp_markers_gaba <- future.apply::future_lapply(1:n_distinct(GABA_GE@meta.data$pca_snn_res.0.5),
                             function(i)
                             {
                               SeuratWrappers::RunPresto(GABA_GE,
                                                         assay = "SCT",
                                                         ident.1 = i,
                                                         recorrect_umi = F)
                             },
                             future.seed = saved_seed
                             )

plan(sequential)
gc()

tmp_markers_gaba <- lapply(tmp_markers_gaba, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})

#write.xlsx(tmp_markers_gaba, file = "GABAmarkers_jqa_103124.xlsx")
```

#### Dendogram 2
```{r}
GABA_GE <- BuildClusterTree(GABA_GE,
                          assay = "SCT",
                          slot  = "data",
                          reduction = "pca",
                          dim = 1:40
                          )

PlotClusterTree(GABA_GE) 

plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)
tree_markers <- future.apply::future_lapply(
  GABA_GE@tools$BuildClusterTree[["edge"]][,1] %>% unique(),
  function(i) {
    SeuratWrappers::RunPresto(GABA_GE,
                              assay = "SCT",
                              ident.1 = GABA_GE@tools$BuildClusterTree,
                              ident.2 = i,
                              recorrect_umi = F)
  },
  future.seed = saved_seed
)
plan(sequential)
gc()

names(tree_markers) <- GABA_GE@tools$BuildClusterTree[["edge"]][,1] %>% unique()

tree_markers <- lapply(tree_markers, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})


#  library(ggtree)
# ggtree(tmp_tree@tools$BuildClusterTree,
#           branch.length="none") + 
#      geom_tiplab(as_ylab = F) + 
#      geom_nodelab()+
#      geom_rootedge(rootedge=1) +
#   geom_text(data = data_frame(x = 1:10, y = 1:10, label = 1:10), aes(x=x,y=y,label=label))
# 
# 
#      #geom_cladelab(node = 19, label = "test", align = T, hjust = -10   theme_tree()

```

```{r}
write.xlsx(tmp_markers_gaba, file = "GABAmarkers_jqa_121924.xlsx")
write.xlsx(tree_markers, file = "GABAtreemarkers_jqa_121924.xlsx")
```

#### Dotplot markers 2


```{r}
tmp_levels <- GABA_GE@tools[["BuildClusterTree"]][["edge"]][GABA_GE@tools[["BuildClusterTree"]][["edge"]][, 2]<20, 2]

tmp_levels <- GABA_GE@tools[["BuildClusterTree"]]$tip.label[tmp_levels]

GABA_GE@meta.data$pca_snn_res.0.5 <- factor(GABA_GE@meta.data$pca_snn_res.0.5,
                                            levels = as.character(tmp_levels))

Idents(GABA_GE) <- "pca_snn_res.0.5"
```

##### Top markers, cluster
```{r}
tmp <- tmp_markers_gaba %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

tmp3 <- tmp %>%
  mutate(cluster = factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>% 
  group_by(cluster) %>%
  slice_max(rank, n = 2) 

tmp_plot_1 <- wrap_plots(
  
  ggtree::ggtree(GABA_GE@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree(),
  
  DotPlot(GABA_GE,
          features = unique(c(tmp3$gene)),
          group.by = c("pca_snn_res.0.5"),
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)

tmp_plot_1

```

##### Top markers, tree
```{r}
tmp <- tree_markers %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

#tmp2 <- tmp2[!(tmp2$gene %in% unique(tmp2$gene[duplicated(tmp2$gene)])), ]

tmp3 <- bind_rows(
  tmp %>%
  group_by(cluster) %>%
  slice_max(rank, n = 10),
    tmp %>%
  group_by(cluster) %>%
  slice_min(rank, n = 10)
) %>% arrange(cluster,rank)

tmp_plot_1 <- wrap_plots(
  
  ggtree::ggtree(GABA_GE@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree(),
  
  DotPlot(GABA_GE,
          features = (unique(tmp3 %>% filter(cluster == "25") %>% pull(gene))),
          group.by = c("pca_snn_res.0.5"),
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)

tmp_plot_1

```
##### Chosen Markers

```{r}
tmp <- tmp_markers_gaba %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

tmp_features <- c(
  "Rai14", "Chrm2",
  "Foxp2", "Prkcd",
  "Ntsr1", "Nts",
  "Tacr1", "Tac1","Isl1",
  "Calcrl", "Adora2a",
  "Drd2", "Drd1",
  "Cd44", 
  "Nos1","Cemip", "Hapln1",
  "Vipr2","Vip",
  "Npy1r","Npy",
  "Crh", "Crhr1",
  "Sst",  

  "Scn4b","Pdyn","Cck","Pnoc","Penk","Nr2f2",
  "Chat", "Slc5a7"
)

#tmp_features <- chosen_markers
#tmp_features <- tmp_features[tmp_features %in% rownames(GABA_GE)]

tmp_features <- tmp %>%
  dplyr::filter(gene %in% tmp_features) %>%
  group_by(gene) %>%
  slice_max(rank, n = 1) %>%
  dplyr::mutate(cluster =
                  factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>%
  arrange(cluster) %>%
  pull(gene)

tmp_plot_1 <- wrap_plots(
  
  ggtree::ggtree(GABA_GE@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree(),
  
  DotPlot(GABA_GE,
          features = tmp_features,
          group.by = c("pca_snn_res.0.5"),
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)

tmp_plot_1

```

##### BICCN Neuropiptides
https://www.nature.com/articles/s41586-023-06812-z#Sec10

```{r}
tmp <- tmp_markers_gaba %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

tmp_features <- c(
  "Adcyap1",	"Npff",
"Adm",	"Nps",
"Agrp",	"Npvf",
"Agt",	"Npw",
"Apln",	"Npy",
"Avp",	"Nts",
"Calca",	"Oxt",
"Calcb",	"Pdyn",
"Cartpt",	"Penk",
"Cck",	"Pmch",
"Cort",	"Pnoc",
"Crh",	"Pomc",
"Edn1",	"Prlh",
"Edn3",	"Pthlh",
"Gal",	"Pyy",
"Gcg",	"Rln3",
"Ghrh",	"Sst",
"Gnrh1",	"Tac1",
"Grp",	"Tac2",
"Hcrt",	"Trh",
"Kiss1",	"Ucn",
"Nmb",	"Ucn3",
"Nms",	"Uts2b",
"Nmu",	"Vip",
"Npb"
)

#tmp_features <- chosen_markers
#tmp_features <- tmp_features[tmp_features %in% rownames(GABA_GE)]

tmp_features <- tmp %>%
  dplyr::filter(gene %in% tmp_features) %>%
  group_by(gene) %>%
  slice_max(rank, n = 1) %>%
  dplyr::mutate(cluster =
                  factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>%
  arrange(cluster) %>%
  pull(gene)


DotPlot(GABA_GE,
        features = unique(tmp_features),
        group.by = c("pca_snn_res.0.5"),
        scale = T,
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))

```

##### Hochberg

```{r}
tmp_features <- c(
  "Tshz1",  "Foxp2",  "Ppp1r1b",
  "Prkcd",  "Crh",  "Adora2a",
  "Penk",  "Scn4b",  "Drd1",
  "Pdyn",  "Nts",  "Zfhx3",
  "Vdr",  "Avp",  "Gal",
  "Fign",  "Ucn3",  "Nxph1",
  "Lhx8",  "Th",  "Prlr",
  "Calcr",  "Cartpt",  "Moxd1",
  "Sst",  "Pvalb",  "Lamp5",
  "Sncg"
)

tmp <- tmp_markers_gaba %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(gene) %>%
  slice_max(rank, n = 100) %>%
  ungroup()


tmp_features <- tmp_features[tmp_features %in% unique(tmp$gene)]
tmp_features <- factor(tmp_features, levels = tmp_features)


clusters <- tmp %>%
  dplyr::filter(gene %in% tmp_features) %>%
  group_by(gene) %>%  
  slice_max(rank, n = 1) %>%
  ungroup() %>% group_by(cluster) %>%
  slice_max(rank, n = 1)

clusters2 <- tmp %>%
  dplyr::filter(gene %in% tmp_features,
                !(cluster %in% clusters$cluster)) %>%
  group_by(cluster) %>%
  slice_max(rank, n = 1)

clusters <- bind_rows(clusters, clusters2)

tmp_clusters <- data_frame(
  gene = tmp_features
) %>% left_join(
  clusters %>%
    dplyr::select(gene,cluster)
) %>% drop_na()

GABA_GE@meta.data$cluster <- factor(GABA_GE@meta.data$pca_snn_res.0.5,
                                         levels = tmp_clusters$cluster)

DotPlot(GABA_GE,
        features = rev(tmp_features),
        scale = T,
        group.by = "cluster",
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  xlab("Marker Gene") +
  ylab("Cluster Cell Type") + coord_flip()

```

###### Violin
```{r}
#Simple one...

tmp_features <- c("Chst9", "Vwc2l",
                  "Dlk1", "Prkcd",
                  "Dnah5", "Drd1",
                  "Tac1", "Drd2",
                  "Nek7", "Moxd1",
                  "Hapln1", "Adra1b",
                  "Vip", "Satb1",
                  "Lhx8", "Glra3",
                  "Lrmda", "Sema3c",
                  "Ebf1", "Rmst"
                  )


test <- GABA_GE[["SCT"]]@data[tmp_features, ] %>%
  t() %>% as.data.frame %>%
  rownames_to_column(var = "cells") %>%
  left_join(GABA_GE@meta.data %>% 
              dplyr::select(-cells) %>%
              rownames_to_column(var = "cells"))
#test <- split(test, test$pca_snn_res.0.5)

test <- pivot_longer(test, cols = all_of(tmp_features), names_to = "feature", values_to = "count") %>%
  mutate(feature = factor(feature, levels = tmp_features)) %>%
  dplyr::mutate(pca_snn_res.0.5 =
                  factor(pca_snn_res.0.5,
                         levels = c(6,18,16,4,7,14,9,12,15,8,17,2,3,11,5,13,1,10)))

test %>%
ggplot() +
  geom_violin(aes(y = count,
                  x = pca_snn_res.0.5,
                  group = pca_snn_res.0.5,
                  color = pca_snn_res.0.5),
              fill = "white",
              scale = "width",
              linewidth = 1,
              show.legend = F
              ) +
  ggforce::geom_sina(aes(y = count,
                  x = pca_snn_res.0.5,
                  group = pca_snn_res.0.5,
                  color = pca_snn_res.0.5),
              fill = "white",
              scale = "width",
              bins = 100,
              position = position_jitter(width = 0, height = 0.25),
              jitter_y = F,
              size = 0.5,
              shape = 16,
              show.legend = F) +
  geom_text(data = data_frame(
    y = apply(GABA_GE[["SCT"]]@data[tmp_features, ], 1, max)/2,
    max = apply(GABA_GE[["SCT"]]@counts[tmp_features, ], 1, max),
    feature = factor(names(max), levels = tmp_features),
    x = 0.25),
    aes(x = x, y = y, label = ceiling(max/10)*10)) +
  coord_flip(expand = F, xlim = c(0,19)) +
  facet_wrap(~feature, nrow = 1,
             strip.position = "top",
             scales = "free_x",
             dir = "v") +
  xlab("GABA Subtype") +
  ylab("Max Count") +
  theme_JQA() +
  scale_color_manual(values = hcl.colors(n = 25, palette = "Viridis")[1:18]) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        strip.placement = "output")

```

###### Dotplot
```{r}
DotPlot(GABA_GE,
        features = factor(tmp_features, levels = tmp_features),
        scale = T,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  xlab("Marker Gene") +
  ylab("") +
  guides(color = guide_colorbar(title = "Scaled\nExpression"),
         size = guide_legend(title = "Percent\nExpressed"))


GABA_GE@meta.data %>%
  dplyr::filter(pca_snn_res.0.5 == 22) %>% count(sample)

test %>%
   dplyr::mutate(pca_snn_res.0.5 =
                   factor(pca_snn_res.0.5,
                          levels = c(6,18,16,4,7,14,9,12,15,8,17,2,3,11,5,13,1,10))) %>%
   group_by(pca_snn_res.0.5, feature) %>%
   summarize(perc = sum(count>0)/n()*100,
             avg = mean(count)) %>%
   ungroup() %>%
 ggplot() +
   geom_point(aes(x = feature,
                  y = pca_snn_res.0.5,
                  size = perc,
                  color = avg),
              shape = 16) +
   scale_size_continuous(breaks = c(0,50,100),
                         limits = c(0,100),
                         range = c(0,6)) +
   ylab("") +
   xlab("") +
   theme_JQA() +
   scale_color_continuous_sequential(palette = "Blues 3")
```

```{r}

tmp_plot_1 <- wrap_plots(
    ggtree::ggtree(test2@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
  ggtree::theme_tree(),
  DotPlot(GABA_GE,
        features = factor(tmp_features, levels = tmp_features),
        scale = T,
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  xlab("Marker Gene") +
  ylab("")+
  guides(color = guide_colorbar(title = "Scaled\nExpression"),
         size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)
```

```{r}

tmp_plot_1 <- wrap_plots(
    ggtree::ggtree(tmp_tree@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
  ggtree::theme_tree(),
  DotPlot(GABA_GE,
        features = factor(tmp_features, levels = tmp_features),
        scale = T,
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  xlab("Marker Gene") +
  ylab("")+
  guides(color = guide_colorbar(title = "Scaled\nExpression"),
         size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)
```

##Comparing to HOchgerner
```{r}
hoch <- list()

tmp <- tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Ptk2b", "Wfs1", "Six3", "Meis2",
                                             "Foxp2", "Pax6", "Ntng1", "Tshz1",
                                             "Fmod", "Adra2a", "Col6a1", "Htr1f",
                                             "Enpp2", "Th", "Cyp26a1"))) %>%
  bind_rows(.id = "cluster") 

hoch[["A"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))

view(tmp)
view(hoch[["A"]])

tmp <-tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Ptk2b", "Wfs1", "Six3", "Meis2",
                                             "Ppr1rl1b",
                                             "Prkcd", "Adora2a", "Drd1", "Pdyn",
                                             "Oprk1", "Ezr", "Id4", "Scn4b", "Crh",
                                             "Ebf1", "Nts"))) %>%
  bind_rows(.id = "cluster") 

hoch[["B"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))

view(tmp)
view(hoch[["B"]])

tmp <- tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Ptk2b", "Wfs1", "Six3", "Meis2",
                                             "Tshz2", "Zfhx3",
                                             "Isl1", "Gpr101", "Vdr", "Gal",
                                             "Fign", "Ucn3", "Lrpprc", "Avp", "Nts",
                                             "Gabre", "Tac1", "Aldoc"))) %>%
  bind_rows(.id = "cluster")


hoch[["C"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))
view(tmp)
view(hoch[["C"]])

tmp <- tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Lhx6", "Tshz2", "Stab1", "Prl1",
                                             "Nxph1", "Sox6", "Moxd1", "Sst", "Pvalb",
                                             "Lhx8", "Th", "Cbln4", "Luzp2", "Nxph2",
                                             "Prlr", "Greb1", "Calcr", "Tac1", "St18",
                                             "Satb1", "Chodl", "Fign", "Npy", "Tmtc4",
                                             "Nek7", "Rpb4", "Vwc2", "Crabp1",
                                             "Etv1", "Pthlh"))) %>%
  bind_rows(.id = "cluster") 


hoch[["D"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))
view(tmp)
view(hoch[["D"]])


tmp_markers_gaba[[1]] %>%
  dplyr::filter(
    gene %in% c("Foxp2", "Fmod",
                "Adra2a",
                "Col6a1",
                "Htr1f",
                "Pax6") 
  )

test <- tmp_markers_gaba %>%lapply(function(x) {
  x %>% dplyr::filter(pct.1 > .75)}
  )

```

##### UMAP
```{r}
test <- GABA_GE@reductions$umap@cell.embeddings
test2 <- GABA_GE@meta.data

all.equal(rownames(test),
          rownames(test2))

test <- bind_cols(test,test2)

test2 <- test %>%
  group_by(pca_snn_res.0.5) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(cell_color = factor(paste0("Inh_",formatC(row_number(),
                                                           format = "d", flag = "0", digits = 1))))

test <- left_join(test, test2 %>% dplyr::select(cell_color, pca_snn_res.0.5))



tmp_plot_2 <- 
  ggplot(test %>%
         group_by(pca_snn_res.0.5) %>%
         slice_sample(prop = 0.25)) +
   geom_point(aes(x = umap_1,
                  y = umap_2,
                  color = cell_color),
              shape = 16) +
  ggrepel::geom_text_repel(data = test %>%
                     group_by(pca_snn_res.0.5) %>%
                     summarize(x = median(umap_1),
                               y = median(umap_2),
                               cell_color = unique(cell_color)),
                   aes(x = x, y = y, label = cell_color),
                  size = 5,
                  min.segment.length = 0.1,
                  max.overlaps = 15,
                  box.padding = 2,
                  color = "grey25") +
  theme_JQA() +
  scale_color_manual(
    values = c(
      hcl.colors(19,palette="Dark3")
    ),
  )  +
  xlab("Umap Dimension 1 (arb)") +
  ylab("Umap Dimension 2 (arb)") +
  theme(legend.position = "none") + ggtitle("Inhibitory Neuron Clusters")

tmp_plot_2

```

```{r}
(tmp_plot_2 + tmp_plot_1) + plot_annotation(tag_levels = "A")
```

```{r}
tmp_file = "GABA_seurat_v3.RData"

if(!file.exists(tmp_file)) {
  save(GABA_GE, tree_markers, tmp_markers_gaba, file = tmp_file)
} else {
  load(tmp_file)
}

rm(tmp_file)
```

## Supplementary tables

```{r}
GABA_GE@meta.data %>% 
  count(pca_snn_res.0.5,sample) %>%
  pivot_wider(names_from = sample,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Inh_",formatC(as.numeric(unfactor(pca_snn_res.0.5)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.5) %>% view

GABA_GE@meta.data %>% 
  count(pca_snn_res.0.5,genotype) %>%
  pivot_wider(names_from = genotype,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Inh_",formatC(as.numeric(unfactor(pca_snn_res.0.5)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.5) %>% view


GABA_GE@meta.data %>% 
  count(pca_snn_res.0.5,sex) %>%
  pivot_wider(names_from = sex,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Inh_",formatC(as.numeric(unfactor(pca_snn_res.0.5)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.5) %>% view


```

## Compositional Analysis

### Kip's method
Kip just used a glm model for this coupled with the compositions package
```{r}
pacman::p_load(compositions)

##He drops clusters with lower than 30 mean counts per cell, we may do that later

tmp_cluster_counts <- GABA_GE@meta.data %>% 
  count(pca_snn_res.0.5,sample) %>%
  pivot_wider(names_from = sample,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Inh_",formatC(as.numeric(unfactor(pca_snn_res.0.5)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.5) 

rownames(tmp_cluster_counts) <- tmp_cluster_counts$cluster
clr_cluster_counts <- as.data.frame(compositions::clr(t(tmp_cluster_counts %>% dplyr::select(-cluster))))
colnames(clr_cluster_counts) <- tmp_cluster_counts$cluster

clr_cluster_counts <- merge(sampleInfo, clr_cluster_counts, by.x = "name", by.y = "row.names")
clr_cluster_counts$sex <- as.factor(ifelse(clr_cluster_counts$sex == "F", 1, 0))
clr_cluster_counts$genotype <- as.factor(ifelse(clr_cluster_counts$genotype == "H", 1, 0))

diff_prop_models <- lapply(6:ncol(clr_cluster_counts),
                    function(x){
                      glm(clr_cluster_counts[,x] ~ sex + genotype,
                          data = clr_cluster_counts)
                    })

diff_prop <- lapply(diff_prop_models, summary)
diff_prop_coef <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,1]})))
diff_prop_std <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,2]})))
diff_prop_p <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,4]})))

diff_prop <- as.data.frame(do.call(cbind,list(colnames(clr_cluster_counts[,-1:-5]),diff_prop_coef,diff_prop_std,diff_prop_p)))

colnames(diff_prop) <- c("ID","Coef","StdErr","pvalue")
diff_prop$pvalue <- as.numeric(as.character(diff_prop$pvalue))
diff_prop_results <- diff_prop[order(diff_prop$pvalue),]
diff_prop_results$FDR <- p.adjust(diff_prop_results$pvalue, method = "fdr")
```
### Best practices
Recommended by the best practices article.  

#### scDC
This is linear and probably bad, doesnt do any compositional analysis

```{r}
#pacman::p_install_gh("SydneyBioX/scDC")
pacman::p_load(scDC)

test2 <- callr::r_bg(function(x) {
  scDC::scDC_noClustering(
    cellTypes = x$pca_snn_res.0.5,
    subject  = as.factor(x$id),
    ncores = 18,
    nboot = 50
  )},
  args = list(x = GABA_GE@meta.data)
  )

tmp_scDC <- test2$get_result()

results <- list()
for(i in unique(tmp_scDC$info$cellTypes)) {
  
  test <- cbind(tmp_scDC$info[,1:2], tmp_scDC$thetastar) %>%
    filter(cellTypes == i)
  
  test <- left_join(test, sampleInfo,
                    by = c("subject"= "id"))
  
  test2 <- lapply(1:ncol(tmp_scDC$nstar),
                  function(x) {
                    stats::glm(formula(paste0("`",x,"` ~ sex + genotype + sex:genotype")),
                               data = test)
                  }
  )
  
  
  ret <- mice::pool(test2,  dfcom = 13)
  ret <- summary(ret)
  ret$cluster <- i
  results[[i]] <- ret
}

diff_prop_id <- as.numeric(unlist(lapply(results, function(x){x[3,7]})))
diff_prop_coef <- as.numeric(unlist(lapply(results, function(x){x[3,2]})))
diff_prop_std <- as.numeric(unlist(lapply(results, function(x){x[3,3]})))
diff_prop_p <- as.numeric(unlist(lapply(results, function(x){x[3,6]})))

diff_prop <- as.data.frame(do.call(cbind,list(diff_prop_id,diff_prop_coef,diff_prop_std,diff_prop_p)))

colnames(diff_prop) <- c("id","Coef","StdErr","pvalue")
diff_prop$pvalue <- as.numeric(as.character(diff_prop$pvalue))
diff_prop_results <- diff_prop[order(diff_prop$pvalue),]
diff_prop_results$FDR <- p.adjust(diff_prop_results$pvalue, method = "fdr")

  

```

#### DCATS
```{r}
pacman::p_load("DCATS")

tmp_cluster_counts <- GABA_GE@meta.data %>% 
  count(pca_snn_res.0.5,sample) %>%
  pivot_wider(names_from = sample,
              values_from = n)

tmp_knn_mat <- knn_simMat(GABA_GE@graphs$pca_snn, GABA_GE@meta.data$pca_snn_res.0.5)



tmp_cluster_counts <- tmp_cluster_counts[match(rownames(tmp_knn_mat),
                                         tmp_cluster_counts$pca_snn_res.0.5,), ]
tmp_count_mat <- tmp_cluster_counts %>%
  dplyr::select(-pca_snn_res.0.5) %>%
  dplyr::mutate(across(everything(), function(x) replace_na(x, 0))) %>% 
  as.matrix %>% t
colnames(tmp_count_mat) <- rownames(tmp_knn_mat)

tmp_samples <- data_frame(name = colnames(tmp_cluster_counts[,2:15])) %>%
  left_join(sampleInfo)

tmp_cluster_counts <- table(GABA_GE@meta.data$id, GABA_GE@meta.data$pca_snn_res.0.5)

tmp_cluster_counts <- tmp_cluster_counts[, match(rownames(tmp_knn_mat),
                                         colnames(tmp_cluster_counts))]
tmp_count_mat <- tmp_cluster_counts

tmp_samples <- data_frame(id = rownames(tmp_count_mat)) %>%
  left_join(sampleInfo)

tmp_design <- data.frame(condition = tmp_samples$genotype,
                         sex = tmp_samples$sex)

test <- dcats_GLM(tmp_count_mat,
                  tmp_design,
                  similarity_mat = tmp_knn_mat,
                  base_model = "FULL"
                  )

test <- dcats_GLM(tmp_count_mat,
                  tmp_design,
                  similarity_mat = tmp_knn_mat,
                  base_model = "FULL")

cbind(rownames(tmp_knn_mat), test$fdr)
```

#### sccomp

```{r}
pacman::p_load(sccomp)

cmdstanr::check_cmdstan_toolchain(fix = TRUE) # Just checking system setting
cmdstanr::install_cmdstan()

test <- GABA_GE@meta.data %>%
  group_by(genotype, sex, id, pca_snn_res.0.5) %>%
  summarize(count = n()) %>% ungroup

sccomp_result <- 
  sccomp_estimate(
    .data = test,
    formula_composition = ~ genotype + sex, 
    .sample =  id, 
    .cell_group = pca_snn_res.0.5, 
    .count = count,
    bimodal_mean_variability_association = T,
    variational_inference = F,
    cores = 1  
  )

sccomp_result <- 
  sccomp_estimate(
    .data = GABA_GE,
    formula_composition = ~ genotype + sex, 
    .sample =  id, 
    .cell_group = pca_snn_res.0.5, 
    bimodal_mean_variability_association = T,
    variational_inference = F,
    cores = 1  
  )

sccomp_result <- 
  sccomp_estimate(
    .data = lognormal_GE,
    formula_composition = ~ genotype + sex, 
    .sample =  id, 
    .cell_group = cell_class_broad, 
    bimodal_mean_variability_association = T,
    variational_inference = F,
    cores = 1,
    max_sampling_iterations = 20000,
    mcmc_seed = saved_seed
  )


test<-sccomp_result %>% sccomp_test() 

test <- test %>%
  group_by(factor) %>%
  mutate(c_FDR = p.adjust(c_pH0, method = "fdr")) %>% ungroup

test$c_FDR <- p.adjust(test$c_pH0)

sccomp_result %>% 
  sccomp_test() %>%
  sccomp_proportional_fold_change(
    formula_composition = ~ genotype,
    from = "H",
    to = "L"
  )


plots <- sccomp_result %>% sccomp_test() %>% plot()

plots$boxplot[[2]]

plots$credible_intervals_1D
plots$credible_intervals_2D
```


##DE

### Combined
We want 3/4 from the methods.. so we will run the two fast ones first, as a DE gene has to be in one of those
First we run the fast ones, which are edgeR and DEseq2

#### Pseudo
```{r}
library(Libra)
plan(sequential)
gc()
options(future.globals.maxSize = 80 * 1024 ^ 3)
plan(callr, workers = 2)

DE <- list()

DE[["DESeq2_WALD"]] <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'DESeq2',
             de_type = 'Wald'
             )
  })

DE[["edgeR_QLF"]]  <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'edgeR',
             de_type = 'QLF'
             )
  })

plan(sequential)
gc()
```

Now we will process those

```{r}
GABA_DE <- lapply(which(resolved(DE)),
                  function(x) value(DE[[x]])) %>%
  bind_rows()

 GABA_DE <- lapply(which(resolved(DE)),
                  function(x) value(DE[[x]]))

GABA_DE[[1]] <- split(GABA_DE[[1]], GABA_DE[[1]]$cell_type)
GABA_DE[[2]] <- split(GABA_DE[[2]], GABA_DE[[2]]$cell_type)
 GABA_DE %>% filter(p_val_adj < 0.1) %>% 
       arrange(cell_type, gene) %>% View

```
Okay, so these do not do any of the default filtering on genes... probably not great.

Now we will run the ones that take longer..
#### DEsingle
```{r}
plan(sequential)
gc()
options(future.globals.maxSize = 20 * 1024 ^ 3)
plan(callr, workers = 10)

DE_single <- list()

for(i in unique(GABA_GE@meta.data$pca_snn_res.0.5)) {
  
  tmp_counts <- LayerData(GABA_GE, assay = "SCT", layer = "counts")
  tmp_cells <- GABA_GE@meta.data$pca_snn_res.0.5 == i
  
  tmp_geno <- GABA_GE@meta.data$genotype[tmp_cells] 
  tmp_geno <- factor(tmp_geno, levels = c("H","L"))
  
  tmp_genes <- GABA_DE %>%
    filter(cell_type == i) %>%
    group_by(de_method) %>%
    slice_min(p_val_adj, n = 100) %>%
    ungroup() %>%
    pull(gene) %>% unique
  
  tmp <- tmp_counts[which(rownames(tmp_counts) %in% tmp_genes), tmp_cells]
  
  DE_single[[i]] <- future({
    DEsingle::DEsingle(counts = tmp,
                       group = tmp_geno,
                       parallel = F)
  })
  
}
```


```{r}
 GABA_DE <- lapply(which(resolved(DE_single)),
                  function(x) value(DE_single[[x]]))
```

### Libra
```{r}
library(Libra)

plan(sequential)
gc()
options(future.globals.maxSize = 80 * 1024 ^ 3)
plan(callr, workers = 8)

DE <- list()

DE[["DESeq2_LRT"]] <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'DESeq2',
             de_type = 'LRT'
             )
  })

DE[["DESeq2_WALD"]] <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'DESeq2',
             de_type = 'Wald'
             )
  })

DE[["edgeR_LRT"]]  <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'edgeR',
             de_type = 'LRT'
             )
  })

DE[["edgeR_QLF"]]  <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'edgeR',
             de_type = 'QLF'
             )
  })

DE[["limma_trend"]] <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'limma',
             de_type = 'trend'
             )
  })

DE[["limma_voom"]] <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'pseudobulk',
             de_method = 'limma',
             de_type = 'voom'
             ) 
  })


DE[["MAST"]] <- future({
  run_de(GABA_GE,
             cell_type_col = "pca_snn_res.0.5",
             label_col = "genotype",
             replicate_col = "id",
             de_family = 'singlecell',
             de_method = 'MAST'
             )
  })

plan(sequential)
gc()
```

```{r}
which(resolved(DE))

tmp <- lapply(which(resolved(DE)), function(x) value(DE[[x]]))

tmp2 <- bind_rows(tmp, .id = "method") 

tmp2 %>%
  filter(p_val_adj < 0.5) %>%
  count(method)

View(tmp2 %>%
       filter(p_val_adj < 0.05,
              method %in% c("DESeq2_WALD", "edgeR_QLF", "limma_voom", "MAST")) %>% 
       arrange(cell_type, gene))


tmp2 %>% filter(str_starts(gene, "Hbb")) %>% View
```


### DEsingle
```{r}
tmp <- GABA_GE@meta.data$pca_snn_res.0.5 == "1"
tmp2 <- GABA_GE@meta.data$genotype[tmp] 
tmp2 <- factor(tmp2, levels = c("H","L"))
tmp <- GABA_GE@assays[["SCT"]]$counts[, tmp]

DE5 <- DEsingle::DEsingle(counts = tmp,
                          group = tmp2,
                          parallel = T,
                          BPPARAM = BiocParallel::MulticoreParam(2, RNGseed=saved_seed))

DE5_type <-  DEsingle::DEtype(results = DE5, threshold = 0.05)
```

### Kip's Method

```{r}
#Kip only looks at clusters where the mean count of cells per sample is above 30, I think we will look at all clusters to start.
#He also drops animals which have fewer than 10 cells in a given cluster.
#He uses raw RNA counts, not SCT... so does Libra... but he converts it to lognorm+1 later

tmp_counts <- LayerData(GABA_GE, assay = "SCT", layer = "counts")
tmp_counts <- as.matrix(tmp_counts)
tmp_counts <- as.data.frame(tmp_counts)

i <- "10"

#This is Kip's animal filter, 10 or more cells per animal
tmp_cells <- GABA_GE@meta.data %>%
  rownames_to_column(var = "cells_new") %>%
  filter(pca_snn_res.0.5 == i) %>%
  group_by(id) %>%
  filter(n() >= 10 ) %>%
  column_to_rownames(var = "cells_new")

tmp_counts_2 <- tmp_counts[, which(colnames(tmp_counts) %in% rownames(tmp_cells))]

#Kip filters genes by using a mean number of counts of 5 across all cells.. which is very aggressive. I think I want to maybe do something using the number of non-zero cells for each gene following our logic early borrowed from SDC

tmp_counts_3 <- t(tmp_counts_2) > 0
tmp_counts_3 <- merge(tmp_cells, tmp_counts_3, by.x = "row.names", by.y = "row.names")

#count the number of non-zero counts by ID
tmp_counts_4 <- tmp_counts_3 %>%
  group_by(sex,genotype,id, .drop = F) %>%
  summarize(across(42:ncol(tmp_counts_3)-3, sum)) %>%
  ungroup

#Compute the average number of non-zero counts by group (sex by genotype)
tmp_counts_5 <- tmp_counts_4 %>%
  group_by(sex,genotype, .drop = F) %>% 
  summarize(across(4:ncol(tmp_counts_4)-2, mean)) %>%
  ungroup

#If any group has lower than a mean of 2% of cells in the cluster (or 25) with non-zero cells per sample, lets drop that gene
tmp_min <- max(ncol(tmp_counts_2)/100*2, 25)


tmp_counts_6 <- tmp_counts_5 %>%
  summarize(across(3:ncol(tmp_counts_5), function(x) any(x<tmp_min)))

tmp_keep <- colnames(tmp_counts_6)[which(!tmp_counts_6[1,])]
tmp_counts_2 <- tmp_counts_2[which(rownames(tmp_counts_2) %in% tmp_keep), ]
tmp_counts_2 <- log2(tmp_counts_2 + 1)
rm(tmp_counts_3, tmp_counts_4, tmp_counts_5, tmp_counts_6)


sca <- suppressMessages(MAST::FromMatrix(exprsArray=as.matrix(tmp_counts_2),
                                         cData=tmp_cells,
                                         fData=data.frame(genes=rownames(tmp_counts_2))
                                         )
                        )
cdr2 <- colSums(SummarizedExperiment::assay(sca)>0)
SummarizedExperiment::colData(sca)$ngeneson <- scale(cdr2)
SummarizedExperiment::colData(sca)$Sex <- as.factor(SummarizedExperiment::colData(sca)$sex)
SummarizedExperiment::colData(sca)$AnimalID <- as.factor(SummarizedExperiment::colData(sca)$id)
SummarizedExperiment::colData(sca)$Group <- as.factor(SummarizedExperiment::colData(sca)$genotype)
SummarizedExperiment::colData(sca)$HighPreference <- ifelse(as.factor(SummarizedExperiment::colData(sca)$Group) == "L", 0, 1)
sprintf("MAST Two-part hurdle GLMM - Cluster %s:\nSamples: %d(%s)\nGenes: %d\nCells: %d",
        i,
        length(unique(SummarizedExperiment::colData(sca)$id)),
        paste0(unique(SummarizedExperiment::colData(sca)$id), collapse = ", "),
        nrow(tmp_counts_2),
        ncol(tmp_counts_2)
        ) %>% cat
```

```{r}
pacman::p_load(progressr)
pacman::p_load(future.callr)
```

```{r}
plan(sequential)
gc()
options(future.globals.maxSize = 80 * 1024 ^ 3)
#We use two here because the underlying method appears to use a full socket by default, so if we spawn 2 it will use 100% of rosalind
plan(callr, workers = 10)

handlers(  handler_progress(
    format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta",
    width    = 60,
    complete = "+")
)
```

```{r}
#Job chunking
tmp_n <- nrow(tmp_counts_2)

tmp_jobs <- sample(1:(tmp_n))
tmp_chunk <- 10
worker_jobs <- split(tmp_jobs, ceiling(seq_len(length(tmp_jobs))/tmp_chunk))
njobs <- length(tmp_jobs)
nchunks <- length(worker_jobs)

#Job lists

worker_futures <- list()
worker_info <- list()
tmp_finished <- list()
tmp_results <- list()
tmp_jobid <- 1
tmp_genesprocessed <- 0

#Save frequency (percentages)
tmp_save_freq <- 0

with_progress({
  
  if(tmp_save_freq > 0)
    tmp_save <- tmp_save_freq*floor(njobs/100)
  
  p <- progressor(steps = nchunks)
  
  
  jobLoop <- function(job_i) {
    ret <- MAST::zlm(~ ngeneson + Sex + Group + (1 | AnimalID),
                     sca[job_i, ],
                     method='glmer',
                     ebayes = F,
                     strictConvergence = FALSE,
                     parallel = T)
    ret <- MAST::summary(ret,
                         doLRT='GroupL')
    p()
    return(ret)
  }
  
  message(paste0("Sending out initial jobs to ", nbrOfWorkers() , " workers.\n"))
  p(amount = 0)
  
  for(worker_i in seq_len(nbrOfWorkers())) {
    
    tmp_worker_info <- list(id = as.character(worker_i),
                            job = worker_jobs[[tmp_jobid]],
                            jobid = tmp_jobid)
    
    worker_futures[[tmp_worker_info$id]] <- future({
      jobLoop(tmp_worker_info$job)
    })
    
    tmp_jobid <- tmp_jobid + 1
    
    worker_info[[tmp_worker_info$id]] <- tmp_worker_info
  }
  
  message(paste0("Waiting for jobs to finish.\n"))
  p(amount = 0)
  
  repeat {
    #Check for finished workers
    tmp_resolved <- which(resolved(worker_futures))
    
    if(length(tmp_resolved)) {
      
      for(worker_i in tmp_resolved) {
        
        #only assign jobs if we still have one
        if(tmp_jobid > nchunks) break;
        
        tmp_worker_info <- worker_info[[names(worker_futures)[worker_i]]]
        
        message(paste0("Worker ", tmp_worker_info$id,
                       " finished their chunk, sending a new one (",
                       nchunks - tmp_jobid, " remain).\n"))
        p(amount = 0)
        
        tmp_finished[[as.character(tmp_worker_info$jobid)]] <- worker_futures[[tmp_worker_info$id]]
        tmp_results[[as.character(tmp_worker_info$jobid)]] <- value(worker_futures[[tmp_worker_info$id]])
        
        tmp_worker_info$job <- worker_jobs[[tmp_jobid]]
        tmp_worker_info$jobid <- tmp_jobid
        
        worker_futures[[tmp_worker_info$id]] <- future({
          jobLoop(tmp_worker_info$job)
        })
        worker_info[[tmp_worker_info$id]] <- tmp_worker_info
        
        tmp_jobid <- tmp_jobid + 1
        
      }
      
      if(tmp_jobid > nchunks) break;
    } else {
      
      #We could also in theory process data here if it is a slower process
      Sys.sleep(1)
    }
  }
  message("All out of chunks, waiting for remaining workers to finish.\n")
  p(amount = 0)
  
  repeat {
    if(all(resolved(worker_futures))) {
      
       for(worker_i in 1:length(worker_futures)) {
         tmp_worker_info <- worker_info[[names(worker_futures)[worker_i]]]
         tmp_finished[[as.character(tmp_worker_info$jobid)]] <- worker_futures[[tmp_worker_info$id]] 
         tmp_results[[as.character(tmp_worker_info$jobid)]] <- value(worker_futures[[tmp_worker_info$id]])
       }
       message(paste0("We are done.\n"))
       p(amount = 0)
      break;
    }
    else {Sys.sleep(1)}
  }
})
```


```{r}
save(worker_jobs, tmp_finished, file = "DE_1_test.Rdata")
```

```{r}
tmp <- lapply(tmp_finished, function(x) {
  ret <- value(x)
  ret <- ret$datatable
  ret <- merge(ret[ret$contrast=='GroupL'
                                  & ret$component=='logFC', c(1,7,5,6,8)],
                        ret[ret$contrast=='GroupL'
                                  & ret$component=='H', c(1,4)],
                    by = 'primerid')
  }
  ) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(`Pr(>Chisq)`, "fdr", n = 7000))

```



# GLUT
```{r}
GLUT_GE <- subset(filtered_GE,
                  subset = (id != "HS_3") &
                    (id != "HS_16"))

test4 <- lognormal_GE@meta.data %>%
  dplyr::select(cells, cell_class, cell_class_broad, cell_doublet, cluster_class, cluster_class_broad, cluster_doublet) 
test4$cells <- str_split_i(test4$cells, "[0-9]_", 2)               
test4 <- left_join(GLUT_GE@meta.data, test4)
rownames(test4) <- test4$cells
GLUT_GE@meta.data <- test4

#Same logic as above for GABA
GLUT_GE <- subset(GLUT_GE,
                  subset = 
                    ((cell_class == "glut") |
                       (cluster_class == "glut"))  &
                    #(cell_doublet == F) &
                    (cell_class_broad != "non") &
                    (cell_class != "gaba"))

GLUT_GE

rm(test,test2,test3,test4, tmp, tmp_markers, tmp_meta, tmp_plot_1, tmp_plot_2)
```

We are left with 12068 cells.  We will filter by low expression again, though.

### Cell Filtering
```{r}
counts <- GetAssayData(object = GLUT_GE, layer = "counts")

nonzeros <- Matrix::rowSums(counts>2)

nonzeros <- data_frame(gene = rownames(counts),
                       ncells = nonzeros)

nonzeros <- nonzeros %>%
arrange(ncells) %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) 

sum(nonzeros$ncells > 10)

nonzeros %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) %>%
  ggplot() +
  geom_point(aes(y = idx, x = ncells)) +
  theme_JQA() +
    coord_cartesian(xlim=c(0,10000))

nonzeros[c(
  which.max(nonzeros$idx>10),
  which.max(nonzeros$idx>20),
  which.max(nonzeros$idx>30),
  which.max(nonzeros$idx>40),
  which.max(nonzeros$idx>50),
  which.max(nonzeros$idx>60),
  which.max(nonzeros$idx>70),
  which.max(nonzeros$idx>80),
  which.max(nonzeros$idx>90)), ]


#This leaves us with ~ 14k genes which is great  
to_drop <- nonzeros$ncells < 10

sum(!to_drop)
#that leaves 13k genes, which is fine I think

GLUT_GE <- subset(x = GLUT_GE,
                      features = nonzeros$gene[!to_drop])
GLUT_GE 
rm(counts,nonzeros,to_drop)
```

12068 glut Neurons, 11804 genes

### Normalization
```{r}
GLUT_GE <- SplitObject(GLUT_GE, split.by = "id")
```

```{r}
plan(sequential)
gc()
plan(multisession, workers = 7)
options(future.globals.maxSize = 80 * 1024 ^ 3)
GLUT_GE <- future.apply::future_lapply(GLUT_GE,
                                      function(x) {
                                        SCTransform(x,
                                                    assay = "RNA",
                                                    vars.to.regress = "percent.mt",
                                                    verbose = FALSE,
                                                    seed.use = saved_seed,
                                                    )
                                      },
                                      future.seed = saved_seed
) 

GLUT_GE <- future.apply::future_lapply(GLUT_GE,
                                            function(x) {
                                              NormalizeData(x,
                                                            assay = "RNA",
                                              )
                                            },
                                            future.seed = saved_seed
) 

GLUT_GE <- merge(GLUT_GE[[1]],
              GLUT_GE[2:14],
              add.cell.ids = names(GLUT_GE))

GLUT_GE[["RNA"]] <- JoinLayers(GLUT_GE[["RNA"]])

VariableFeatures(GLUT_GE, assay = "SCT") <- rownames(GLUT_GE[["SCT"]]@scale.data)

# lognormal_GE <- NormalizeData(subset(filtered_GE,
#                                      subset = (id != "HS_3") &
#                                        (id != "HS_16")))
# 
# lognormal_GE <- FindVariableFeatures(lognormal_GE)
# 
# lognormal_GE <- ScaleData(lognormal_GE)
plan(sequential)
gc()
```

```{r}
GLUT_GE <- RunPCA(GLUT_GE,
                       assay = "SCT")
```


```{r}
ElbowPlot(GLUT_GE,ndims = 40)
```

```{r}
DimHeatmap(GLUT_GE,dims=1:10, cells=500)
DimHeatmap(GLUT_GE,dims=11:20, cells=500)
DimHeatmap(GLUT_GE,dims=21:30, cells=500)
DimHeatmap(GLUT_GE,dims=31:40, cells=500)
```

```{r}
GLUT_GE <- RunUMAP(GLUT_GE, 
                        dims = 1:40,
                        reduction = "pca",
                        reduction.name = "umap",
                        seed.use = saved_seed)
```


### Clustering
```{r}
GLUT_GE <- FindNeighbors(GLUT_GE,
                              dims = 1:40,
                              reduction = "pca",
                              graph.name = paste0("pca_", c("nn","snn")))
```

```{r}
library(reticulate)
numpy <- reticulate::import("numpy")
pandas <- reticulate::import("pandas")
leidenalg <- reticulate::import("leidenalg")

plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)

GLUT_GE <- FindClusters(GLUT_GE,
                             method  = "igraph",
                             algorithm = 4,
                             random.seed = saved_seed,
                             graph.name = "pca_snn",
                             resolution = c(0.01,0.05,0.1,0.25,0.5,0.8, 1, 1.5))

plan(sequential)
gc()
```

```{r}
clustree(GLUT_GE, prefix = "pca_snn_res.")
```

Another clean clustering, similar numbers of distinct clusters at each resolution as we see in the gabaergic.. I like 0.25 this time

```{r}
DimPlot(GLUT_GE,
        group.by = c("genotype",
                     "id",
                     "class",
                     "GEX_lane",
                     "pca_snn_res.0.01",
                     "pca_snn_res.0.05",
                     "pca_snn_res.0.1",
                     "pca_snn_res.0.25",
                     "pca_snn_res.0.5"),
        label = T)
```

There are some genoype/lane effects here driven by individual samples - unlike in the gaba.. can we justify doing an integration here and not on gaba?  Or is it easy enough to say just throw away clusters which are dominated by a single sample (and that they are likely the result of disection)

```{r}
ggplot(GLUT_GE@meta.data) +
    geom_bar(aes(x=pca_snn_res.0.25, fill=(cell_class)), position=position_fill()) + theme_JQA()

```

```{r}
GLUT_GE@meta.data %>% dplyr::count(pca_snn_res.0.25)
```

No clusters below 100

```{r}
DimPlot(GLUT_GE,
        group.by = c("pca_snn_res.0.25"),
        label = T)
```

```{r}
GLUT_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.25, fill=(class)), position=position_fill())+ theme_JQA()
#really only two mixed clusters on cell_class

GLUT_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.25, fill=(id)), position=position_fill())+ theme_JQA()
#12 and 9 are completely dominated by hs_10.. 13 is only a few (5)..

GLUT_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.25, fill=(genotype)), position=position_fill())+ theme_JQA()

#We have two fewer HP mice now.

GLUT_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.25, fill=(sex)), position=position_fill())+ theme_JQA()
```

```{r}
Idents(object = GLUT_GE) <- "pca_snn_res.0.25"
```

```{r}
DimPlot(GLUT_GE,
        group.by = "pca_snn_res.0.25",
        label = T)
```

```{r}
test <- GLUT_GE@meta.data

test2 <- sapply(unique(test$id), function(sample) {
  sapply(as.numeric(levels(test$pca_snn_res.0.25)), function(cluster) {
  fisher.test(test$id == sample,
              test$pca_snn_res.0.25 == cluster,
              alternative = "g")$p.value
  })
})

test3 <- -log10(test2)
test3[test3 == Inf] <- 100
test3[test3 == NA] <- 0
test3[test3 > 25] <- 25
rownames(test3) <- as.numeric(levels(test$pca_snn_res.0.25))

test4 <- data_frame(id = colnames(test3)) %>%
  left_join(sampleInfo)

Heatmap(test3,
        cluster_rows = T,
        row_title_side = "right",
        top_annotation = HeatmapAnnotation(Sex = test4$sex,
                                           Genotype = test4$genotype),
        col = colorRamp2(c(0, 1.3,
                     3.475, 5.650, 7.825, 10, 25), c("black", "white",
                                                 "#00B785","#53CC67","#B2DC3C","#FDE333", "red")))

```
Some clusters are only from HS_12 and HS_10.. see 10 and 5 in partciular

```{r}
FeaturePlot(
  GLUT_GE,
  reduction = "umap",
  features = c("nFeature_RNA",
               "nCount_RNA",
               "nCount_SCT",
               "nFeature_SCT",
               "percent.hge",
               "percent.mt",
               "percent.rp"),
  order = TRUE,
  pt.size = 0.4,
  min.cutoff = 'q10',
  label = T
)
```


```{r}
GLUT_GE@meta.data %>%
  dplyr::count(pca_snn_res.0.25)

VlnPlot(GLUT_GE, features = "nCount_RNA", group.by = "pca_snn_res.0.25", assay = "RNA") /
VlnPlot(GLUT_GE, features = "nFeature_RNA", group.by = "pca_snn_res.0.25", assay = "RNA")/
VlnPlot(GLUT_GE, features = "nCount_SCT", group.by = "pca_snn_res.0.25", assay = "RNA") /
VlnPlot(GLUT_GE, features = "nFeature_SCT", group.by = "pca_snn_res.0.25", assay = "RNA")  
#15 and 21 seem to have double the feature count..

VlnPlot(GLUT_GE,
        assay = "SCT",
        layer = "data",
        features = c("Gad1", "Gad2", "Slc32a1", "Slc17a7", "Slc17a6"), group.by = "pca_snn_res.0.25", ncol = 2)
```
10 and 13 are weird (driven by a single sample)


```{r}
GLUT_GE_orig <- GLUT_GE
```

```{r}
plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)
GLUT_GE <- PrepSCTFindMarkers(GLUT_GE)
plan(sequential)
gc()

#9302, ~2x that (almost) of the gaba cells
```

## Cluster DE Markers
```{r}
Idents(object = GLUT_GE) <- "pca_snn_res.0.25"

plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)

tmp_markers_glut <- future.apply::future_lapply(1:13,
                             function(i)
                             {
                               SeuratWrappers::RunPresto(GLUT_GE,
                                                         assay = "SCT",
                                                         ident.1 = i)
                             },
                             future.seed = saved_seed
                             )

plan(sequential)
gc()


tmp_markers_glut <- lapply(tmp_markers_glut, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})

#write.xlsx(tmp_markers_gaba, file = "GABAmarkers_jqa_103124.xlsx")
```

#### Dendogram
```{r}
tmp_tree <- BuildClusterTree(GLUT_GE,
                          assay = "SCT",
                          slot  = "scale.data",
                          dim = 1:40
                          )

PlotClusterTree(tmp_tree, 
                direction = "rightwards")
```

```{r}
GLUT_GE@meta.data$pca_snn_res.0.25 %>% table
```


#### Dotplot markers

#### Dotplot doublets
Let's check for doublets against our known types, since we didn't drop any earlier
```{r}
(DotPlot(GLUT_GE,
        features = c(tmp_gaba, tmp_glut, tmp_non),
        group.by = c("pca_snn_res.0.25"),
        scale = F,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))) %>%
  plotly::ggplotly()

#18 is gaba+glut
#19 is gaba+opc
#17 is a little noisy, maybe olig

DimPlot(GLUT_GE,
        split.by = "pca_snn_res.0.25",
        reduction = "umap",
        ncol = 5)
```


```{r}
tmp <- tmp_markers_glut %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

#the number of markers we want per cluster
tmp2 <- tmp %>%
  #mutate(cluster = if_else(cluster %in% c("6", "9", "10"), "6", cluster)) %>%
  group_by(cluster) %>%
  slice_max(rank, n = 10) %>%
  #dplyr::filter(rank > 2.8) %>%
  #dplyr::select(-n) %>%
  ungroup()
  
tmp2 %>% count(cluster) %>% arrange(n)

(DotPlot(GLUT_GE,
        features = unique(tmp2$gene),
        group.by = c("pca_snn_res.0.25"),
        scale = T,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))) %>%
  plotly::ggplotly()
```
We can merge 13,1 and 6... but clusters 10 and 13 are completely dominated by a single sample.. so we are just going to get rid of them, 11 for being dominated by 2, as well.  These are all smaller, outlier clusters.

```{r}
GLUT_GE <- subset(GLUT_GE, 
                    subset = (pca_snn_res.0.25 != "11") &
                      (pca_snn_res.0.25 != "10") &
                      (pca_snn_res.0.25 != "13")
                      )
```

```{r}
tmp <- GLUT_GE@meta.data

tmp <- tmp %>%
    mutate(pca_snn_res.0.25 = if_else(pca_snn_res.0.25 %in% c("6"), "1", pca_snn_res.0.25)
           )

tmp2 <- tmp %>%
  count(pca_snn_res.0.25) %>%
  arrange(desc(n)) %>%
  mutate(new = factor(row_number()))

tmp2

tmp <- left_join(tmp, tmp2)
tmp$pca_snn_res.0.25 <- tmp$new
tmp <- tmp %>%
  dplyr::select(-new, -n)
rownames(tmp) <- GLUT_GE[["SCT"]]@data@Dimnames[[2]]

GLUT_GE@meta.data <- tmp

#Do we need to do this?
Idents(object = GLUT_GE) <- "pca_snn_res.0.25"
```

Now lets go re-run the cluster thing, we have to set recorrect_umi to false here, since we subset a SCT stack that we previoulsly ran prepsctmarkers on

## Cluster DE Markers 2
```{r}
plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)

tmp_markers_glut <- future.apply::future_lapply(1:9,
                             function(i)
                             {
                               SeuratWrappers::RunPresto(GLUT_GE,
                                                         assay = "SCT",
                                                         ident.1 = i,
                                                         recorrect_umi = F)
                             },
                             future.seed = saved_seed
                             )

plan(sequential)
gc()

tmp_markers_glut <- lapply(tmp_markers_glut, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})

#write.xlsx(tmp_markers_gaba, file = "GABAmarkers_jqa_103124.xlsx")
```

#### Dendogram 2
```{r}
GLUT_GE <- BuildClusterTree(GLUT_GE,
                          assay = "SCT",
                          slot  = "data",
                          dim = 1:40
                          )

PlotClusterTree(GLUT_GE) 

plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)
tree_markers_glut <- future.apply::future_lapply(
  GLUT_GE@tools$BuildClusterTree[["edge"]][,1] %>% unique(),
  function(i) {
    SeuratWrappers::RunPresto(GLUT_GE,
                              assay = "SCT",
                              ident.1 = GLUT_GE@tools$BuildClusterTree,
                              ident.2 = i,
                              recorrect_umi = F)
  },
  future.seed = saved_seed
)
plan(sequential)
gc()

names(tree_markers_glut) <- GLUT_GE@tools$BuildClusterTree[["edge"]][,1] %>% unique()

tree_markers_glut <- lapply(tree_markers_glut, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})


#  library(ggtree)
# ggtree(tmp_tree@tools$BuildClusterTree,
#           branch.length="none") + 
#      geom_tiplab(as_ylab = F) + 
#      geom_nodelab()+
#      geom_rootedge(rootedge=1) +
#   geom_text(data = data_frame(x = 1:10, y = 1:10, label = 1:10), aes(x=x,y=y,label=label))
# 
# 
#      #geom_cladelab(node = 19, label = "test", align = T, hjust = -10   theme_tree()

```

```{r}
write.xlsx(tmp_markers_glut, file = "GLUTmarkers_jqa_121924.xlsx")
write.xlsx(tree_markers_glut, file = "GLUTtreemarkers_jqa_121924.xlsx")
```

#### Dotplot markers 2


```{r}
tmp_levels <- GLUT_GE@tools[["BuildClusterTree"]][["edge"]][GLUT_GE@tools[["BuildClusterTree"]][["edge"]][, 2]<10, 2]

#tmp_levels <- GLUT_GE@tools[["BuildClusterTree"]]$tip.label

GLUT_GE@meta.data$pca_snn_res.0.25 <- factor(GLUT_GE@meta.data$pca_snn_res.0.25,
                                            levels = as.character(tmp_levels))

Idents(GLUT_GE) <- "pca_snn_res.0.25"
```

##### Top markers, cluster
```{r}
tmp <- tmp_markers_glut %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

tmp3 <- tmp %>%
  mutate(cluster = factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>% 
  group_by(cluster) %>%
  slice_max(rank, n = 2) 

tmp_plot_1 <- wrap_plots(
  
  ggtree::ggtree(GLUT_GE@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree(),
  
  DotPlot(GLUT_GE,
          features = (unique(tmp3$gene)),
          group.by = c("pca_snn_res.0.25"),
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)

tmp_plot_1

```

##### Top markers, tree
```{r}
tmp <- tree_markers_glut %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

#tmp2 <- tmp2[!(tmp2$gene %in% unique(tmp2$gene[duplicated(tmp2$gene)])), ]

tmp3 <- bind_rows(
  tmp %>%
  group_by(cluster) %>%
  slice_max(rank, n = 3),
    tmp %>%
  group_by(cluster) %>%
  slice_min(rank, n = 3)
) %>% arrange(cluster,rank)

tmp_plot_1 <- wrap_plots(
  
  ggtree::ggtree(GLUT_GE@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree(),
  
  DotPlot(GLUT_GE,
          features = (unique(tmp3$gene)),
          group.by = c("pca_snn_res.0.25"),
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)

tmp_plot_1

```

##### BICCN Neuropiptides
https://www.nature.com/articles/s41586-023-06812-z#Sec10

```{r}
tmp <- tmp_markers_glut %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

tmp_features <- c(
  "Adcyap1",	"Npff",
"Adm",	"Nps",
"Agrp",	"Npvf",
"Agt",	"Npw",
"Apln",	"Npy",
"Avp",	"Nts",
"Calca",	"Oxt",
"Calcb",	"Pdyn",
"Cartpt",	"Penk",
"Cck",	"Pmch",
"Cort",	"Pnoc",
"Crh",	"Pomc",
"Edn1",	"Prlh",
"Edn3",	"Pthlh",
"Gal",	"Pyy",
"Gcg",	"Rln3",
"Ghrh",	"Sst",
"Gnrh1",	"Tac1",
"Grp",	"Tac2",
"Hcrt",	"Trh",
"Kiss1",	"Ucn",
"Nmb",	"Ucn3",
"Nms",	"Uts2b",
"Nmu",	"Vip",
"Npb"
)

#tmp_features <- chosen_markers
#tmp_features <- tmp_features[tmp_features %in% rownames(GABA_GE)]

tmp_features <- tmp %>%
  dplyr::filter(gene %in% tmp_features) %>%
  group_by(gene) %>%
  slice_max(rank, n = 1) %>%
  dplyr::mutate(cluster =
                  factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>%
  arrange(cluster) %>%
  pull(gene)


DotPlot(GLUT_GE,
        features = unique(tmp_features),
        group.by = c("pca_snn_res.0.25"),
        scale = T,
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))

```

##Comparing to HOchgerner
```{r}

#The vglut1 clusters..
hoch <- list()

tmp <- tmp_markers_glut %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Slc17a7",
                                             "Zic4","Trp73","Rspo2","Sema3e",
                                             "Sema5a","Dcn","Wfs1","Sorc3",
                                             "Cbln1","Coch","Thrsp","Lamp5",
                                             "Bok","Pamr1","Sim1","Ciql1","Gpx3",
                                             "Rxfp1","Ermn","Mpped1"))) %>%
  bind_rows(.id = "cluster") 

hoch[["A"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))

view(tmp)
view(hoch[["A"]])

tmp <-tmp_markers_glut %>%
  lapply(function(x) x %>% dplyr::filter(gene %in%  
                                           c("Slc17a7",
                                             "Cd36","Calb2","Gpr101","Grem1",
                                             "Dcn","C1ql2","Cartpt","Fam46a",
                                             "Fmo1","Rxfp3","Il33","Fbln1","Eps8",
                                             "Cd44", "Fibcd1", "Vit", "Car12",
                                             "Celsr1", "Trh", "Mdga1"
                                            ))) %>%
  bind_rows(.id = "cluster") 

hoch[["B"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))

view(tmp)
view(hoch[["B"]])

tmp <- tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Slc17a7",
                                             "Trh", "Rxfp1", "Medag", "Kit",
                                             "Slc23a", "Plxd3", "Reln",
                                             "Tac1","Igfbp5","Oasl2","Ifit1",
                                             "Grp","Cpne8", "St8sia2", "Cald1",
                                             "Gsg1l", "Prox1", "Igfn1", "Ndst5",
                                             "Mid1","Cdh22"))) %>%
  bind_rows(.id = "cluster")


hoch[["C"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))
view(tmp)
view(hoch[["C"]])

tmp <- tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Lhx6", "Tshz2", "Stab1", "Prl1",
                                             "Nxph1", "Sox6", "Moxd1", "Sst", "Pvalb",
                                             "Lhx8", "Th", "Cbln4", "Luzp2", "Nxph2",
                                             "Prlr", "Greb1", "Calcr", "Tac1", "St18",
                                             "Satb1", "Chodl", "Fign", "Npy", "Tmtc4",
                                             "Nek7", "Rpb4", "Vwc2", "Crabp1",
                                             "Etv1", "Pthlh"))) %>%
  bind_rows(.id = "cluster") 


hoch[["D"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))
view(tmp)
view(hoch[["D"]])


tmp_markers_gaba[[1]] %>%
  dplyr::filter(
    gene %in% c("Foxp2", "Fmod",
                "Adra2a",
                "Col6a1",
                "Htr1f",
                "Pax6") 
  )

test <- tmp_markers_gaba %>%lapply(function(x) {
  x %>% dplyr::filter(pct.1 > .75)}
  )

```

##### UMAP
```{r}
test <- GLUT_GE@reductions$umap@cell.embeddings
test2 <- GLUT_GE@meta.data

all.equal(rownames(test),
          rownames(test2))

test <- bind_cols(test,test2)

test2 <- test %>%
  group_by(pca_snn_res.0.25) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(cell_color = factor(paste0("Exc_",formatC(row_number(),
                                                           format = "d", flag = "0", digits = 1))))

test <- left_join(test, test2 %>% dplyr::select(cell_color, pca_snn_res.0.25))

tmp_plot_2 <- ggplot(test %>%
         group_by(pca_snn_res.0.25) %>%
         slice_sample(prop = 0.25)) +
  geom_point(aes(x = umap_1,
                 y = umap_2,
                 color = cell_color),
             shape = 16) +
  ggrepel::geom_text_repel(data = test %>%
                     group_by(pca_snn_res.0.25) %>%
                     summarize(x = median(umap_1),
                               y = median(umap_2),
                               cell_color = unique(cell_color)),
                   aes(x = x, y = y, label = cell_color),
                  size = 5,
                  min.segment.length = 0.1,
                  max.overlaps = 15,
                  box.padding = 2,
                  color = "grey25") +
  theme_JQA() +
  scale_color_manual(
    values = c(
      hcl.colors(9,palette="Dark3")
    ),
  )  +
  xlab("Umap Dimension 1 (arb)") +
  ylab("Umap Dimension 2 (arb)") +
  theme(legend.position = "none") + ggtitle("Excitatory Neuronal Clusters")

tmp_plot_2

```

```{r}
(tmp_plot_2 + tmp_plot_1) + plot_annotation(tag_levels = "A")
```


```{r}
tmp_file = "Glut_seurat_v3.RData"

if(!file.exists(tmp_file)) {
  save(GLUT_GE, tree_markers_glut, tmp_markers_glut, file = tmp_file)
} else {
  load(tmp_file)
}

rm(tmp_file)
```

## Supplementary tables

```{r}
GLUT_GE@meta.data %>% 
  count(pca_snn_res.0.25,sample) %>%
  pivot_wider(names_from = sample,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Exc_",formatC(as.numeric(unfactor(pca_snn_res.0.25)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.25) %>% view

GLUT_GE@meta.data %>% 
  count(pca_snn_res.0.25,genotype) %>%
  pivot_wider(names_from = genotype,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Exc_",formatC(as.numeric(unfactor(pca_snn_res.0.25)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.25) %>% view


GLUT_GE@meta.data %>% 
  count(pca_snn_res.0.25,sex) %>%
  pivot_wider(names_from = sex,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Exc_",formatC(as.numeric(unfactor(pca_snn_res.0.25)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.25) %>% view


```
## Compositional Analysis

### Kip's method
Kip just used a glm model for this coupled with the compositions package
```{r}
pacman::p_load(compositions)

##He drops clusters with lower than 30 mean counts per cell, we may do that later

tmp_cluster_counts <- GLUT_GE@meta.data %>% 
  count(pca_snn_res.0.25,sample) %>%
  pivot_wider(names_from = sample,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(paste0("Inh_",formatC(as.numeric(unfactor(pca_snn_res.0.25)),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  dplyr::select(-pca_snn_res.0.25) 

rownames(tmp_cluster_counts) <- tmp_cluster_counts$cluster
clr_cluster_counts <- as.data.frame(compositions::clr(t(tmp_cluster_counts %>% dplyr::select(-cluster))))
colnames(clr_cluster_counts) <- tmp_cluster_counts$cluster

clr_cluster_counts <- merge(sampleInfo, clr_cluster_counts, by.x = "name", by.y = "row.names")
clr_cluster_counts$sex <- as.factor(ifelse(clr_cluster_counts$sex == "F", 1, 0))
clr_cluster_counts$genotype <- as.factor(ifelse(clr_cluster_counts$genotype == "H", 1, 0))

diff_prop_models <- lapply(6:ncol(clr_cluster_counts),
                    function(x){
                      glm(clr_cluster_counts[,x] ~ sex + genotype,
                          data = clr_cluster_counts)
                    })

diff_prop <- lapply(diff_prop_models, summary)
diff_prop_coef <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,1]})))
diff_prop_std <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,2]})))
diff_prop_p <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,4]})))

diff_prop <- as.data.frame(do.call(cbind,list(colnames(clr_cluster_counts[,-1:-5]),diff_prop_coef,diff_prop_std,diff_prop_p)))

colnames(diff_prop) <- c("ID","Coef","StdErr","pvalue")
diff_prop$pvalue <- as.numeric(as.character(diff_prop$pvalue))
diff_prop_results <- diff_prop[order(diff_prop$pvalue),]
diff_prop_results$FDR <- p.adjust(diff_prop_results$pvalue, method = "fdr")
```


# NON
```{r}
NON_GE <- subset(filtered_GE,
                  subset = (id != "HS_3") &
                    (id != "HS_16"))

test4 <- lognormal_GE@meta.data %>%
  dplyr::select(cells, cell_class, cell_class_broad, cell_doublet, cluster_class, cluster_class_broad, cluster_doublet) 
test4$cells <- str_split_i(test4$cells, "[0-9]_", 2)               
test4 <- left_join(NON_GE@meta.data, test4)
rownames(test4) <- test4$cells
NON_GE@meta.data <- test4

# test4$cell_class_broad %>% table
# test4 %>%
#   dplyr::filter(cluster_class_broad == "non" | cell_class_broad == "non",
#                 cell_doublet == F,
#                 cell_class != "glut",
#                 cell_class != "gaba") %>%
#   dim()
# 
# #Same logic as above
# NON_GE <- subset(NON_GE,
#                   subset = 
#                     ((cluster_class_broad == "non") |
#                        (cell_class_broad == "non"))  &
#                     (cell_doublet == F) &
#                     (cell_class != "gaba")&
#                     (cell_class != "glut"))

NON_GE <- subset(NON_GE,
                  subset = 
                    ((cell_class_broad == "non") |
                       (cluster_class_broad == "non") |
                       (cell_class_broad == "other"))  &
                    #(cell_doublet == F) &
                    (cell_class != "glut") &
                    (cell_class != "gaba"))


lognormal_GE@meta.data$cluster_class_broad %>% table

NON_GE

rm(test,test2,test3,test4, tmp, tmp_markers, tmp_meta, tmp_plot_1, tmp_plot_2)
```

We are left with 16814 cells.  We will filter by low expression again, though.

### Cell Filtering
```{r}
counts <- GetAssayData(object = NON_GE, layer = "counts")

nonzeros <- Matrix::rowSums(counts>2)

nonzeros <- data_frame(gene = rownames(counts),
                       ncells = nonzeros)

nonzeros <- nonzeros %>%
arrange(ncells) %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) 

sum(nonzeros$ncells > 10)

nonzeros %>%
  #dplyr::filter(nonzeros < 10000) %>%
  dplyr::mutate(idx = row_number()/n()*100) %>%
  ggplot() +
  geom_point(aes(y = idx, x = ncells)) +
  theme_JQA() +
    coord_cartesian(xlim=c(0,10000))

nonzeros[c(
  which.max(nonzeros$idx>10),
  which.max(nonzeros$idx>20),
  which.max(nonzeros$idx>30),
  which.max(nonzeros$idx>40),
  which.max(nonzeros$idx>50),
  which.max(nonzeros$idx>60),
  which.max(nonzeros$idx>70),
  which.max(nonzeros$idx>80),
  which.max(nonzeros$idx>90)), ]


#This leaves us with ~ 14k genes which is great  
to_drop <- nonzeros$ncells < 10

sum(!to_drop)
#that leaves 13k genes, which is fine I think

NON_GE <- subset(x = NON_GE,
                      features = nonzeros$gene[!to_drop])
NON_GE 
rm(counts,nonzeros,to_drop)
```

16814 non-neuronal cells, 9907 genes

### Normalization
```{r}
NON_GE <- SplitObject(NON_GE, split.by = "id")
```

```{r}
plan(sequential)
gc()
plan(multisession, workers = 7)
options(future.globals.maxSize = 80 * 1024 ^ 3)
NON_GE <- future.apply::future_lapply(NON_GE,
                                      function(x) {
                                        SCTransform(x,
                                                    assay = "RNA",
                                                    vars.to.regress = "percent.mt",
                                                    verbose = FALSE,
                                                    seed.use = saved_seed,
                                                    )
                                      },
                                      future.seed = saved_seed
) 

NON_GE <- future.apply::future_lapply(NON_GE,
                                            function(x) {
                                              NormalizeData(x,
                                                            assay = "RNA",
                                              )
                                            },
                                            future.seed = saved_seed
) 

NON_GE <- merge(NON_GE[[1]],
              NON_GE[2:14],
              add.cell.ids = names(NON_GE))

NON_GE[["RNA"]] <- JoinLayers(NON_GE[["RNA"]])

VariableFeatures(NON_GE, assay = "SCT") <- rownames(NON_GE[["SCT"]]@scale.data)

# lognormal_GE <- NormalizeData(subset(filtered_GE,
#                                      subset = (id != "HS_3") &
#                                        (id != "HS_16")))
# 
# lognormal_GE <- FindVariableFeatures(lognormal_GE)
# 
# lognormal_GE <- ScaleData(lognormal_GE)
plan(sequential)
gc()
```

```{r}
NON_GE <- RunPCA(NON_GE,
                       assay = "SCT")
```


```{r}
ElbowPlot(NON_GE,ndims = 40)
```

```{r}
DimHeatmap(NON_GE,dims=1:10, cells=500)
DimHeatmap(NON_GE,dims=11:20, cells=500)
DimHeatmap(NON_GE,dims=21:30, cells=500)
DimHeatmap(NON_GE,dims=31:40, cells=500)
```

```{r}
NON_GE <- RunUMAP(NON_GE, 
                        #dims = 1:10,
                        dims = 1:20,
                        reduction = "pca",
                        reduction.name = "umap",
                        seed.use = saved_seed)
```


### Clustering
```{r}
NON_GE <- FindNeighbors(NON_GE,
                              #dims = 1:10,
                              dims = 1:20,
                              reduction = "pca",
                              graph.name = paste0("pca_", c("nn","snn")))
```

```{r}
library(reticulate)
numpy <- reticulate::import("numpy")
pandas <- reticulate::import("pandas")
leidenalg <- reticulate::import("leidenalg")

plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)

NON_GE <- FindClusters(NON_GE,
                             method  = "igraph",
                             algorithm = 4,
                             random.seed = saved_seed,
                             graph.name = "pca_snn",
                             resolution = c(0.01,0.05,0.1,0.25,0.5,0.8, 1, 1.5))

plan(sequential)
gc()
```

```{r}
clustree(NON_GE, prefix = "pca_snn_res.")
```

Another clean clustering, they dont split at all until 0.1

```{r}
DimPlot(NON_GE,
        group.by = c("genotype",
                     "id",
                     "class",
                     "GEX_lane",
                     "pca_snn_res.0.01",
                     "pca_snn_res.0.05",
                     "pca_snn_res.0.1",
                     "pca_snn_res.0.25",
                     "pca_snn_res.0.5",
                     "pca_snn_res.0.8",
                     "pca_snn_res.1",
                     "pca_snn_res.1.5"),
        label = T)
```

There are some genoype/lane effects here driven by individual samples - just like in glut.. though it doesn't look like there will ever be any clusters that really isoalte them..

we need to do 0.5 and recombine, i think, to get all the smaller clusters (which are probably doublets)

```{r}
ggplot(NON_GE@meta.data) +
    geom_bar(aes(x=pca_snn_res.0.5, fill=(cell_class)), position=position_fill()) + theme_JQA()

```
cells are clustering pretty damn well by type


```{r}
NON_GE@meta.data %>% dplyr::count(pca_snn_res.0.5)
```

all the small clusters are below 100 (11-15)

```{r}
DimPlot(NON_GE,
        group.by = c("pca_snn_res.0.5"),
        label = T)
```

```{r}
NON_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.5, fill=(class)), position=position_fill())+ theme_JQA()
#really only two mixed clusters on cell_class

NON_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.5, fill=(id)), position=position_fill())+ theme_JQA()

NON_GE@meta.data %>%
  ggplot() +
  geom_bar(aes(x=pca_snn_res.0.5, fill=(genotype)), position=position_fill())+ theme_JQA()
#12, 8, 3 is domniated by a few samples
```

```{r}
Idents(object = NON_GE) <- "pca_snn_res.0.5"
```


```{r}
test <- NON_GE@meta.data

test2 <- sapply(unique(test$id), function(sample) {
  sapply(as.numeric(levels(test$pca_snn_res.0.5)), function(cluster) {
  fisher.test(test$id == sample,
              test$pca_snn_res.0.5 == cluster,
              alternative = "g")$p.value
  })
})

test3 <- -log10(test2)
test3[test3 == Inf] <- 100
test3[test3 == NA] <- 0
test3[test3 > 25] <- 25
rownames(test3) <- as.numeric(levels(test$pca_snn_res.0.5))

test4 <- data_frame(id = colnames(test3)) %>%
  left_join(sampleInfo)

Heatmap(test3,
        cluster_rows = T,
        row_title_side = "right",
        top_annotation = HeatmapAnnotation(Sex = test4$sex,
                                           Genotype = test4$genotype),
        col = colorRamp2(c(0, 1.3,
                     3.475, 5.650, 7.825, 10, 25), c("black", "white",
                                                 "#00B785","#53CC67","#B2DC3C","#FDE333", "red")))

```


```{r}
FeaturePlot(
  NON_GE,
  reduction = "umap",
  features = c("nFeature_RNA",
               "nCount_RNA",
               "nCount_SCT",
               "nFeature_SCT",
               "percent.hge",
               "percent.mt",
               "percent.rp"),
  order = TRUE,
  pt.size = 0.4,
  min.cutoff = 'q10',
  label = T
)
```


```{r}
NON_GE@meta.data %>%
  dplyr::count(pca_snn_res.0.5)

VlnPlot(NON_GE, features = "nCount_RNA", group.by = "pca_snn_res.0.5", assay = "RNA") /
VlnPlot(NON_GE, features = "nFeature_RNA", group.by = "pca_snn_res.0.5", assay = "RNA")/
VlnPlot(NON_GE, features = "nCount_SCT", group.by = "pca_snn_res.0.5", assay = "RNA") /
VlnPlot(NON_GE, features = "nFeature_SCT", group.by = "pca_snn_res.0.5", assay = "RNA")  
#15 and 21 seem to have double the feature count..

VlnPlot(NON_GE,
        assay = "SCT",
        layer = "data",
        features = c("Gad1", "Gad2", "Slc32a1", "Slc17a7", "Slc17a6"), group.by = "pca_snn_res.0.5", ncol = 2)
```


```{r}
NON_GE_orig <- NON_GE
```

```{r}
plan(sequential)
gc()
plan(multisession, workers = 4)
options(future.globals.maxSize = 80 * 1024 ^ 3)
NON_GE <- PrepSCTFindMarkers(NON_GE)
plan(sequential)
gc()

#1300, 1/4 that (almost) of the gaba cells
```

## Cluster DE Markers
```{r}
Idents(object = NON_GE) <- "pca_snn_res.0.5"

plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)

tmp_markers_non <- future.apply::future_lapply(1:15,
                             function(i)
                             {
                               SeuratWrappers::RunPresto(NON_GE,
                                                         assay = "SCT",
                                                         ident.1 = i)
                             },
                             future.seed = saved_seed
                             )

plan(sequential)
gc()


tmp_markers_non <- lapply(tmp_markers_non, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})

#write.xlsx(tmp_markers_gaba, file = "GABAmarkers_jqa_103124.xlsx")
```

#### Dendogram
```{r}
tmp_tree_non <- BuildClusterTree(NON_GE,
                          assay = "SCT",
                          slot  = "scale.data",
                          dim = 1:20
                          )

PlotClusterTree(tmp_tree_non, 
                direction = "rightwards")
```

```{r}
NON_GE@meta.data$pca_snn_res.0.5 %>% table
```


#### Dotplot markers
```{r}
tmp <- tmp_markers_non %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

#the number of markers we want per cluster
tmp2 <- tmp %>%
  #mutate(cluster = if_else(cluster %in% c("6", "9", "10"), "6", cluster)) %>%
  group_by(cluster) %>%
  slice_max(rank, n = 10) %>%
  #dplyr::filter(rank > 2.8) %>%
  #dplyr::select(-n) %>%
  ungroup()
  
tmp2 %>% count(cluster) %>% arrange(n)

(DotPlot(NON_GE,
        features = unique(c(tmp_neuronal, tmp2$gene)),
        group.by = c("pca_snn_res.0.5"),
        scale = T,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))) %>%
  plotly::ggplotly()


(DotPlot(NON_GE,
        features = unique(c(tmp_neuronal, tmp_non)),
        group.by = c("pca_snn_res.0.5"),
        scale = T,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))) %>%
  plotly::ggplotly()
```
Doublets: 14,12,9,6,7,13
10 needs to be split (it is endothelial plus the bit near 1 is micro)

merge(8,3,4) as oligo

```{r}
test <- FindSubCluster(
  NON_GE,
  cluster = "10",
  graph.name = "pca_snn",
  subcluster.name = "sub.cluster",
  resolution = 0.1,
  algorithm = 4
)

(DotPlot(test,
        features = unique(c(tmp_neuronal, tmp_non)),
        group.by = c("sub.cluster"),
        scale = T,
        cluster.idents = T) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))) %>%
  plotly::ggplotly()

```

not worth it, barely get anyting useable out


```{r}
#Doublets: 14,12,9,6,7,13,10
NON_GE <- subset(NON_GE, subset = (pca_snn_res.0.5 != "14") & 
                   (pca_snn_res.0.5 != "15") &
                   (pca_snn_res.0.5 != "13") &
                   (pca_snn_res.0.5 != "11") &
                   (pca_snn_res.0.5 != "10") &
                   (pca_snn_res.0.5 != "9")&
                   (pca_snn_res.0.5 != "6")&
                   (pca_snn_res.0.5 != "7"))
```

```{r}
tmp <- NON_GE@meta.data

tmp <- tmp %>%
    mutate(pca_snn_res.0.5 = if_else(pca_snn_res.0.5 %in% c("8", "4"), "3", pca_snn_res.0.5)
           )

tmp2 <- tmp %>%
  count(pca_snn_res.0.5) %>%
  arrange(desc(n)) %>%
  mutate(new = factor(row_number()))

tmp2

tmp <- left_join(tmp, tmp2)
tmp$pca_snn_res.0.5 <- tmp$new
tmp <- tmp %>%
  dplyr::select(-new, -n)
rownames(tmp) <- NON_GE[["SCT"]]@data@Dimnames[[2]]

NON_GE@meta.data <- tmp

#Do we need to do this?
Idents(object = NON_GE) <- "pca_snn_res.0.5"
```


Now lets go re-run the cluster thing, we have to set recorrect_umi to false here, since we subset a SCT stack that we previoulsly ran prepsctmarkers on

## Cluster DE Markers 2
```{r}
plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)

tmp_markers_non <- future.apply::future_lapply(1:5,
                             function(i)
                             {
                               SeuratWrappers::RunPresto(NON_GE,
                                                         assay = "SCT",
                                                         ident.1 = i,
                                                         recorrect_umi = F)
                             },
                             future.seed = saved_seed
                             )

plan(sequential)
gc()

tmp_markers_non <- lapply(tmp_markers_non, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})

#write.xlsx(tmp_markers_gaba, file = "GABAmarkers_jqa_103124.xlsx")
```

```{r}
DimPlot(NON_GE,
        group.by = c("pca_snn_res.0.5"),
        label = T)
```

#### Dendogram 2
```{r}
NON_GE <- BuildClusterTree(NON_GE,
                          assay = "SCT",
                          slot  = "data",
                          dim = 1:20
                          )

PlotClusterTree(NON_GE) 

plan(sequential)
gc()
plan(multisession, workers = 6)
options(future.globals.maxSize = 40 * 1024 ^ 3)
tree_markers_non <- future.apply::future_lapply(
  NON_GE@tools$BuildClusterTree[["edge"]][,1] %>% unique(),
  function(i) {
    SeuratWrappers::RunPresto(NON_GE,
                              assay = "SCT",
                              ident.1 = NON_GE@tools$BuildClusterTree,
                              ident.2 = i,
                              recorrect_umi = F)
  },
  future.seed = saved_seed
)
plan(sequential)
gc()

names(tree_markers_non) <- NON_GE@tools$BuildClusterTree[["edge"]][,1] %>% unique()

tree_markers_non <- lapply(tree_markers_non, function(i) {
  i %>% 
 mutate(logFC_rank = percent_rank(avg_log2FC),
         abund = log2((pct.1+0.0001) / (pct.2+0.0001)),
         abund_rank = percent_rank(abund),
         auc_rank = percent_rank(auc),
         neg_logFC_rank = percent_rank(-avg_log2FC),
         neg_abund = log2((pct.2+0.0001) / (pct.1+0.0001)),
         neg_abund_rank = percent_rank(neg_abund),
         rank = if_else((avg_log2FC < 0) | (abund < 0),
                        -neg_abund_rank-neg_logFC_rank,
                        logFC_rank+abund_rank+auc_rank),
         ) %>%
  arrange(desc(rank)) %>%
  dplyr::select(-ends_with("_rank"), -neg_abund) %>%
    rownames_to_column(var = "gene")
})


#  library(ggtree)
# ggtree(tmp_tree@tools$BuildClusterTree,
#           branch.length="none") + 
#      geom_tiplab(as_ylab = F) + 
#      geom_nodelab()+
#      geom_rootedge(rootedge=1) +
#   geom_text(data = data_frame(x = 1:10, y = 1:10, label = 1:10), aes(x=x,y=y,label=label))
# 
# 
#      #geom_cladelab(node = 19, label = "test", align = T, hjust = -10   theme_tree()

```

```{r}
write.xlsx(tmp_markers_non, file = "NONmarkers_jqa_121924.xlsx")
write.xlsx(tree_markers_non, file = "NONtreemarkers_jqa_121924.xlsx")
```

#### Dotplot markers 2


```{r}
tmp_levels <- NON_GE@tools[["BuildClusterTree"]][["edge"]][NON_GE@tools[["BuildClusterTree"]][["edge"]][, 2]<6, 2]

NON_GE@meta.data$pca_snn_res.0.5 <- factor(NON_GE@meta.data$pca_snn_res.0.5,
                                            levels = as.character(tmp_levels))
Idents(NON_GE) <- "pca_snn_res.0.5"
```

We can actually just rename these since we know what they are
```{r}
NON_GE@meta.data$pca_snn_res.0.5 <- factor(NON_GE@meta.data$pca_snn_res.0.5,
                                           levels = as.character(1:5),
                                           labels = c(
                                             "Oligo",
                                             "Micro",
                                             "Astro",
                                             "OPC",
                                             "Immature")
                                           )
Idents(NON_GE) <- "pca_snn_res.0.5"
```



##### Top markers, cluster
```{r}
tmp <- tmp_markers_non %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

tmp3 <- tmp %>%
  mutate(cluster = factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>% 
  group_by(cluster) %>%
  slice_max(rank, n = 5) 

tmp_plot_1 <- wrap_plots(
  
  ggtree::ggtree(NON_GE@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree(),
  
  DotPlot(NON_GE,
          features = (unique(tmp3$gene)),
          group.by = "pca_snn_res.0.5",
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)

tmp_plot_1

tmp_plot_1 <- 
  DotPlot(NON_GE,
          features = (unique(c(tmp3$gene))),
          group.by = "pca_snn_res.0.5",
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed"))


tmp_plot_1

```

##### Top markers, tree
```{r}
tmp <- tree_markers_non %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

#tmp2 <- tmp2[!(tmp2$gene %in% unique(tmp2$gene[duplicated(tmp2$gene)])), ]

tmp3 <- bind_rows(
  tmp %>%
  group_by(cluster) %>%
  slice_max(rank, n = 3),
    tmp %>%
  group_by(cluster) %>%
  slice_min(rank, n = 3)
) %>% arrange(cluster,rank)

tmp_plot_1 <- wrap_plots(
  
  ggtree::ggtree(NON_GE@tools$BuildClusterTree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree(),
  
  DotPlot(NON_GE,
          features = (unique(tmp3$gene)),
          group.by = c("pca_snn_res.0.5"),
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed")),
  nrow =1 ,
  widths = c(0.1,1)
)

tmp_plot_1

```

Cluster 5 is still weird, is it just a doublet?  It has 3 equal populations of gaba/glut/non if we look at cluster class, and it is just a lot of random non-neuronal things if we look at cell class...

1.5 at least splits it into three clusters.. we could also dig into the lognormal meta which has the weightings, i think
```{r}
NON_GE@meta.data %>% dplyr::filter(pca_snn_res.0.5 == "???") %>% count(cluster_class_broad)

NON_GE@meta.data %>% dplyr::filter(pca_snn_res.0.5 == "???") %>% count(cell_class)
NON_GE@meta.data %>% dplyr::filter(pca_snn_res.0.5 == "5") %>% count(pca_snn_res.1.5)

tmp <- lognormal_GE@meta.data %>%
  filter(cells %in% rownames(NON_GE@meta.data %>% dplyr::filter(pca_snn_res.0.5 == "???")))

View(tmp)

tmp %>% count(cluster_class, cluster_doublet)
tmp %>% count(pca_snn_res.0.8)
```

I think 

```{r}
tmp <- NON_GE@meta.data

tmp <- tmp %>%
    mutate(pca_snn_res.0.5 = if_else(pca_snn_res.0.5 %in% c("8", "6", "5"), "5", pca_snn_res.0.5))
    #mutate(pca_snn_res.0.5 = if_else(pca_snn_res.0.5 =="5", pca_snn_res.1.5, pca_snn_res.0.5))

tmp2 <- tmp %>%
  count(pca_snn_res.0.5) %>%
  arrange(desc(n)) %>%
  mutate(new = factor(row_number()))

tmp2

tmp <- left_join(tmp, tmp2)
tmp$pca_snn_res.0.5 <- tmp$new
tmp <- tmp %>%
  dplyr::select(-new, -n)
rownames(tmp) <- NON_GE[["SCT"]]@data@Dimnames[[2]]

NON_GE@meta.data <- tmp

#Do we need to do this?
Idents(object = NON_GE) <- "pca_snn_res.0.5"
```



##### BICCN Neuropiptides
https://www.nature.com/articles/s41586-023-06812-z#Sec10

```{r}
tmp <- tmp_markers_non %>%
  lapply(function(x) {
    x %>% 
      dplyr::filter(!str_starts(gene, "Gm"),
                    !str_starts(gene, "ENS"),
                    !str_ends(gene, "Rik"))
    }
    ) %>%
  bind_rows(.id = "cluster") %>%
  group_by(cluster) %>%
  arrange(cluster, desc(rank)) %>%
  mutate(n = row_number()) %>%
  ungroup()

tmp_features <- c(
  "Adcyap1",	"Npff",
"Adm",	"Nps",
"Agrp",	"Npvf",
"Agt",	"Npw",
"Apln",	"Npy",
"Avp",	"Nts",
"Calca",	"Oxt",
"Calcb",	"Pdyn",
"Cartpt",	"Penk",
"Cck",	"Pmch",
"Cort",	"Pnoc",
"Crh",	"Pomc",
"Edn1",	"Prlh",
"Edn3",	"Pthlh",
"Gal",	"Pyy",
"Gcg",	"Rln3",
"Ghrh",	"Sst",
"Gnrh1",	"Tac1",
"Grp",	"Tac2",
"Hcrt",	"Trh",
"Kiss1",	"Ucn",
"Nmb",	"Ucn3",
"Nms",	"Uts2b",
"Nmu",	"Vip",
"Npb"
)

#tmp_features <- chosen_markers
#tmp_features <- tmp_features[tmp_features %in% rownames(GABA_GE)]

tmp_features <- tmp %>%
  dplyr::filter(gene %in% tmp_features) %>%
  group_by(gene) %>%
  slice_max(rank, n = 1) %>%
  dplyr::mutate(cluster =
                  factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>%
  arrange(cluster) %>%
  pull(gene)


DotPlot(NON_GE,
        features = unique(tmp_features),
        group.by = c("pca_snn_res.0.5"),
        scale = T,
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2))

```

##Comparing to HOchgerner
```{r}

#The vglut1 clusters..
hoch <- list()

tmp <- tmp_markers_glut %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Slc17a7",
                                             "Zic4","Trp73","Rspo2","Sema3e",
                                             "Sema5a","Dcn","Wfs1","Sorc3",
                                             "Cbln1","Coch","Thrsp","Lamp5",
                                             "Bok","Pamr1","Sim1","Ciql1","Gpx3",
                                             "Rxfp1","Ermn","Mpped1"))) %>%
  bind_rows(.id = "cluster") 

hoch[["A"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))

view(tmp)
view(hoch[["A"]])

tmp <-tmp_markers_glut %>%
  lapply(function(x) x %>% dplyr::filter(gene %in%  
                                           c("Slc17a7",
                                             "Cd36","Calb2","Gpr101","Grem1",
                                             "Dcn","C1ql2","Cartpt","Fam46a",
                                             "Fmo1","Rxfp3","Il33","Fbln1","Eps8",
                                             "Cd44", "Fibcd1", "Vit", "Car12",
                                             "Celsr1", "Trh", "Mdga1"
                                            ))) %>%
  bind_rows(.id = "cluster") 

hoch[["B"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))

view(tmp)
view(hoch[["B"]])

tmp <- tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Slc17a7",
                                             "Trh", "Rxfp1", "Medag", "Kit",
                                             "Slc23a", "Plxd3", "Reln",
                                             "Tac1","Igfbp5","Oasl2","Ifit1",
                                             "Grp","Cpne8", "St8sia2", "Cald1",
                                             "Gsg1l", "Prox1", "Igfn1", "Ndst5",
                                             "Mid1","Cdh22"))) %>%
  bind_rows(.id = "cluster")


hoch[["C"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))
view(tmp)
view(hoch[["C"]])

tmp <- tmp_markers_gaba %>%
  lapply(function(x) x %>% dplyr::filter(gene %in% 
                                           c("Lhx6", "Tshz2", "Stab1", "Prl1",
                                             "Nxph1", "Sox6", "Moxd1", "Sst", "Pvalb",
                                             "Lhx8", "Th", "Cbln4", "Luzp2", "Nxph2",
                                             "Prlr", "Greb1", "Calcr", "Tac1", "St18",
                                             "Satb1", "Chodl", "Fign", "Npy", "Tmtc4",
                                             "Nek7", "Rpb4", "Vwc2", "Crabp1",
                                             "Etv1", "Pthlh"))) %>%
  bind_rows(.id = "cluster") 


hoch[["D"]] <- tmp %>%
  group_by(cluster) %>%
  summarize(rank = sum(rank),
            neg_rank = sum(neg_rank))
view(tmp)
view(hoch[["D"]])


tmp_markers_gaba[[1]] %>%
  dplyr::filter(
    gene %in% c("Foxp2", "Fmod",
                "Adra2a",
                "Col6a1",
                "Htr1f",
                "Pax6") 
  )

test <- tmp_markers_gaba %>%lapply(function(x) {
  x %>% dplyr::filter(pct.1 > .75)}
  )

```

##### UMAP
```{r}
test <- NON_GE@reductions$umap@cell.embeddings
test2 <- NON_GE@meta.data

all.equal(rownames(test),
          rownames(test2))

test <- bind_cols(test,test2)

test2 <- test %>%
  group_by(pca_snn_res.0.5) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(cell_color = factor(paste0("Non_",formatC(row_number(),
                                                           format = "d", flag = "0", digits = 1))))

test <- left_join(test, test2 %>% dplyr::select(cell_color, pca_snn_res.0.5))

tmp_plot_2 <- ggplot(test %>%
         group_by(pca_snn_res.0.5) %>%
         slice_sample(prop = 0.25)) +
  geom_point(aes(x = umap_1,
                 y = umap_2,
                 color = cell_color),
             shape = 16) +
  ggrepel::geom_text_repel(data = test %>%
                     group_by(pca_snn_res.0.5) %>%
                     summarize(x = median(umap_1),
                               y = median(umap_2),
                               cell_color = unique(cell_color)),
                   aes(x = x, y = y, label = cell_color),
                  size = 5,
                  min.segment.length = 0.1,
                  max.overlaps = 15,
                  box.padding = 2,
                  color = "grey25") +
  theme_JQA() +
  scale_color_manual(
    values = c(
      hcl.colors(5,palette="Dark3")
    ),
  )  +
  xlab("Umap Dimension 1 (arb)") +
  ylab("Umap Dimension 2 (arb)") +
  theme(legend.position = "none") + ggtitle("Non-Neuronal Neuronal Clusters")

tmp_plot_2

```

```{r}
(tmp_plot_2 + tmp_plot_1) + plot_annotation(tag_levels = "A")
```


```{r}
tmp_file = "Non_seurat_v3.RData"

if(!file.exists(tmp_file)) {
  save(NON_GE, tree_markers_non, tmp_markers_non, file = tmp_file)
} else {
  load(tmp_file)
}

rm(tmp_file)
```

## Supplementary tables

```{r}
NON_GE@meta.data %>% 
  count(pca_snn_res.0.5,sample) %>%
  pivot_wider(names_from = sample,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(pca_snn_res.0.5,
                                 levels = c(1,2,3,4,5,6),
                                 labels = c("Oligo",
                                            "Microglia",
                                            "Astrocyte",
                                            "OPC",
                                            "NA",
                                            "Endo"))) %>%
                  dplyr::select(-pca_snn_res.0.5) %>% view

NON_GE@meta.data %>% 
  count(pca_snn_res.0.5,genotype) %>%
  pivot_wider(names_from = genotype,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(pca_snn_res.0.5,
                                 levels = c(1,2,3,4,5,6),
                                 labels = c("Oligo",
                                            "Microglia",
                                            "Astrocyte",
                                            "OPC",
                                            "NA",
                                            "Endo"))) %>%
                  dplyr::select(-pca_snn_res.0.5) %>% view


NON_GE@meta.data %>% 
  count(pca_snn_res.0.5,sex) %>%
  pivot_wider(names_from = sex,
              values_from = n) %>%
  dplyr::mutate(cluster = factor(pca_snn_res.0.5,
                                 levels = c(1,2,3,4,5,6),
                                 labels = c("Oligo",
                                            "Microglia",
                                            "Astrocyte",
                                            "OPC",
                                            "NA",
                                            "Endo"))) %>%
                  dplyr::select(-pca_snn_res.0.5) %>% view


```


## Compositional Analysis

### Kip's method
Kip just used a glm model for this coupled with the compositions package
```{r}
pacman::p_load(compositions)

##He drops clusters with lower than 30 mean counts per cell, we may do that later

tmp_cluster_counts <- NON_GE@meta.data %>% 
  count(pca_snn_res.0.5,sample) %>%
  pivot_wider(names_from = sample,
              values_from = n) %>%
  dplyr::mutate(cluster = pca_snn_res.0.5) %>%
  dplyr::select(-pca_snn_res.0.5) 

rownames(tmp_cluster_counts) <- tmp_cluster_counts$cluster
clr_cluster_counts <- as.data.frame(compositions::clr(t(tmp_cluster_counts %>% dplyr::select(-cluster))))
colnames(clr_cluster_counts) <- tmp_cluster_counts$cluster

clr_cluster_counts <- merge(sampleInfo, clr_cluster_counts, by.x = "name", by.y = "row.names")
clr_cluster_counts$sex <- as.factor(ifelse(clr_cluster_counts$sex == "F", 1, 0))
clr_cluster_counts$genotype <- as.factor(ifelse(clr_cluster_counts$genotype == "H", 1, 0))

diff_prop_models <- lapply(6:ncol(clr_cluster_counts),
                    function(x){
                      glm(clr_cluster_counts[,x] ~ sex + genotype,
                          data = clr_cluster_counts)
                    })

diff_prop <- lapply(diff_prop_models, summary)
diff_prop_coef <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,1]})))
diff_prop_std <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,2]})))
diff_prop_p <- as.numeric(unlist(lapply(diff_prop, function(x){stats::coef(x)[3,4]})))

diff_prop <- as.data.frame(do.call(cbind,list(colnames(clr_cluster_counts[,-1:-5]),diff_prop_coef,diff_prop_std,diff_prop_p)))

colnames(diff_prop) <- c("ID","Coef","StdErr","pvalue")
diff_prop$pvalue <- as.numeric(as.character(diff_prop$pvalue))
diff_prop_results <- diff_prop[order(diff_prop$pvalue),]
diff_prop_results$FDR <- p.adjust(diff_prop_results$pvalue, method = "fdr")
```
