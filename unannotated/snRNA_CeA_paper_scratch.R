---
title: "PARC HP/LP snMultiome paper scratch"
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
# GEO Annotation
JQA 4/1/2025 - this is to make a cell annotation table for the GEO upload, as a processed file.




#Supplemental tables and exports

## Composition
```{r}

tmp <- list(Inh = comp_gaba,
            Exc = comp_glut,
            Non = comp_non)

write.xlsx(tmp, file = "Composition_final.xlsx")

```

## Markers

```{r}
names(markers_gaba) <- paste0("Inh_",formatC(1:18,
                                              format = "d", flag = "0", digits = 1))
names(markers_glut) <- paste0("Exc_",formatC(1:10,
                                              format = "d", flag = "0", digits = 1))
names(markers_non) <- c("Oligo", "Micro", "Astro", "OPC")


write.xlsx(markers_gaba, file = "GABAmarkers_final.xlsx")
write.xlsx(markers_glut, file = "GLUTmarkers_final.xlsx")
write.xlsx(markers_non, file = "NONmarkers_final.xlsx")
```


## DE

```{r}

tmp <- list(Inh = DE_gaba,
            Exc = DE_glut,
            Non = DE_non)

DE_out <- lapply(names(tmp), function(cell_type) {
  
  tmp_DE <- tmp[[cell_type]]
  
  names(tmp_DE) <- c("Pseudo", "Single")
  
  #First we standardize the fdr/p/logfc/mean numbers from each of the different methods
  tmp_DE[["Single"]] <- lapply(tmp_DE[["Single"]],
                               function(x) {
                                 if(!is.null(x[["MAST"]]))
                                   x[["MAST"]] <- x[["MAST"]] %>%
                                     mutate(gene = primerid,
                                            mean = 0,
                                            logfc =0,
                                            p = `Pr(>Chisq)`,
                                            fdr = fdr,
                                            .keep = "none")
                                 if(!is.null(x[["DEsingle"]]))
                                   x[["DEsingle"]] <- x[["DEsingle"]] %>%
                                     rownames_to_column(var = "gene") %>%
                                     mutate(gene = gene,
                                            mean = 0,
                                            logfc =0,
                                            p = `pvalue`,
                                            fdr = pvalue.adj.FDR,
                                            .keep = "none")
                                 return(x)
                               })
  
  tmp_DE[["Pseudo"]] <- lapply(tmp_DE[["Pseudo"]],
                               function(x) {
                                 x[["voom"]] <- NULL
                                 x[["edgeR"]] <- x[["edgeR"]] %>%
                                   rownames_to_column(var = "gene") %>%
                                   mutate(gene = gene,
                                          mean = logCPM,
                                          logfc = logFC,
                                          p = PValue,
                                          fdr = FDR,
                                          .keep = "none")
                                 x[["deseq2"]] <- x[["deseq2"]] %>%
                                   rownames_to_column(var = "gene") %>%
                                   mutate(gene = gene,
                                          mean = baseMean,
                                          logfc = log2FoldChange,
                                          p = pvalue,
                                          fdr = padj,
                                          .keep = "none")
                                 return(x)
                               })
  #Bind all the rows
  tmp_DE <- lapply(tmp_DE, function(x) {
    x[sapply(x, is.null)] <- NULL
    x <- lapply(x, function(x) bind_rows(x, .id = "test_name"))
    return (bind_rows(x, .id = "cluster"))
  }) %>% bind_rows(.id = "type")
  
  tmp_DE <- tmp_DE %>%
    dplyr::select(-type) %>%
    filter(p < 0.05) %>%
    pivot_wider(names_from = c(test_name),
                values_from = c(mean,logfc,p,fdr)) %>%
    dplyr::select(-mean_deseq2, -mean_DEsingle, -mean_MAST, -logfc_deseq2, -logfc_DEsingle, -logfc_MAST) 
  
  #Counts of fdr < 0.05 at the gene level
  tmp_DE$n_DE <- apply((tmp_DE %>% dplyr::select(starts_with("fdr"))) < 0.05,
                 1,
                 function(x) sum(x, na.rm = T))
  tmp_DE <- tmp_DE %>%
    arrange(desc(n_DE)) %>%
    mutate(cluster = paste0(cell_type,"_",formatC(as.numeric(cluster),
                                              format = "d", flag = "0", digits = 1)))
  return(tmp_DE)
})

names(DE_out) <- c("Inh", "Exc", "Non")

DE_out[["Non"]] <- DE_out[["Non"]] %>%
  mutate(cluster = case_match(cluster,
                              "Non_01" ~ "Oligo",
                              "Non_02" ~ "Micro",
                              "Non_03" ~ "Astro",
                              "Non_04" ~ "OPC"
                              ))


                                    
write.xlsx(DE_out, file = "DE_final.xlsx")
```

### Glut plot
```{r}
y2 <- AverageExpression(GLUT_GE,
                         return.seurat =F,
                         assay = "SCT",
                        layer = "counts",
                        group.by= c("sample", "pca_snn_res.0.25"))
y2 <- y2$SCT %>% as.data.frame()

y2 <- AggregateExpression(GLUT_GE,
                         return.seurat =T,
                         assay = "SCT",
                        group.by= c("sample", "pca_snn_res.0.25"))
y2 <-LayerData(y2, assay = "SCT", layer = "data")
y2 <- y2 %>% as.data.frame()


y3 <- y2 %>% t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(name = str_split_i(sample,"_",1),
         cluster = str_split_i(sample,"_",2)) %>%
  left_join(sampleInfo)

y4 <- y3 %>% 
  dplyr::select(name, cluster, genotype, Cck, Cdh8, Nrxn3, Pex5l, Fign, ENSMUSG00000095041) %>%
  pivot_longer(cols = c(Cck, Cdh8, Nrxn3, Pex5l, Fign, ENSMUSG00000095041),
               names_to = "gene",
               values_to = "exp") %>%
  filter(
    (gene %in% c("Cck", "Cdh8", "Nrxn3", "Pex5l") & cluster == 9) |
    (gene == "Fign" & cluster == 8) |
    (gene == "ENSMUSG00000095041" & cluster == 6)
  )
  

tmp_plot <- grouped_ggbetweenstats(y4 %>%
                         mutate(gene = paste0(
                           gene,
                           "\n(",
                           paste0("Exc_",formatC(as.numeric(cluster),
                                              format = "d", flag = "0", digits = 1)),
                           ")")),
                       x = genotype,
                       y = exp,
                       xlab = "Line",
                       ylab = "Psuedo-expression (SCT)",
                       grouping.var = gene,
                       results.subtitle = F,
                       bf.message = F,
                       var.equal = F,
                       annotation.args = list(tag_levels = "A")
                       )

tmp_plot

ggsave("Glut_DE.svg",
       plot = tmp_plot,
       width=9,
       height=5,
       units="in")


```

## Animal Annotatnoi
```{r}
write.xlsx(left_join(sampleInfo, tmp_cuttoffs), file = "Annotation_final.xlsx")
```


# Methods

## Bimodal class counts
```{r}
tmp <- filtered_GE@meta.data

tmp %>%
  group_by(class) %>%
  dplyr::summarize(mean_c = mean(nCount_RNA),
                   sd_c = sd(nCount_RNA),
                   mean_f = mean(nFeature_RNA),
                   sd_f = sd(nFeature_RNA),
                   )

tmp %>%
  dplyr::select(nCount_low, nCount_high) %>%
  distinct() %>%
  dplyr::summarize(mean_low = mean(nCount_low),
                   sd_low = sd(nCount_low),
                   mean_high = mean(nCount_high),
                   sd_high = sd(nCount_high),
                   )

tmp %>%
  dplyr::select(nFeature_low, nFeature_high) %>%
  distinct() %>%
  dplyr::summarize(mean_low = mean(nFeature_low),
                   sd_low = sd(nFeature_low),
                   mean_high = mean(nFeature_high),
                   sd_high = sd(nFeature_high),
                   )

```

#Results

## Coarse Clustering

```{r}
lognormal_GE@meta.data %>%
  filter((cell_doublet == T) |
         (cell_class == "other")) %>%
  dim

 (lognormal_GE@meta.data %>%
  filter((cell_doublet == T) |
         (cell_class == "other")) %>%
  nrow) / nrow(lognormal_GE@meta.data)

tmp_GE@meta.data %>%
  group_by(cluster_class_dotplot) %>%
  dplyr::summarize(mean_c = mean(nCount_RNA),
                   sd_c = sd(nCount_RNA),
                   mean_f = mean(nFeature_RNA),
                   sd_f = sd(nFeature_RNA),
                   )

tmp_GE@meta.data %>%
  dplyr::count(cluster_class_dotplot) %>%
  mutate(p = n/sum(n))

48+23+29
```

# Figures

## Coarse clustering

First we need to do a final filtering

```{r}
tmp_GE <-  subset(lognormal_GE,
                  subset = (cell_doublet == F) &
                    (cell_class != "other"))
tmp_GE@assays$RNA <- NULL
tmp_GE@assays$SCT$counts <- NULL
tmp_GE@assays$SCT$scale.data <- NULL

tmp_GE
```


```{r}
save(tmp_GE, file = "overall_plot_data.RData")
```

```{r}
tmp_gaba <- c("Gad1", "Gad2", "Dlx6os1", "Slc32a1")
tmp_glut <-c("Slc17a7", "Slc17a6", "Nrn1")
tmp_astro <- c("Gja1", "S1pr1", "Fgfr3", "Gli2", "Ntsr2")
tmp_micro <- c("Cx3cr1", "Ikzf1", "C1qa", "C1qb", "C1qc")
tmp_olig <- c("Mal", "Mog", "Mag", "Prr5l", "Plp1")
tmp_OPC <- c("Pdgfra", "Cspg4", "Olig1", "Olig2", "C1ql1")
```


### Dotplot

```{r}
tmp_GE@meta.data$cluster_class_dotplot <- tmp_GE@meta.data$cluster_class_broad
tmp_GE@meta.data$cluster_class_dotplot[tmp_GE@meta.data$doublet] <- "doublet"
tmp_GE@meta.data$cluster_class_dotplot[tmp_GE@meta.data$cluster_class=="imm"] <- "non"

tmp_GE@meta.data$cluster_class_dotplot <- factor(tmp_GE@meta.data$cluster_class_broad,
                                                       levels = c("gaba", "glut","non"))

#table(tmp_GE@meta.data$cell_class)
#Idents(object = lognormal_GE) <- "cell_class_dotplot"
#Idents(object = lognormal_GE) <- "pca_snn_res.0.5"

tmp_plot_2 <-
  DotPlot(tmp_GE,
        features = c(tmp_gaba, tmp_glut, tmp_astro, tmp_micro, tmp_olig, tmp_OPC),
        group.by = "cluster_class_dotplot",
        split.by = "cluster_class_dotplot",
        assay = "SCT",
        cols = c("#0072B4",
                 "#CC1C2F",
                 "#188B41")[c(3,1,2)] ,
        scale = T,
        cluster.idents = F) +
  theme_JQA() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
  ylab("Cluster Cell Type") +
  ggtitle("Marker Gene Expresion") +
  scale_y_discrete(labels = c("Inhibitory\n(17 clusters with\n22,221 cells)",
                              "Excitatory\n(11 clusters with\n10,287 cells)",
                              "Non-Neuronal\n(11 clusters with\n13,731 cells)")
                   ) +
  guides(size = guide_legend(title = "Percent\nExpressed")) +
  theme(axis.title.x = element_blank())
  
  table(tmp_GE@meta.data$cell_class_broad)  

  tmp_GE@meta.data %>%
  dplyr::select(cluster_class_broad, pca_snn_res.0.8) %>%
  distinct() %>% dplyr::count(cluster_class_broad)


tmp_plot_2
```


### UMAP
```{r}
test <- tmp_GE@reductions$umap@cell.embeddings
test2 <- tmp_GE@meta.data

all.equal(rownames(test),
          rownames(test2))

test <- bind_cols(test,test2)

test2 <- test %>%
  mutate(cluster_class_dotplot = factor(cluster_class_broad,
                                        labels = c("Inh", "Exc", "Non"))) %>%
  group_by(pca_snn_res.0.8) %>%
  summarise(n = n(),
            cluster_class_dotplot = unique(cluster_class_dotplot)) %>%
  group_by(cluster_class_dotplot) %>%
  arrange(desc(n)) %>%
  mutate(cluster_color = factor(paste0(str_to_title(cluster_class_dotplot),"_",formatC(row_number(),
                                                           format = "d", flag = "0", digits = 1)))) %>%
  ungroup()

test <- left_join(test, test2 %>% dplyr::select(cluster_color, pca_snn_res.0.8))

test$pca_snn_res.0.8 %>% table

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
      hcl.colors(15,palette="Greens3")[1:11]
    ),
  ) +
  scale_shape_manual(
                     values = c(16,16,16),
                     labels = c("Inhibitory",
                                "Excitatory",
                                "Non-Neuronal")
                     ) +
  guides(color = "none",
         shape = guide_legend("Cluster Cell Type",
                              override.aes = list(
                                size = 5,
                                shape = c(15,15,15),
                                color = c("#0072B4",
                                          "#CC1C2F",
                                          "#188B41")
                                )
                              )
  ) +
  ggtitle("UMAP") +
  xlab("Dim 1 (arb)") +
  ylab("Dim 2 (arb)")

tmp_plot_1
```

### Composition

```{r}
tmp_plot_3 <- test %>%
  mutate(cluster_class_dotplot = factor(cluster_class_dotplot,
                                        labels = c("Inhibitory", "Excitatory", "Non-Neuronal"))) %>%
  group_by(id) %>%
  dplyr::count(cluster_class_dotplot, .drop=F) %>% 
  left_join(sampleInfo) %>% 
  mutate(n = replace_na(n, 0),
         total = sum(n),
         prop = replace_na(n/total,0)) %>%
  ggplot(aes(x = genotype, y = prop, color = genotype)) +
  geom_violin(scale = "width",
              show.legend = F) +
  ggforce::geom_sina(scale = "width",
                     jitter_y = F,
                     size = 2,
                     shape = 16,)+
  scale_color_discrete_qualitative("dark3", labels = c("HP", "LP")) +
  scale_x_discrete(labels = c("HP", "LP") ) +
  theme_JQA()+
  xlab("Line")+
  ylab("Proportion")+
  ggtitle("Overall Line Proportions")+
  guides(color = guide_legend(title = "Line")) +
  expand_limits(y=0) +
  facet_wrap(~cluster_class_dotplot, scales = "free", ncol = 3)

tmp_plot_3
```

### Export
```{r}

tmp_plot <- wrap_plots(
  tmp_plot_1,
  (tmp_plot_2/tmp_plot_3),
  widths = c(3,4),
  guides = "collect"
) + plot_annotation(tag_levels = "A")

tmp_plot

ggsave("Overall.svg",
       plot = tmp_plot,
       width=18,
       height=6,
       units="in")
```



## Inhibitory Neurons

### Composition

```{r}
tmp_plot_1 <- GABA_GE@meta.data %>%
  mutate(pca_snn_res.0.5 = factor(paste0("Inh_", as.numeric(unfactor(pca_snn_res.0.5))),
                levels = paste0("Inh_", 1:18)),
         genotype = factor(genotype,
                           levels = c("H","L"),
                           labels = c("HP", "LP"))) %>%
  group_by(id) %>%
  dplyr::count(pca_snn_res.0.5, .drop=F) %>% 
  left_join(sampleInfo) %>% 
  mutate(n = replace_na(n, 0),
         total = sum(n),
         prop = replace_na(n/total,0)) %>% 
  ggplot(aes(x = genotype, y = prop, color = genotype)) +
  geom_violin(scale = "width",
              show.legend = F) +
  ggforce::geom_sina(scale = "width",
                     jitter_y = F,
                     size = 2,
                     shape = 16,)+
  scale_color_discrete_qualitative("dark3") +
  scale_x_discrete(labels = c("HP", "LP")) +
  theme_JQA()+
  xlab("Line")+
  ylab("Proportion")+
  ggtitle("Inhibitory Neurons Line Proportions")+
  guides(color = guide_legend(title = "Line")) +
  expand_limits(y=0) +
  facet_wrap(~pca_snn_res.0.5, scales = "free", ncol = 6)

tmp_plot_1
```

#### Export
```{r}
ggsave("Inhibitory_composition.svg",
       plot = tmp_plot_1,
       width=10,
       height=5,
       units="in")
```

### UMAP
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
  xlab("Dim 1 (arb)") +
  ylab("Dim 2 (arb)") +
  theme(legend.position = "none") + ggtitle("Inhibitory Neuron Clusters")

tmp_plot_2
```



### Dendro

```{r}
gaba_tree <- GABA_GE@tools$BuildClusterTree
gaba_tree$tip.label <- paste0("Inh_", gaba_tree$tip.label)
gaba_tree <-  (ggtree::ggtree(gaba_tree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree()) %>% 
  ggtree::rotate(19) %>%
  ggtree::rotate(20) %>%
  ggtree::rotate(21)

```


### Dotplot

```{r}
tmp_levels <- GABA_GE@tools[["BuildClusterTree"]][["edge"]][GABA_GE@tools[["BuildClusterTree"]][["edge"]][, 2]<19, 2]

tmp_levels <- GABA_GE@tools[["BuildClusterTree"]]$tip.label[tmp_levels]


GABA_GE@meta.data$pca_snn_res.0.5 <- factor(GABA_GE@meta.data$pca_snn_res.0.5,
                                            levels = as.character(tmp_levels))

Idents(GABA_GE) <- "pca_snn_res.0.5"


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

  "Foxp2", "Prkcd",
  "Ntsr1", "Nts",
  "Tacr1", "Tac1","Isl1",
  "Calcrl", "Adora2a",
  "Drd2", "Drd1",

  "Nos1","Cemip",
  "Vipr2","Vip",
  "Npy1r","Npy",
  "Crh", "Crhr1",
  "Sst",  
  "Rmst", "Glra3", "Grik1",
  #"Slc5a7","Ptk2b","Hapln1", "Cd44", "Chrm2",  "Rai14"
  "Scn4b","Pdyn","Cck","Pnoc","Penk","Nr2f2",
  "Chat", 
  "Wfs1", "Six3", "Meis2","Maf", "Zeb2",
  
  "Oprk1", "Sv2b", "Pde4b", "Chrm3", "Htr4", "Grik3", "Htr2c", "Chrm2", "Scn4b","Tac2", "Grin3a", "Lamp5"
) %>% unique

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
  
  gaba_tree,
  
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

### Export

```{r}
tmp_plot <- wrap_plots(
  tmp_plot_2 ,
  wrap_plots(wrap_elements(tmp_plot_1),
    (gaba_tree +
  ggtree::layout_dendrogram() +
  ggtree::geom_tiplab(angle = 30)),
  ncol = 1,
  heights = c(4,2)),
  widths = c(2,4),
  guides = "collect"
) + plot_annotation(tag_levels = "A")

tmp_plot <- wrap_plots(
  wrap_elements(tmp_plot_1+theme(axis.title.x=element_blank())) ,
  wrap_plots(tmp_plot_2,
    (gaba_tree +
  ggtree::layout_dendrogram() +
  ggtree::geom_tiplab(angle = 30)),
  ncol = 2,
  widths = c(3,4)),
  heights = c(4,3),
  guides = "collect"
) + plot_annotation(tag_levels = "A")

ggsave("Inhibitory.svg",
       plot = tmp_plot,
       width=16,
       height=9,
       units="in")
```



## Excitatory Neurons

### Composition

```{r}
tmp_plot_1 <- GLUT_GE@meta.data %>%
  mutate(pca_snn_res.0.25 = factor(paste0("Exc_", as.numeric(unfactor(pca_snn_res.0.25))),
                levels = paste0("Exc_", 1:10)),
         genotype = factor(genotype,
                           levels = c("H","L"),
                           labels = c("HP", "LP"))) %>%
  group_by(id) %>%
  dplyr::count(pca_snn_res.0.25, .drop=F) %>% 
  left_join(sampleInfo) %>% 
  mutate(n = replace_na(n, 0),
         total = sum(n),
         prop = replace_na(n/total,0)) %>% 
  ggplot(aes(x = genotype, y = prop, color = genotype)) +
  geom_violin(scale = "width",
              show.legend = F) +
  ggforce::geom_sina(scale = "width",
                     jitter_y = F,
                     size = 2,
                     shape = 16)+
  scale_color_discrete_qualitative("dark3", labels = c("HP", "LP")) +
  scale_x_discrete(labels = c("HP", "LP")) +
  theme_JQA()+
  xlab("Line")+
  ylab("Proportion")+
  ggtitle("Excitatory Neurons Line Proportions")+
  guides(color = guide_legend(title = "Line")) +
  expand_limits(y=0) +
  facet_wrap(~pca_snn_res.0.25, scales = "free", ncol = 5)

tmp_plot_1
```

#### Export
```{r}
ggsave("Excitatory_composition.svg",
       plot = tmp_plot_1,
       width=10*5/6,
       height=5*2/3,
       units="in")
```

### UMAP
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
      hcl.colors(10,palette="Dark3")
    ),
  )  +
  xlab("Dim 1 (arb)") +
  ylab("Dim 2 (arb)") +
  theme(legend.position = "none") + ggtitle("Excitatory Neuron Clusters")

tmp_plot_2
```



### Dendro

```{r}
glut_tree <- GLUT_GE@tools$BuildClusterTree
glut_tree$tip.label <- paste0("Exc_", glut_tree$tip.label)
glut_tree <-  (ggtree::ggtree(glut_tree,
         branch.length="none") + 
    ggtree::geom_rootedge(rootedge=1) +
    ggtree::theme_tree())

```


### Dotplot

```{r}
tmp_levels <- GLUT_GE@tools[["BuildClusterTree"]][["edge"]][GLUT_GE@tools[["BuildClusterTree"]][["edge"]][, 2]<11, 2]

tmp_levels <- GLUT_GE@tools[["BuildClusterTree"]]$tip.label[tmp_levels]


GLUT_GE@meta.data$pca_snn_res.0.25 <- factor(GLUT_GE@meta.data$pca_snn_res.0.25,
                                            levels = as.character(tmp_levels))
Idents(GLUT_GE) <- "pca_snn_res.0.25"
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

tmp3 <- tmp %>%
  mutate(cluster = factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>% 
  group_by(cluster) %>%
  slice_max(rank, n = 2) 


tmp_features <- c(
  "Satb2", "Satb1", "Tacr1", "Nts", "Bdnf",  "Oprk1",
                  "Nr2f2", "Rspo2", "Grik4", "Hgf", "Rtn4r")

tmp_features <- tmp %>%
  dplyr::filter(gene %in% unique(c(tmp3$gene, tmp_glut, tmp_features))) %>%
  group_by(gene) %>%
  slice_max(rank, n = 1) %>%
  dplyr::mutate(cluster =
                  factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>%
  arrange(cluster) %>%
  pull(gene)

tmp_plot_1 <- wrap_plots(
  
  glut_tree,
  
  DotPlot(GLUT_GE,
          features = (unique(tmp_features)),
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

### Export

```{r}
tmp_plot <- wrap_plots(
  tmp_plot_2 ,
  wrap_plots(wrap_elements(tmp_plot_1),
    (glut_tree +
  ggtree::layout_dendrogram() +
  ggtree::geom_tiplab(angle = 30)),
  ncol = 1,
  heights = c(4,2)),
  widths = c(2,4),
  guides = "collect"
) + plot_annotation(tag_levels = "A")

tmp_plot <- wrap_plots(
  wrap_elements(tmp_plot_1+theme(axis.title.x=element_blank())) ,
  wrap_plots(tmp_plot_2,
    (glut_tree +
  ggtree::layout_dendrogram() +
  ggtree::geom_tiplab(angle = 30)),
  ncol = 2,
  widths = c(3,4)),
  heights = c(4,3),
  guides = "collect"
) + plot_annotation(tag_levels = "A")

tmp_plot

ggsave("Excitatory.svg",
       plot = tmp_plot,
       width=16,
       height=9,
       units="in")
```


## Non Neurons


### Composition

```{r}
tmp_plot_1 <- NON_GE@meta.data %>%
  mutate(pca_snn_res.0.5 = factor(pca_snn_res.0.5,
                                           levels = as.character(1:4),
                                           labels = c(
                                             "Oligo",
                                             "Micro",
                                             "Astro",
                                             "OPC")),
         genotype = factor(genotype,
                           levels = c("H","L"),
                           labels = c("HP", "LP"))) %>%
  group_by(id) %>%
  dplyr::count(pca_snn_res.0.5, .drop=F) %>% 
  left_join(sampleInfo) %>% 
  mutate(n = replace_na(n, 0),
         total = sum(n),
         prop = replace_na(n/total,0)) %>% 
  ggplot(aes(x = genotype, y = prop, color = genotype)) +
  geom_violin(scale = "width",
              show.legend = F) +
  ggforce::geom_sina(scale = "width",
                     jitter_y = F,
                     size = 2,
                     shape = 16)+
  scale_color_discrete_qualitative("dark3", labels = c("HP", "LP")) +
  scale_x_discrete(labels = c("HP", "LP")) +
  theme_JQA()+
  xlab("Line")+
  ylab("Proportion")+
  ggtitle("Non-Neuronal Line Proportions")+
  guides(color = guide_legend(title = "Line")) +
  expand_limits(y=0) +
  facet_wrap(~pca_snn_res.0.5, scales = "free", ncol = 5)

tmp_plot_1
```

#### Export
```{r}
ggsave("Non_composition.svg",
       plot = tmp_plot_1,
       width=10*4/6,
       height=5/3,
       units="in")
```

### UMAP
```{r}
test <- NON_GE@reductions$umap@cell.embeddings
test2 <- NON_GE@meta.data

all.equal(rownames(test),
          rownames(test2))

test <- bind_cols(test,test2)

test2 <- test %>%
  group_by(pca_snn_res.0.5) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

test2$cell_color = c(
  "Oligo",
  "Micro",
  "Astro",
  "OPC")

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
      hcl.colors(4,palette="Dark3")
    ),
  )  +
  xlab("Dim 1 (arb)") +
  ylab("Dim 2 (arb)") +
  theme(legend.position = "none") + ggtitle("Non-Neuronal Clusters")

tmp_plot_2
```

```{r}
NON_GE@meta.data$pca_snn_res.0.5 <- factor(NON_GE@meta.data$pca_snn_res.0.5,
                                           levels = as.character(1:4),
                                           labels = c(
                                             "Oligo",
                                             "Micro",
                                             "Astro",
                                             "OPC")
                                           )
Idents(NON_GE) <- "pca_snn_res.0.5"
```


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

tmp_features <- c(tmp_OPC, tmp_olig, tmp_astro, tmp_micro)

tmp_features <- tmp %>%
  dplyr::filter(gene %in% unique(c(tmp_features))) %>%
  group_by(gene) %>%
  slice_max(rank, n = 1) %>%
  dplyr::mutate(cluster =
                  factor(cluster,
                         levels = rev(as.character(tmp_levels)))) %>%
  arrange(cluster) %>%
  pull(gene)

tmp_plot_1 <- 
  DotPlot(NON_GE,
          features = (unique(tmp_features)),
          group.by = c("pca_snn_res.0.5"),
          scale = T,
          cluster.idents = F) +
    theme_JQA() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust = 1.2)) +
    ylab("")+
    guides(color = guide_colorbar(title = "Scaled\nExpression"),
           size = guide_legend(title = "Percent\nExpressed"))

tmp_plot_1
```

### Export

```{r}
tmp_plot <- wrap_plots(tmp_plot_2,
                       tmp_plot_1,
                       nrow = 1,
                       guides = "collect",
                       widths = c(3,4)
                       ) + plot_annotation(tag_levels = "A")
tmp_plot

ggsave("Non.svg",
       plot = tmp_plot,
       width=16,
       height=4.5,
       units="in")
```


