library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

setwd("E:/TRAVAIL/William/Lyon/GSE134722_RAW")

directorylist = c(
  "GSM3964166_FirstInstarLarvalBrainNormalCondition_sample1_10X",
  "GSM3964166_FirstInstarLarvalBrainNormalCondition_sample2_10X",
  "GSM3964166_FirstInstarLarvalBrainNormalCondition_sample3_10X",
  "GSM3964166_FirstInstarLarvalBrainStarvationCondition_sample4_10X",
  "GSM3964166_FirstInstarLarvalBrainStarvationCondition_sample5_10X",
  "GSM3964166_FirstInstarLarvalBrainNormalCondition_sample6_10X",
  "GSM3964166_FirstInstarLarvalBrainStarvationCondition_sample7_10X"
)

for (i in 1:length(directorylist)) {
  path=paste0("./",directorylist[i])
  dir.create(path)
}
listcso=list()
for (i in 1:length(directorylist)) {
  path=paste0("./",directorylist[i])

  from1 = list.files(path)
  to1 = c("barcodes.tsv","genes.tsv","matrix.mtx")
  
  setwd(path)
  for (j in 1:3) {
    file.rename(from1[j],to1[j])
  }
  setwd("E:/TRAVAIL/William/Lyon/GSE134722_RAW")
  
  data <- Read10X(data.dir = path)
  
  name = paste0(ifelse(stringr::str_detect(string=path,pattern="NormalCondition"),"NC","SC"),i)
  print(name)
  listcso[[i]] <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
}

rm(from1,to1,path)

for (i in 1:length(listcso)) {
  listcso[[i]] <- NormalizeData(listcso[[i]], verbose = FALSE)
  listcso[[i]] <- FindVariableFeatures(listcso[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}
#FirstInstarLarvalBrainNormalCondition = NC
#FirstInstarLarvalBrainStarvationCondition = SC
reference = c("NC1", "NC2", "NC3","SC4","SC5","NC6","SC7")
names(listcso) = reference

anchors <- FindIntegrationAnchors(object.list = listcso[reference], dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

integrated <- FindNeighbors(integrated, dims = 1:10)
integrated <- FindClusters(integrated, resolution = 0.5)

markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

p1 <- DimPlot(integrated, reduction = "umap", group.by = "orig.ident")
p2 = FeaturePlot(integrated, features = c("Rfx"))
p1+p2


topmarkers = markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 10, wt = -p_val_adj)
p3 = VlnPlot(integrated, topmarkers$gene, combine=F)



annotation = readxl::read_excel("./elife-50354-fig1-data1-v2.xlsx")
annotation = unique(na.omit(annotation))
annotation = annotation[c(1:10,12:dim(annotation)[1]),]


tmp = full_join(markers,annotation,by="gene")

#load("./snrna_scripts.RData")
save.image("./devinprogress23082020.Rdata")

integrated@meta.data$orig.ident 

