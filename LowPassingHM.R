LowPassingHM <- function(seurat_object, beta=5, markers = NULL, assay = "RNA",graph.name = "RNA_snn"){
  require("Seurat")
  
  DefaultAssay(seurat_object) <- assay
  expr <- t(GetAssayData(seurat_object))
  graph <- as.matrix(seurat_object@graphs[[graph.name]])[rownames(expr),rownames(expr)]
  L <- diag(rowSums(graph)) - graph
  
  expr <- LapSmooth(L,beta,expr[,markers])
  seurat_object[["hp_exp"]] <- CreateAssayObject(t(expr))
  DefaultAssay(seurat_object) <- "hp_exp"
  seurat_object <- ScaleData(seurat_object)
  DoHeatmap(seurat_object, features = markers) + NoLegend()
}

LapSmooth <- function(L,beta,GEP){
  cellID <- rownames(GEP)
  geneID <- colnames(GEP)
  inv <- solve(diag(nrow(L)) + beta*L)
  GEP <- inv %*% GEP
  rownames(GEP) <- cellID
  colnames(GEP) <- geneID
  rm(inv,geneID,cellID);gc()
  GEP
}