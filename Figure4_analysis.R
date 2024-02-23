# this script creates feature plot, violin plot and dot plot for bone marrow mouse dataset 
# published in Leimkuehler et al, 2021
library(Seurat)
library(ggplot2)

## datasets

# Leimkuehler et al 2021
mLeimDatPath <- "data/Leimkuehleretal2021/Myelofibrosis scRNA-seq Supplementary Data/Robjects&Markdown/"


fPathLeimkuehler2021d1 <- paste0(mLeimDatPath, "tpo.umap.seurat2.rds")
fPathLeimkuehler2021d2 <- paste0(mLeimDatPath, "jak2.umap.seurat2.rds")
fPathLeimkuehler2021d3 <- paste0(mLeimDatPath, "tpo.gli1.umap.rds")

## results location

mResPath <- "figures/"
rPath1 <- paste0(mResPath, "mouseTHPO/")
rPath2 <- paste0(mResPath, "mouseJAK2/")
rPath3 <- paste0(mResPath, "mouseTHPO_Gli1/")
#bacinPath <- paste0(mResPath, "Baccinetal/")

# Baccin et al 2020
#fPathBaccin2020 <- "C:/CUBA_Projects/BM_stroma/data/Baccinetal/NicheData10x.rda"

## gene lists

goiL1 <- c("Nrp2","Ncam1","Gli1","Cxcl12","Pdgfra","Pdgfrb","Nes",
           "Cd34","Lepr","Cspg4","Cxcr4","Pdpn","Col1a1","Ngfr")

goiL2 <- c("Nrp2","Ncam1","Alpl","Bglap","Col1a1","Ibsp","Runx2","Sp7","Tnfrsf11b")

goiL3 <- c("Nrp1","Nrp2","Ncam1","Sema3c","Sema3g","Sema3f","Vegfa","Vegfb",
           "Vegfc","Figf","Tgfb1","Angptl4","Kdr","Flt1","Flt4","Cxcr4",
           "Tgfbr1","Tgfbr2","Tgfbr3")

goiL4 <- c("Nrp2","Ncam1", "Erg","Hif1a","Zeb1","Mki67","Sox18","Klf4","Rasd1","Rcan1",
           "Angptl2","Esm1", "Cd274")

dotPltList <- unique(c(goiL1,goiL2,goiL3,goiL4))
dotPltList <- dotPltList[c(1:2,21,3:20,22:length(dotPltList))]

## functions

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

plotGenesDp <- function(sObj, dpList, ncolValue,resPath,dataType){
  heightValDP <- 24
  widthValDP <- 30

  dp <- DotPlot(sObj, features = dpList, group.by = "modCellClstName") +
    scale_y_reverse()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.title = element_blank())

  ggsave(dp, filename = paste0(resPath,"DotPlot_all",".pdf"),
         device = "pdf", units = "cm",
         height = heightValDP,width = widthValDP)

  if(dataType == "thpo"){
    ctrl = "EV"
    trt = "TPO"
    sObjSub1 <- subset(sObj,subset = orig.ident==ctrl)
    sObjSub2 <- subset(sObj,subset = orig.ident==trt)

  }else{
    ctrl = "control"
    trt = "exp"
    sObjSub1 <- subset(sObj,subset = rep==ctrl)
    sObjSub2 <- subset(sObj,subset = rep==trt)
  }

  dp1 <- DotPlot(sObjSub1, features = dpList, group.by = "modCellClstName") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.title = element_blank())

  ggsave(dp1, filename = paste0(resPath,"DotPlot_",ctrl,".pdf"),
         device = "pdf", units = "cm",
         height = heightValDP,width = widthValDP)

  dp2 <- DotPlot(sObjSub2, features = dpList, group.by = "modCellClstName") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.title = element_blank())

  ggsave(dp2, filename = paste0(resPath,"DotPlot_",trt,".pdf"),
         device = "pdf", units = "cm",
         height = heightValDP,width = widthValDP)

  stats <- "finished"
  return(stats)
}


plotGenes <- function(sObj, geneList, ncolValue,resPath,dataType){

  heightValVP <- 17
  widthValVP <- 38

  fp <- FeaturePlot(sObj, features = geneList)
  ggsave(fp, filename = paste0(resPath,"FeaturePlot_",paste0(geneList,collapse = "_"),".pdf"),
         device = "pdf", units = "cm",
         height = heightValFP,width = widthValFP)

  vp <- VlnPlot(object = sObj, features = geneList,group.by = "cellClstName",
                ncol = ncolValue) &
    theme(axis.text.x = element_text(size = 18),
          axis.title = element_blank(),
          plot.title = element_text(size = 18))

  ggsave(vp, filename = paste0(resPath,"ViolinPlot_",paste0(geneList,collapse = "_"),".pdf"),
         device = "pdf", units = "in",
         height = heightValVP,width = widthValVP+10)

  if(dataType == "thpo"){
    ctrl = "EV"
    trt = "TPO"
    sObjSub1 <- subset(sObj,subset = orig.ident==ctrl)
    sObjSub2 <- subset(sObj,subset = orig.ident==trt)

  }else{
    ctrl = "control"
    trt = "exp"
    sObjSub1 <- subset(sObj,subset = rep==ctrl)
    sObjSub2 <- subset(sObj,subset = rep==trt)
  }

  vp1 <- VlnPlot(object = sObjSub1, features = geneList,
                 group.by = "cellClstName",
                 ncol = ncolValue) &
    theme(axis.text.x = element_text(size = 18),
          axis.title = element_blank(),
          plot.title = element_text(size = 18))

  ggsave(vp1, filename = paste0(resPath,"ViolinPlot_",ctrl,"_",
                                paste0(geneList,collapse = "_"),".pdf"),
         device = "pdf", units = "in",
         height = heightValVP,width = widthValVP+5)

  vp2 <- VlnPlot(object = sObjSub2, features = geneList,
                 group.by = "cellClstName",
                 ncol = ncolValue) &
    theme(axis.text.x = element_text(size = 18),
          axis.title = element_blank(),
          plot.title = element_text(size = 18))

  ggsave(vp2, filename = paste0(resPath,"ViolinPlot_",trt,"_",
                                paste0(geneList,collapse = "_"),".pdf"),
         device = "pdf", units = "in",
         height = heightValVP,width = widthValVP+5)

  stats <- "finished"
  return(stats)
}

genesDataset <- function(dataPath,ncols,resultLoc,isLiemkuhler,type){
  umap.seurat2 <- loadRData(dataPath)
  if(isLiemkuhler == 1 & type == "thpo"){
    umap.seurat2 <- UpdateSeuratObject(umap.seurat2)
    umap.seurat2@meta.data$name <- factor(umap.seurat2@meta.data$name,
                                          levels = c("EV_early1","EV_late1",
                                                     "EV_early2","EV_late2",
                                                     "TPO_early1","TPO_late1",
                                                     "TPO_early2","TPO_late2"))
    tmpAnno <- gsub("[0-9], ","",umap.seurat2$anno)
    newAnno <- tmpAnno
    newAnno[grep("MSC",newAnno)] <- "MSC"
    newAnno[grep("SCPs",newAnno)] <- "SCPs"
    names(newAnno) <- names(umap.seurat2$name)
    umap.seurat2[["modAnno"]] <- newAnno

    modNme <- gsub("[1-2]","",umap.seurat2$name)
    names(modNme) <- names(umap.seurat2$name)
    umap.seurat2[["modName"]] <- modNme

    combAnno <- paste0(newAnno,"_",umap.seurat2$modName)
    names(combAnno) <- names(umap.seurat2$modName)
    umap.seurat2[["cellClstName"]] <- combAnno

    numVal <- table(umap.seurat2$cellClstName)
    modCombAnno <- paste0(umap.seurat2$cellClstName,
                          " (n=",numVal[as.character(umap.seurat2$cellClstName)],")")
    umap.seurat2[["modCellClstName"]] <- modCombAnno

    cellAnno <- paste0(tmpAnno," ",umap.seurat2$name)
    names(cellAnno) <- names(umap.seurat2$name)
    umap.seurat2[["modCellAnno"]] <- cellAnno
    save(umap.seurat2,file = paste0(rPath1,"seuratTPObj.RData"))
    }
  else if(isLiemkuhler == 1 & type == "thpo_gli1"){
    umap.seurat2 <- UpdateSeuratObject(umap.seurat2)
    umap.seurat2@meta.data$name <- factor(umap.seurat2@meta.data$name,
                                          levels = c("EV_early1","EV_late1",
                                                     "EV_early2","EV_late2",
                                                     "TPO_early1","TPO_late1",
                                                     "TPO_early2","TPO_late2","CK80"))
  }else{
    umap.seurat2 <- UpdateSeuratObject(umap.seurat2)

    tmpAnno <- gsub("[0-9], ","",umap.seurat2$anno)
    newAnno <- tmpAnno
    newAnno[grep("MSC",newAnno)] <- "MSC"
    newAnno[grep("SCPs",newAnno)] <- "SCPs"
    names(newAnno) <- names(umap.seurat2$name)
    umap.seurat2[["modAnno"]] <- newAnno

    combAnno <- paste0(newAnno,"_",umap.seurat2$name)
    names(combAnno) <- names(umap.seurat2$name)
    umap.seurat2[["cellClstName"]] <- combAnno

    cellAnno <- paste0(tmpAnno," ",umap.seurat2$name)
    names(cellAnno) <- names(umap.seurat2$name)
    umap.seurat2[["modCellAnno"]] <- cellAnno

  }

  statsL1 <- plotGenes(sObj = umap.seurat2,geneList = goiL1,
                       ncolValue = ncols,resPath = resultLoc,
                       dataType = type)

  statsL2 <- plotGenes(sObj = umap.seurat2,geneList = goiL2,
                       ncolValue = ncols,resPath = resultLoc,
                       dataType = type)

  statsL3 <- plotGenes(sObj = umap.seurat2,geneList = goiL3,
                       ncolValue = ncols,resPath = resultLoc,
                       dataType = type)

  statsL4 <- plotGenes(sObj = umap.seurat2,geneList = goiL4,
                       ncolValue = ncols,resPath = resultLoc,
                       dataType = type)

  statsL5 <- plotGenesDp(sObj = umap.seurat2,dpList = ?dotPltList,
                          ncolValue = ncols,resPath = resultLoc,
                          dataType = type)
  return(statsL5)
  return(c(statsL1,statsL2,statsL3,statsL4,statsL5))
}


fStats1 <- genesDataset(dataPath = fPathLeimkuehler2021d1, ncols = 5,
                        resultLoc = rPath1, isLiemkuhler = 1, type = "thpo")
fStats2 <- genesDataset(dataPath = fPathLeimkuehler2021d2, ncols = 5,
                        resultLoc = rPath2, isLiemkuhler = 1, type = "jak2")
fStats3 <- genesDataset(dataPath = fPathLeimkuehler2021d3, ncols = 5,
                        resultLoc = rPath3, isLiemkuhler = 1, type = "thpo_gli1")


