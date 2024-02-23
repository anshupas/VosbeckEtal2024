# this script creates dot plot for bone marrow mouse dataset published in Baccin et al, 2020
library(Seurat)
library(ggplot2)

## datasets

# Baccin et al 2021

mBaccinDatPath <- "data/Baccinetal/NicheData10x.rda"

## results location
mResPath <- "figures/"
bacinPath <- mResPath
#bacinPath <- paste0(mResPath, "Baccinetal/")

# Baccin et al 2020
fPathBaccin2020 <- mBaccinDatPath

## gene lists

goiL1 <- c("Nrp2","Ncam1","Gli1","Cxcl12","Pdgfra","Pdgfrb","Nes",
           "Cd34","Lepr","Cspg4","Cxcr4","Pdpn","Col1a1","Ngfr")

goiL2 <- c("Nrp2","Ncam1","Alpl","Bglap","Col1a1","Ibsp","Runx2","Sp7","Tnfrsf11b")

goiL3 <- c("Nrp1","Nrp2","Ncam1","Sema3c","Sema3g","Sema3f","Vegfa","Vegfb",
           "Vegfc","Figf","Tgfb1","Angptl14","Kdr","Flt1","Flt4","Cxcr4",
           "Tgfbr1","Tgfbr2","Tgfbr3") 

goiL4 <- c("Nrp2","Ncam1", "Erg","Hif1a","Zeb1","Mki67","Sox18","Klf4","Rasd1","Rcan1",
           "Angptl2","Esm1", "Cd274") 
geneListF <- unique(c(goiL1, goiL2, goiL3, goiL4))


## functions

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

plotGenes <- function(sObj, geneList, ncolValue,resPath){
  heightValFP <- 24
  widthValFP <- 26
  
  heightValVP <- 30
  widthValVP <- 55
  
  dp <- DotPlot(sObj, features = geneList) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.title = element_blank())
  
  ggsave(dp, filename = paste0(resPath,"DotPlot_genesOfInterest_v1.pdf"),
         device = "pdf", units = "cm", 
         height = heightValVP,width = widthValVP)
  
  stats <- "finished"
  return(stats)
}

genesDataset <- function(dataPath,ncols,resultLoc,isLiemkuhler,subsetInfo){
  umap.seurat2 <- loadRData(dataPath)
  if(isLiemkuhler == 1){
    umap.seurat2 <- UpdateSeuratObject(umap.seurat2)
  }
  
  else{
    umap.seurat2 <- subset(umap.seurat2, 
                           idents = subsetInfo)
    umap.seurat2@active.ident <- factor(x = umap.seurat2@active.ident, 
                                        levels = subsetInfo)
  }
  
  
  stats <- plotGenes(sObj = umap.seurat2,geneList = geneListF,
                       ncolValue = ncols,resPath = resultLoc)
  
  return(stats)
}

allCells <- c("Ng2+ MSCs", "Adipo-CAR","Osteo-CAR","Osteoblasts", 
              "Fibro/Chondro p.", "Chondrocytes", "Endosteal fibro.", 
              "Arteriolar fibro.", "Stromal fibro.", "Myofibroblasts",
              "Arteriolar ECs", "Sinusoidal ECs", "Mk prog.")

fStats4 <- genesDataset(dataPath = mBaccinDatPath, ncols = 5,
                        resultLoc = bacinPath, isLiemkuhler = 0,
                        subsetInfo = allCells)

subNiche <- subset(NicheData10x,idents = allCells)
subNiche@active.ident <- factor(x = subNiche@active.ident,levels = allCells)


