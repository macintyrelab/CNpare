## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,comment = "#>")
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)

## PACKAGES
library(dplyr)
library(CNpare)

# DATA DIRECTORY
path_to_data<-"C:/Users/bhernando/Desktop/CNIO/Projects/BlasProject/quality_tool/data/"

## ----input, include=FALSE-----------------------------------------------------
cells_segcn<-cells_segcn
cells_mapping<-cells_mapping
CCLE_metadata<-CCLE_metadata
CCLE_metadata$Cell.line.primary.name<-toupper(CCLE_metadata$Cell.line.primary.name)

## ----pos_bins-----------------------------------------------------------------
allchr=c(1:22) #Add 23 if want to include chrX
lengthChr<-lengthChr
posBins <- lapply(allchr,function(chr) 
    getBinsStartsEnds(window=500000, chr, lengthChr[chr]))

## ----cnbins-------------------------------------------------------------------
ccle<-as.data.frame(cbind(cellid=cells_mapping$cellid[cells_mapping$study=="CCLE"],sample=cells_mapping$fileid[cells_mapping$study=="CCLE"]))
cells_segcn <- cells_segcn[cells_segcn$sample%in%ccle$sample,]
cells_segcn <- getCINProfiles(segcn=cells_segcn, samples=unique(cells_segcn$sample))
cells_segcn <- left_join(cells_segcn,ccle, by="sample")[,c(1:4,6)]
colnames(cells_segcn)[5]<-"sample"

#get copy-number data per bin
ccle_cn <- getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample))

## ----comparisons, warning=FALSE-----------------------------------------------
exp_cell=as.matrix(ccle_cn[,colnames(ccle_cn)=="COV362"])
colnames(exp_cell)<-"COV362"

measures<-getSimilarities(dat1=exp_cell, dat2=ccle_cn, method="all")

## ----mclosest_cn, echo=FALSE--------------------------------------------------
measures<-measures[order(measures$manhattan),]
head(measures,5)

## ----plot_top1, fig1, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
exp_cell=cells_segcn[cells_segcn$sample=="COV362",]
mod1_cell=cells_segcn[cells_segcn$sample=="KCI-MOH1",]
CNPlot_events(exp_cell,mod1_cell)

## ----pclosest_cn, echo=FALSE--------------------------------------------------
measures<-measures[order(-measures$r),]
head(measures,5)

## ----plot_top2, fig2, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
exp_cell=cells_segcn[cells_segcn$sample=="COV362",]
mod2_cell=cells_segcn[cells_segcn$sample=="HCC1428",]
CNPlot_events(exp_cell,mod2_cell)

## ----plot_differences1, fig3, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
CNPlot_events(exp_cell,mod1_cell,plot_diff = TRUE)

## ----plot_differences2, fig4, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
CNPlot_events(exp_cell,mod2_cell,plot_diff = TRUE)

## ----tissue_origin------------------------------------------------------------
sample="COV362"
cells.tissue<-as.character(na.omit(CCLE_metadata$Cell.line.primary.name[CCLE_metadata$Site.Primary == CCLE_metadata$Site.Primary[CCLE_metadata$Cell.line.primary.name == sample]]))
tissue_segcn<-cells_segcn[cells_segcn$sample %in% cells.tissue,]
tissue.profiles<-getProfiles(segcn=tissue_segcn, samples=unique(tissue_segcn$sample))

## ----quantify_signs-----------------------------------------------------------
component_parameters<-readRDS(paste0(path_to_data,"component_parameters.rds"))

CN_features <- extractCopynumberFeatures(tissue.profiles)
sample_by_component <- generateSampleByComponentMatrix(CN_features, all_components=component_parameters)
signature_quantification <- quantifySignatures(sample_by_component)

## ----plot_k_optimal, fig5, fig.height = 3, fig.width = 6, fig.align = "center", eval=TRUE, echo=FALSE, message=FALSE----
factoextra::fviz_nbclust(t(signature_quantification), kmeans, method = "wss")

## ----plot_clusters, fig6, fig.height = 6, fig.width = 6, fig.align = "center", eval=TRUE, echo=FALSE, message=FALSE----
plotClusters(matrix = t(signature_quantification),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             k=4)

