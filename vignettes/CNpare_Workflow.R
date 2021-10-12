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
exp_cell=as.matrix(ccle_cn[,colnames(ccle_cn)=="OVKATE"])
colnames(exp_cell)<-"OVKATE"

measures<-getSimilarities(dat1=exp_cell, dat2=ccle_cn, method="all")

## ----mclosest_cn, echo=FALSE--------------------------------------------------
measures<-measures[order(measures$manhattan),]
head(measures,5)

## ----plot_top1, fig1, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
exp_cell=cells_segcn[cells_segcn$sample=="OVKATE",]
mod1_cell=cells_segcn[cells_segcn$sample=="Panc 02.03",]
CNPlot_events(exp_cell,mod1_cell,method_diff="non-normalized")

## ----pclosest_cn, echo=FALSE--------------------------------------------------
measures<-measures[order(-measures$r),]
head(measures,5)

## ----plot_probability, fig2, fig.height = 4, fig.width = 4, fig.align = "center", eval=TRUE, echo=FALSE----
plot_simdensity(measures, method="pearson")

## ----plot_differences1, fig3, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
CNPlot_events(exp_cell,mod1_cell, method_diff = "non-normalized", plot_diff = TRUE)

## ----tissue_origin------------------------------------------------------------
sample="OVKATE"
cells.tissue<-as.character(na.omit(CCLE_metadata$Cell.line.primary.name[CCLE_metadata$Site.Primary == CCLE_metadata$Site.Primary[CCLE_metadata$Cell.line.primary.name == sample]]))
tissue_segcn<-cells_segcn[cells_segcn$sample %in% cells.tissue,]
tissue.profiles<-getProfiles(segcn=tissue_segcn, samples=unique(tissue_segcn$sample))

## ----quantify_signs, message=FALSE--------------------------------------------
component_parameters<-readRDS(paste0(path_to_data,"component_parameters.rds"))

CN_features <- extractCopynumberFeatures(tissue.profiles)
sample_by_component <- generateSampleByComponentMatrix(CN_features, all_components=component_parameters)
signature_quantification <- quantifySignatures(sample_by_component)

## ----split_signs--------------------------------------------------------------
exp_cell=as.matrix(signature_quantification[,colnames(signature_quantification)=="OVKATE"])
colnames(exp_cell)<-"OVKATE"
signature_quantification=signature_quantification[,colnames(signature_quantification)!="OVKATE"]

## ----closest_cluster----------------------------------------------------------
samples<-getClusterSamples(matrix=signature_quantification, cell=exp_cell)
print(samples)

## ----plot_clusters, fig4, fig.height = 6, fig.width = 6, fig.align = "center", eval=TRUE, echo=FALSE, message=FALSE----
m<-cbind(signature_quantification,exp_cell)
plotClusters(matrix=m, samples=samples)

