## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE,comment = "#>")
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)

## PACKAGES NEEDED
library(dplyr)
library(CNpare)

## ----input, include=FALSE-----------------------------------------------------
cells_segcn<-cells_segcn
cells_mapping<-cells_mapping

## ----pos_bins-----------------------------------------------------------------
allchr=c(1:22) #Add 23 if you want to include chrX
lengthChr=lengthChr
posBins <- lapply(allchr,function(chr) 
    getBinsStartsEnds(window=500000, chr, lengthChr[chr]))

## ----select_samples, echo=FALSE-----------------------------------------------
cells_segcn<-getCINProfiles(segcn=cells_segcn, samples=unique(cells_segcn$sample))
print(paste0("Number of cell lines with detectable CIN: ",length(unique(cells_segcn$sample))))

## ----cnbins-------------------------------------------------------------------
#cell lines with CN profile generated with CCLE and GDSC data
ccle <- cells_mapping[cells_mapping$study=="CCLE",c(1,4)]
ccle <- ccle[ccle$fileid %in% cells_segcn$sample,] #with ASCAT CN profile
gdsc <- cells_mapping[cells_mapping$study=="GDSC",c(1,4)] #in GDSC cellid and fileid is the same
gdsc <- gdsc[gdsc$fileid %in% cells_segcn$sample,] #with ASCAT CN profile
both <- ccle[ccle$cellid %in% gdsc$fileid,] 

# Get copy-number data per bin
ccle_cn <- getCNbins(posBins=posBins, data=cells_segcn, samples=both$fileid)
gdsc_cn <- getCNbins(posBins=posBins, data=cells_segcn, samples=both$cellid)

## ----similarities-------------------------------------------------------------
measures<-getSimilarities(dat1=ccle_cn, dat2=gdsc_cn, method="all")
measures<-left_join(measures,both,by="fileid")

## ----top_lists----------------------------------------------------------------
list.tophits<-getTopHit(samples=both$cellid, measure=measures[,c(7,2:6)], method="all")

## ----allmatch, echo=FALSE-----------------------------------------------------
#Top-1 in all measures --> equal profile with equal ploidy
list.tophits$equal<-list.tophits[,1]==list.tophits[,2] & list.tophits[,1]==list.tophits[,3] & list.tophits[,1]==list.tophits[,4] & list.tophits[,1]==list.tophits[,5]
top1<-list.tophits[list.tophits$equal,] 
print("Number of samples in the Top-1 for all comparison measures") 
print(paste0(nrow(top1),"/304 samples"))

## ----plot_allmatch, fig1, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
b<-both[both$cellid=="A172",]
exp_cell=cells_segcn[cells_segcn$sample==b$cellid,]
mod_cell=cells_segcn[cells_segcn$sample==b$fileid,]
exp_cell$sample<-paste0(b$cellid, "-GDSC")
mod_cell$sample<-paste0(b$cellid, "-CCLE")
CNPlot_events(exp_cell, mod_cell, plot_diff = FALSE)

## ----non-allmatch, echo=FALSE-------------------------------------------------
nontop<-list.tophits[!list.tophits$equal,]
nontop$equal<-nontop[,1]==nontop[,2] & nontop[,1]==nontop[,3]
diffploidy.top1<-nontop[nontop$equal,]
print("Number of samples in the Top-1 for only pearson's r and cosine similarity") 
print(paste0(nrow(diffploidy.top1),"/304 samples"))

## ----plot_non-allmatch, fig2, fig.height = 4, fig.width = 8, fig.align = "center", eval=TRUE, echo=FALSE----
b<-both[both$cellid=="22RV1",]
exp_cell=cells_segcn[cells_segcn$sample==b$cellid,]
mod_cell=cells_segcn[cells_segcn$sample==b$fileid,]
exp_cell$sample<-paste0(b$cellid, "-GDSC")
mod_cell$sample<-paste0(b$cellid, "-CCLE")
CNPlot_events(exp_cell, mod_cell, plot_diff = FALSE)

## ----calculate_allpercentages-------------------------------------------------
differences<-c()
for (i in 1:nrow(both)){
  # Get profile pairs
  gdsc<-both[i,1]
  ccle<-both[i,2]
  
  #Calculate the % of difference
  unify<-unifySegments(cells_segcn[cells_segcn$sample==gdsc,],cells_segcn[cells_segcn$sample==ccle,])
  diff<-getDifference(unify)
  differences<-rbind(differences,c(gdsc,round(diff,2)))
}
colnames(differences)<-c("cellid","percDiff")
differences<-data.frame(differences,stringsAsFactors = F)
differences$percDiff<-as.numeric(differences$percDiff)

## ----plot_percentage_distribution, fig2, fig.height = 4, fig.width = 4, fig.align = "center", eval=TRUE, echo=FALSE----
plot_diffdensity(differences)

