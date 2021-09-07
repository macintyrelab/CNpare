###########################################
## FUNCTIONS FOR QUANTIFYING DIFFERENCES ##
###########################################


#' @title Alignment of genomic positions of segments
#' @description This function unifies boundaries of copy-number segments
#' @name unifySegments
#'
#' @param posSeg dataframe with segments of the reference sample
#' @param data dataframe with segments of the experimental model
#'
#' @return dataframe with segments of the model aligned to those in reference
#' @export
#' @examples

unifySegments<-function(posSeg,data){
    #generate output matrix of unified segments of copy number
    cn_new<-c()
    for(j in 1:nrow(posSeg))
    {
        chr<-posSeg[j,1]
        start<-posSeg[j,2]
        end<-posSeg[j,3]
        curr_cn<-data$segVal
        ind<-which(data$chromosome==chr & data$start>=start & data$end<=end)
        if(length(ind)!=0){
            cn_new<-c(cn_new,median(curr_cn[c(ind[1]:ind[length(ind)])],na.rm=T))
        } else {
            ind<-which(data$chromosome==chr & data$start<=start & data$end<=end)
            if(length(ind)!=0){
                cn_new<-c(cn_new,median(curr_cn[c(ind[1]:ind[length(ind)])],na.rm=T))
            } else {
                ind<-which(data$chromosome==chr & data$start>=start & data$end>=end)
                if(length(ind)!=0){
                    cn_new<-c(cn_new,median(curr_cn[c(ind[1]:ind[length(ind)])],na.rm=T))
                } else {
                cn_new<-c(cn_new,NA)
                }
            }
        }
    }
    cn_out<-cbind(posSeg,cn_new)
    colnames(cn_out)<-c("chromosome","start","end","segVal","sample","segVal_B")
    cn_out<-data.frame(cn_out,stringsAsFactors = F)
    #round absolute copy-number
    cn_out$segVal<-round(cn_out$segVal)
    cn_out$segVal_B<-round(cn_out$segVal_B)
    #return
    cn_out
}


#' @title Calculation of the extent of differences between two profiles
#' @description This function computes the % genome difference
#' @name getDifference
#'
#' @param unify dataframe with copy numbers in segments of two profiles.
#' Obtained from unifySegments
#'
#' @return % genome difference between two copy-number profiles
#' @export
#' @examples

getDifference<- function(unify){
    #get genome size
    genome_size<-CNpare:::chr_sizes
    genome_size<-sum(genome_size$length)

    diff<-c()
    #Identify similar and different segments
    sameSegs<-unify[which(unify$segVal==unify$segVal_B),]
    difSegs<-unify[which(unify$segVal!=unify$segVal_B),]

    #Calculate the % of difference
    sum_diff<-sum(as.numeric(difSegs$end-difSegs$start))
    diff<-(sum_diff/genome_size)*100
}
